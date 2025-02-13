      !> Profile normalization
      module normalizer_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Updated to the new and completely
!                             different way of storing the norm and
!                             profile variables (TdPA)
!                           - Changed how to deal with the
!                             normalization data regarding the final
!                             formal solutions (TdPA)
!                           - Added new routines to deal with the
!                             new structures and to the normalization
!                             that is now needed for PRD (TdPA)
!                           - The option to store absorption profiles
!                             in files has been removed (TdPA)
!                           - Removed OpenMP support (TdPA)
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
!  normalization
!    Manage the normalization and precalculation of absorption
!  profiles for the quadrature of the self-consistent problem
!
!  normalize
!    Normalize the absorption profiles for the polarization problem
!
!  normalize_PRD
!    Normalize the absorption profiles for the first order profiles
!  in the PRD emissivity for the polarization problem
!
!  normalizeI:
!    Normalize the absorption profiles for the intensity problem
!
!  normalizeI_PRD
!    Normalize the absorption profiles for the first order profiles
!  in the PRD emissivity for the intensity problem
!
!  writebadbound
!    Write warning about bad normalization of Voigt profiles in a file
!
!  normalize_cle
!    Dummy initialization of the normalization structure for the CLE
!  synthesis mode
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use cram_mod
      use free_mod
      use funnj_mod
      use parameters_mod , only : cZero, TINYB, TINYJS, TINYN, &
                                  sqrt3, PI, c, BADNORM , TINYO, &
                                  B2LK, TINYVEL
      use profile_mod
      use setmpi_mod
      use types_mod

      ! To output the normalization problems
      logical:: obadnorm = .True.

      ! CPU threshold to ask the master for permission
      integer:: cpulimit = 16

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage the normalization and precalculation of absorption
      !! profiles for the quadrature of the self-consistent problem\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!        Bstrength(double(:)): Magnetic field strength\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!          Input(Input_class): Structure with configuration
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!              rlimw(logical): If can write RAM limit message\n
      !!       polarization(logical): Normalizing profiles for
      !!                              polarization problem
      subroutine normalization(Atom,LTElines,Atmo,Bstrength,Geom, &
                               Frec,Red,Input,Flgsg,MPID,rlimw, &
                               polarization)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Input_class), intent(in):: Input
      type(Fctsg_class), intent(inout):: Flgsg
      type(MPI_class), intent(inout):: MPID
      logical, intent(in):: polarization
      logical, intent(inout):: rlimw
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: ofram


      ! Global master, verbose
      if (gpid.eq.0) then
        umsg = ' - Normalizing profiles'
        call verbose
      end if

      ! Polarization
      if (polarization) then

        ! Normalize profiles
        call normalize(Atom,LTElines,Atmo,Bstrength,Geom,MPID, &
                       Frec,Red,Flgsg,Input%folder,ofram,0,0,.False.)
        if (laborted) return

        ! Global master verbose
        if (gpid.eq.0) then
          umsg = ' - Profiles normalized '
          call verbose
        end if

        ! If storing Voigt profiles
        if (VPRAM) then

          ! If CPU went above RAM
          if (ofram.and.rlimw) then

            ! Issue warning
            write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                ' reached the limit of profile allocations.'
            call verbose

            ! No more messages about this
            rlimw = .False.

          end if ! RAM limit reached
        end if ! Storing Voigt profiles

      ! Intensity
      else

        ! Call normalization
        call normalizeI(Atom,LTElines,Atmo,Geom,MPID,Frec,Red, &
                        Input%folder,ofram,0,0,.False.)
        if (laborted) return

        ! Global master verbose
        if (gpid.eq.0) then
          umsg = ' - Profiles normalized '
          call verbose
        end if

        ! If storing Voigt profiles
        if (VIRAM) then

          ! If CPU went above RAM
          if (ofram.and.rlimw) then

            ! Issue warning
            write(umsg,'(A,1x,i4,1x,A)') ' # Processor',pid, &
                  ' reached the limit of profile allocations.'
            call verbose

            ! No more messages about this
            rlimw = .False.

          end if ! RAM limit reached
        end if ! Storing Voigt profiles
      end if ! Intensity

      end subroutine normalization

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize the absorption profiles for the polarization
      !! problem\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!        Bstrength(double(:)): Magnetic field strength\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!        folder(character(:)): Output folder path\n
      !!              ofram(logical): If reached the RAM limit\n
      !!                ith(integer): Polar index if LOS\n
      !!                iph(integer): Azimuth index if LOS\n
      !!                LOS(logical): If normalizing for a LOS
      !!                              direction
      subroutine normalize(Atom,LTElines,Atmo,Bstrength,Geom,MPID, &
                           Frec,Red,Flgsg,folder,ofram,ith,iph,LOS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(inout):: Flgsg
      type(MPI_class), intent(inout):: MPID
      character(len=500), intent(in):: folder
      logical, intent(in):: LOS
      integer, intent(in):: ith,iph
      double precision, dimension(:), intent(in):: Bstrength
      logical, intent(out):: ofram

      ! Local

      logical:: extracomm,field,lVPRAM,lvel
      logical, dimension(:,:,:), allocatable:: tosend

      integer:: ia,jtran,ktran,itermf,itermu,iJf,iJu
      integer:: jdir,iz,jj,njdir,indx
      integer:: ith1,iph1,ifreq,if0,if1
      integer:: indU,indF,indK
      integer:: nMu,nMf,iMu,iMf,iU,mF
      integer:: ios,lcheckram,finished
      integer, dimension(5):: id, id1_b
      integer, dimension(0:nproc-1):: nbf1
      integer, dimension(:), allocatable:: outofbound
      integer, dimension(:,:,:), allocatable:: size1
      integer, dimension(:,:,:), allocatable:: checkram

      double precision:: num,RAM,lRAM,d1,W0,W1,rLu,rLf,S
      double precision:: rJumax,rJfmax,rJu,rJf,rMu,rMf,f62
      double precision:: el,eu,dnubw,au,af,auf,atuf,Dfreqw
      double precision:: DwT,Dw,iDw,vfac,vfacw,ct,st,cc,sc,vel
      double precision, dimension(:),allocatable:: buff1

      complex(kind=8):: prof


      ! Routine name
      urou = 'normalize'

      ! Initialize
      ofram = .False.

      ! Initialize extra communication flag
      if (nproc.gt.cpulimit) then
        extracomm = .True.
      else
        extracomm = .False.
      end if

      ! Get real size of the direction dimension
      if (dyn) then
        njdir = Geom%njdir
      else
        njdir = 1
      end if

      ! Estimate memory in just normalization constants
      call cram_estimate_norm(Atom,LTElines,Atmo,Bstrength, &
                              njdir,.True.,lRAM)

      ! Get current size
      RAM = cram_add(1)

      ! Check if going over limits already
      if ((RAM+lRAM).gt.RLIM) then

        ! We cannot store profiles at all
        lVPRAM = .False.

      ! If there is space
      else

        ! Store if allowed by user
        lVPRAM = VPRAM

      end if

      ! Initialize velocity flag
      lvel = .False.

      !
      ! Normalize active atoms
      !

      ! For each atom
      do ia=1,nA

        ! Allocate and initialize size for MPI
        allocate(size1(Atom(ia)%ntran,Rz0:Rz1,njdir))
        size1 = 0

        ! Allocate and initialize bool for MPI
        allocate(tosend(Atom(ia)%ntran,Rz0:Rz1,njdir))
        ! Initialize tosend
        tosend = .False.

        ! If Master doing MPI
        if (pid.eq.0.and.MPID%mpi) then

          ! Allocate the check for the master
          allocate(checkram(Atom(ia)%ntran,Rz0:Rz1,njdir))

          ! Initialize if storing profiles
          if (lVPRAM) then
            checkram = 1
          else
            checkram = -1
          end if

        end if


        !
        ! Allocate the norm array
        !

        ! For each height
        do iz=Rz0,Rz1

          ! If dynamic
          if (dyn) then

            ! Check local velocity
            vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))
            lvel = vel.gt.TINYVEL

          end if ! Dynamic

          ! Magnetic field?
          field = Bstrength(iz).gt.TINYB

          ! For each direction
          do jdir=1,njdir

            ! Skip static multiple directions
            if (.not.lvel.and.jdir.gt.1) cycle

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Get size depending on the local magnetic field
              if (field) then
                size1(jtran,iz,jdir) = Atom(ia)%trano(jtran)%ncomB
              else
                size1(jtran,iz,jdir) = Atom(ia)%trano(jtran)%ncomNB
              end if

              ! If absent transition, skip
              if (Atom(ia)%fflag(jtran)%absent) cycle

              ! Identify the terms
              itermf = Atom(ia)%fst(jtran)%iterml
              itermu = Atom(ia)%fst(jtran)%itermu

              ! Rolling index
              ktran = jtran + Atom(ia)%tshift

              ! Get index
              indx = Red%idzao(ktran,iz,jdir)

              ! Save in shorter variable
              jj = size1(jtran,iz,jdir)

              ! Allocate and initialize norm
              allocate(Red%dzao(indx)%Norm(jj))
              Red%dzao(indx)%Norm = 0d0

              ! Master doing MPI skips
              if (MPID%mpi.and.pid.eq.0) cycle

              ! If allowed to allocate
              if (lVPRAM) then

                ! Prediction
                ! The subtraction is because the normalization has
                ! already been considered in the count
                d1 = 16d-6*dble((Atom(ia)%if1(jtran) - &
                                 Atom(ia)%if0(jtran) + 1)*jj) - &
                     8d-6*dble(jj)

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored and flag that we are full in RAM
                  Red%dzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%dzao(indx)%VRAM = .True.

                  ! Allocate
                  allocate(Red%dzao(indx)%cp(Atom(ia)%if0(jtran): &
                                             Atom(ia)%if1(jtran),jj))
                  ! Update RAM
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%dzao(indx)%VRAM = .False.

              end if ! Storing

            end do ! transitions
          end do ! directions
        end do ! heights

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
        ! MASTER doing MPI
        !
        if (MPID%mpi.and.pid.eq.0) then

          ! Initialize finished
          finished = 1

          !
          ! Calculate normalization
          !

          ! Until done with all processes
          do while (finished.lt.nproc)

            ! Receive metadata
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

              ! Add to finished and restart
              finished = finished + 1
              cycle

            end if ! Finished CPU

            ! Only in extra mode receive another ping
            if (extracomm) then
              do while (.True.)
                call MPI_SEND(id(1),1,MPI_INTEGER,id(1),id(1), &
                              MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do
            end if

            ! Receive the buffer with the integrals
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

            ! Store checkram if the CPU could not keep this profile
            if (lVPRAM.and.id(5).lt.0) &
              checkram(id(2),id(3),id(4)) = -1

            ! Line index
            ktran = id(2) + Atom(ia)%tshift

            ! Get norm index
            indx = Red%idzao(ktran,id(3),id(4))

            ! Reset index
            jj = 0

            ! Run over components
            do jj=1,size1(id(2),id(3),id(4))

              ! Accumulate the sub-integrals
              Red%dzao(indx)%Norm(jj) = Red%dzao(indx)%Norm(jj) + &
                                        buff1(jj)
            end do ! Components
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
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

            ! Magnetic field?
            field = Bstrength(iz).gt.TINYB

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! Dopple shift init
              vfac = 1d0

              ! If velocity
              if (lvel) then

                ! If emergent
                if (LOS) then

                  ! Trigonometry
                  ct = Geom%L_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph))
                  sc = sin(Geom%L_phi(iph))

                ! If quadrature
                else

                  ! Recover the indexes
                  ith1 = Geom%ithv(jdir)
                  iph1 = Geom%iphv(jdir)

                  ! Trigonometry
                  ct = Geom%V_mu(ith1)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph1)
                  sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

                end if ! LOS or quadrature

                ! Get Doppler shift
                vfac = 1d0 - Atmo%vx(iz)*st*cc - &
                             Atmo%vy(iz)*st*sc - &
                             Atmo%vz(iz)*ct

              end if ! Velocity

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! Skip absent
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! Rolling index
                ktran = jtran + Atom(ia)%tshift

                ! Get index
                indx = Red%idzao(ktran,iz,jdir)

                ! Flag to send
                tosend(jtran,iz,jdir) = .True.

                ! Actual Doppler width
                Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                                Atmo%vmi(iz)**2d0)
                iDw = 1d0/Dw

                ! Get terms
                itermu = Atom(ia)%fst(jtran)%itermu
                itermf = Atom(ia)%fst(jtran)%iterml

                ! Get contributions to damping parameter
                au = Atom(ia)%damp(itermu,iz)
                af = Atom(ia)%damp(itermf,iz)
                auf = Atom(ia)%ldamp(jtran,iz)
                atuf = (au + af + auf)*iDw

                ! Get atomic quantities
                S = Atom(ia)%Sval(itermu)

                ! Upper term
                rLu = Atom(ia)%rLval(itermu)
                rJumax = rLu+S
                nMu = nint(2d0*rJumax+1d0)

                ! Lower term
                rLf = Atom(ia)%rLval(itermf)
                rJfmax = rLf + S
                nMf = nint(2d0*rJfmax+1d0)

                ! Get indexes
                if0 = Atom(ia)%if0(jtran)
                if1 = Atom(ia)%if1(jtran)

                ! Get weights
                W0 = Atom(ia)%W0(jtran)
                W1 = Atom(ia)%W1(jtran)

                ! Scaled Doppler width and normalization
                vfacw = vfac*iDw
                d1 = 1d-5*iDw/sqrt(PI)

                !
                ! Proper normalization
                !

                ! No magnetic field
                if (.not.field) then

                  ! Run over Ju
                  do iJu=1,Atom(ia)%nJ(itermu)

                    ! Get energy, Ju, and index
                    eu = Atom(ia)%FSfreq(iJu,itermu)
                    rJu = Atom(ia)%rJval(iJu,itermu)
                    indU = Atom(ia)%irho(itermu)%irho_ij(iJu)

                    ! Run over Jl
                    do iJf=1,Atom(ia)%nJ(itermf)

                      ! Get energy, Jl, and indexes
                      el = Atom(ia)%FSfreq(iJf,itermf)
                      rJf = Atom(ia)%rJval(iJf,itermf)
                      indF = Atom(ia)%irho(itermf)%irho_ij(iJf)
                      indK = Atom(ia)%trano(jtran)%indNB(indF,indU)

                      ! 6-j
                      f62 = fun6j(rLu,rLf,1d0,rJf,rJu,S,Flgsg)

                      ! Check zero
                      if (abs(f62).lt.TINYJS) cycle

                      !
                      ! Calculate profile
                      !

                      ! Common quantities
                      Dfreqw = (eu - el)*iDw

                      ! Lower boundary

                      ! Voigt
                      call voigt(Dfreqw - Frec%omega(if0)*vfacw, &
                                 atuf,prof)

                      ! Add to norm
                      Red%dzao(indx)%Norm(indK) = dble(prof)*W0*d1

                      ! Store profile
                      if (Red%dzao(indx)%VRAM) &
                        Red%dzao(indx)%cp(if0,indK) = prof

                      ! For each frequency
                      do ifreq=if0+1,if1-1

                        ! Voigt
                        call voigt(Dfreqw - Frec%omega(ifreq)*vfacw, &
                                   atuf,prof)

                        ! Add to the integral
                        Red%dzao(indx)%Norm(indK) = &
                                      Red%dzao(indx)%Norm(indK) + &
                                      dble(prof)*Frec%W_freq(ifreq)*d1

                        ! Store profile
                        if (Red%dzao(indx)%VRAM) &
                          Red%dzao(indx)%cp(ifreq,indK) = prof

                      end do ! frequencies

                      !
                      ! Upper Boundary

                      ! Voigt
                      call voigt(Dfreqw - Frec%omega(if1)*vfacw, &
                                 atuf,prof)

                      ! Add to the integral
                      Red%dzao(indx)%Norm(indK) = &
                                     Red%dzao(indx)%Norm(indK) + &
                                     dble(prof)*W1*d1

                      ! Store profile
                      if (Red%dzao(indx)%VRAM) &
                        Red%dzao(indx)%cp(if1,indK) = prof

                    end do ! Jf
                  end do ! Ju

                ! Yes magnetic field
                else

                  ! Run over Mu
                  do iMu=1,nMu

                    ! Get Mu
                    rMu = -rJumax + dble(iMu-1)

                    ! Run over mu_u
                    do iU=1,Atom(ia)%nblk(iMu,itermu)

                      ! Get energy and index
                      eu = Atom(ia)%eval(iU,iMu,itermu,iz)
                      indU = Atom(ia)%irho(itermu)%jM(iU,iMu)

                      ! Run over Mf
                      do iMf=1,nMf

                        ! Get Mf
                        rMf = -rJfmax + dble(iMf-1)

                        ! Selection rules
                        if (nint(abs(rMu-rMf)).gt.1) cycle

                        ! Run over mu_l
                        do mF=1,Atom(ia)%nblk(iMf,itermf)

                          ! Get energy and index
                          el = Atom(ia)%eval(mF,iMf,itermf,iz)
                          indF = Atom(ia)%irho(itermf)%jM(mF,iMf)
                          indK = Atom(ia)%trano(jtran)%indB(indF,indU)

                          !
                          ! Calculate profile
                          !

                          ! Common quantities
                          Dfreqw = (eu - el + &
                                    Atom(ia)%Dfreq(jtran))*iDw

                          !
                          ! Lower Boundary

                          ! Voigt
                          call voigt(Dfreqw - Frec%omega(if0)*vfacw, &
                                     atuf,prof)

                          ! Add to norm
                          Red%dzao(indx)%Norm(indK) = dble(prof)*W0*d1

                          ! Store profile
                          if (Red%dzao(indx)%VRAM) &
                            Red%dzao(indx)%cp(if0,indK) = prof

                          ! For each frequency
                          do ifreq=if0+1,if1-1

                            ! Voigt
                            call voigt(Dfreqw - &
                                       Frec%omega(ifreq)*vfacw, &
                                       atuf,prof)

                            ! Add to the integral
                            Red%dzao(indx)%Norm(indK) = &
                                      Red%dzao(indx)%Norm(indK) + &
                                      dble(prof)*Frec%W_freq(ifreq)*d1

                            ! Store
                            if (Red%dzao(indx)%VRAM) &
                              Red%dzao(indx)%cp(ifreq,indK) = prof

                          end do ! frequencies

                          !
                          ! Upper boundary

                          ! Voigt
                          call voigt(Dfreqw - Frec%omega(if1)*vfacw, &
                                     atuf,prof)

                          ! Add to norm
                          Red%dzao(indx)%Norm(indK) = &
                                         Red%dzao(indx)%Norm(indK) + &
                                         dble(prof)*W1*d1

                          ! Store
                          if (Red%dzao(indx)%VRAM) &
                            Red%dzao(indx)%cp(if1,indK) = prof

                        end do ! iL
                      end do ! Ml
                    end do ! iU
                  end do ! Mu

                end if ! Magnetic field presence

                ! If doing MPI and there is information to send
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

                  ! Prepare metadata
                  if (Red%dzao(indx)%VRAM) then
                    id1_b = (/ pid, jtran, iz, jdir,  1 /)
                  else
                    id1_b = (/ pid, jtran, iz, jdir, -1 /)
                  end if

                  ! Send metadata
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

                  ! If extra communication, send another ping
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
                      call MPI_SEND(Red%dzao(indx)%Norm(1), &
                                    nbf1(pid), MPI_DOUBLE_PRECISION, &
                                    0, pid, MPI_COMM_RT, ierr)
                      if (ierr.eq.0) exit
                    end do
                  else
                    call MPI_ISEND(Red%dzao(indx)%Norm(1), &
                                   nbf1(pid), MPI_DOUBLE_PRECISION, &
                                   0, pid, MPI_COMM_RT, &
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

            ! Prepare finished metadata
            id = (/ -1, -1, -1, -1, -1 /)

            ! Send metadata
            if (extracomm) then
              do while (.True.)
                call MPI_SEND(id(1),5,MPI_INTEGER,0,0, &
                               MPI_COMM_RT,ierr)
                if (ierr.eq.0) exit
              end do
            else
              call MPI_ISEND(id(1),5,MPI_INTEGER,0,0, &
                              MPI_COMM_RT,MPID%request1,ierr)
            end if ! Type of comm
          end if ! MPI
        end if ! Master/slave

        ! If MPI
        if (MPID%mpi) then

          !
          ! Broadcast the results
          !

          ! Slaves allocating profiles
          if (pid.gt.0.and.lVPRAM) then

            ! Allocate checkram
            allocate(checkram(Atom(ia)%ntran,Rz0:Rz1,njdir))

          end if

          ! For each height
          do iz=Rz0,Rz1

            ! For each direction
            do jdir=1,njdir

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! If no size, skip
                if (size1(jtran,iz,jdir).le.0) cycle

                ! Slave without line
                if (pid.gt.0.and.Atom(ia)%fflag(jtran)%absent) then

                  ! Get dummy
                  call MPI_BCAST(buff1(1), &
                                 size1(jtran,iz,jdir), &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                ! Master or slave with line
                else

                  ! Rolling index
                  ktran = jtran + Atom(ia)%tshift

                  ! Get index
                  indx = Red%idzao(ktran,iz,jdir)

                  ! Share Norm
                  call MPI_BCAST(Red%dzao(indx)%Norm(1), &
                                 size1(jtran,iz,jdir), &
                                 MPI_DOUBLE_PRECISION, 0, &
                                 MPI_COMM_RT, ierr)

                end if ! Absent line

              end do ! transitions
            end do ! directions
          end do ! heights

          ! If storing in RAM
          if (lVPRAM) then

            ! Send checkram
            lcheckram = Atom(ia)%ntran*Rnz*njdir

            ! Share Norm
            call MPI_BCAST(checkram(1,Rz0,1),lcheckram, &
                           MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

            ! Slaves deal with it
            if (pid.gt.0) then

              ! For each direction
              do jdir=1,njdir

                ! For each height
                do iz=Rz0,Rz1

                  ! For each transition
                  do jtran=1,Atom(ia)%ntran

                    ! Skip absent
                    if (Atom(ia)%fflag(jtran)%absent) cycle

                    ! Check size
                    if (size1(jtran,iz,jdir).le.0) cycle

                    ! Rolling index
                    ktran = jtran + Atom(ia)%tshift

                    ! Get index
                    indx = Red%idzao(ktran,iz,jdir)

                    ! If was saving but cannot
                    if (Red%dzao(indx)%VRAM.and. &
                        checkram(jtran,iz,jdir).lt.0) then

                      ! Not storing now
                      Red%dzao(indx)%VRAM = .False.

                      ! Identify the terms
                      itermf = Atom(ia)%fst(jtran)%iterml
                      itermu = Atom(ia)%fst(jtran)%itermu

                      ! Deallocate profile
                      RAM = RAM - 1d-6*dble( &
                                   sizeof(Red%dzao(indx)%cp) + &
                                   sizeof(Red%dzao(indx)%Norm))
                      deallocate(Red%dzao(indx)%cp)

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

        ! Allocate and initialize out of bounds counter
        allocate(outofbound(Atom(ia)%ntran))
        outofbound = 0

        ! If MPI, master does not need this
        if (MPID%mpi.and.pid.eq.0) then

          ! Run over indexes
          do indx=1,Red%ndzao

            ! Deallocate
            deallocate(Red%dzao(indx)%Norm)

          end do

        ! Serial or slave
        else

          ! For each height
          do iz=Rz0,Rz1

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if ! Dynamic

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! Skip absent
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! Rolling index
                ktran = jtran + Atom(ia)%tshift

                ! Get index
                indx = Red%idzao(ktran,iz,jdir)

                ! Skip
                if (indx.le.0) cycle

                ! If storing
                if (Red%dzao(indx)%VRAM) then

                  ! Run over components
                  do jj=1,size(Red%dzao(indx)%Norm)

                    ! Easier to write variable
                    d1 = Red%dzao(indx)%Norm(jj)

                    ! If the norm is not zero
                    if (d1.gt.TINYN) then

                      ! Check close to 1
                      if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then

                        ! Add to count
                        outofbound(jtran) = outofbound(jtran) + 1

                        ! If outputting bad norm info, write it
                        if (obadnorm) &
                          call writebadbound(folder, &
                                             Atom(ia)%Element, &
                                             iz,jdir,.True.,jtran, &
                                             jj,d1)
                      end if

                      ! Normalize profile
                      Red%dzao(indx)%cp(:,jj) = &
                        dcmplx(dble(Red%dzao(indx)%cp(:,jj))/d1, &
                               dimag(Red%dzao(indx)%cp(:,jj)))

                    end if ! Small norm

                  end do ! Components

                  ! Deallocate the norm because it is not needed
                  deallocate(Red%dzao(indx)%Norm)

                ! Not storing
                else

                  ! Run over components
                  do jj=1,size(Red%dzao(indx)%Norm)

                    ! Easier to write variable
                    d1 = Red%dzao(indx)%Norm(jj)

                    ! If the norm is not zero
                    if (d1.gt.TINYN) then

                      ! Check close to 1
                      if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then

                        ! Add to count
                        outofbound(jtran) = outofbound(jtran) + 1

                        ! If outputting bad norm info, write it
                        if (obadnorm) &
                          call writebadbound(folder, &
                                             Atom(ia)%Element, &
                                             iz,jdir,.True.,jtran, &
                                             jj,d1)
                      end if

                      ! Get inverse
                      Red%dzao(indx)%Norm(jj) = 1d0/d1

                    end if ! Small norm

                  end do ! Components

                end if ! Storing

              end do ! transitions
            end do ! directions
          end do ! heights

        end if ! Master and MPI

        ! If MPI
        if (MPID%mpi) then

          ! Add together all bad normalization data

          ! Master
          if (pid.eq.0) then

            ! Master in place
            call MPI_REDUCE(MPI_IN_PLACE,outofbound, &
                            Atom(ia)%ntran,MPI_INTEGER, &
                            MPI_SUM,0,MPI_COMM_RT,ierr)
          ! Slave
          else

            ! Slave just send
            call MPI_REDUCE(outofbound,outofbound, &
                            Atom(ia)%ntran,MPI_INTEGER, &
                            MPI_SUM,0,MPI_COMM_RT,ierr)

          end if ! Master/slave
        end if ! MPI


        ! Master
        if (pid.eq.0) then

          ! Check bad limits
          if (maxval(outofbound).gt.0) then

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Check line is affected
              if (outofbound(jtran).gt.0) then

                ! Write number of  bad normalizations
                write(umsg,'(A,i4,4A,i6,A)') &
                  ' - Warning: transition ',jtran,' in ', &
                  Atom(ia)%Element,' has bad normalization', &
                  ' for the chosen width in ',outofbound(jtran), &
                  ' heights, directions, and components.'
                call verbose

              end if ! If bad normalization

            end do ! Transition

          end if ! Any transition had bad normalization
        end if ! Master

        ! Free
        deallocate(buff1)
        deallocate(size1)
        deallocate(tosend)
        deallocate(outofbound)
        if (allocated(checkram)) deallocate(checkram)

      end do ! Atoms

      !
      ! LTE lines
      !

      ! If LTE lines present
      if (allocated(LTElines)) then

        ! For each LTE line
        do ia=1,size(LTElines)

          ! Skip if absent
          if (LTElines(ia)%absent) cycle

          !
          ! Allocate the norm array
          !

          ! For each height
          do iz=LTElines(ia)%Rz0,Rz1

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if ! Dynamic

            ! Magnetic field?
            field = Bstrength(iz).gt.TINYB

            ! Get size depending on the local magnetic field
            if (field) then
              jj = LTElines(ia)%ncom
            else
              jj = 1
            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! Skip slave
              if (MPID%mpi.and.pid.eq.0) cycle

              ! Get index
              indx = Red%idzao(nxtran+ia,iz,jdir)

              ! Allocate profile itself if storing and it is present
              if (lVPRAM) then

                ! Prediction
                ! No subtraction because the norm itself is not
                ! stored in LTE lines
                d1 = 16d-6*dble(LTElines(ia)%if1 - &
                                LTElines(ia)%if0 + 1)*jj 

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored
                  Red%dzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%dzao(indx)%VRAM = .True.

                  ! Allocate depending on the local magnetic field
                  if (field) then
                    allocate(Red%dzao(indx)%cp(LTElines(ia)%if0: &
                                               LTElines(ia)%if1,jj))
                  else
                    allocate(Red%dzao(indx)%p(LTElines(ia)%if0: &
                                              LTElines(ia)%if1))
                  end if

                  ! Update RAM
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%dzao(indx)%VRAM = .False.

              end if ! Storing

            end do ! directions
          end do ! heights

          !
          ! SLAVE OR SINGLE PROCESSOR
          !
          if (.not.MPID%mpi.or.pid.gt.0) then

            !
            ! Calculate profile
            !

            ! For each height
            do iz=LTElines(ia)%Rz0,Rz1

              ! Thermal part of the Doppler width
              DwT = sqrt(LTElines(ia)%cDopp*LTElines(ia)%cDopp* &
                         Atmo%T(iz) + &
                         Atmo%vmi(iz)*Atmo%vmi(iz))

              ! If dynamic
              if (dyn) then

                ! Check local velocity
                vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                           Atmo%vy(iz)*Atmo%vy(iz) + &
                           Atmo%vz(iz)*Atmo%vz(iz))
                lvel = vel.gt.TINYVEL

              end if ! Dynamic

              ! Magnetic field?
              field = Bstrength(iz).gt.TINYB

              ! Get size depending on the local magnetic field
              if (field) then
                jj = LTElines(ia)%ncom
              else
                jj = 1
              end if

              ! For each direction
              do jdir=1,njdir

                ! Skip static multiple directions
                if (.not.lvel.and.jdir.gt.1) cycle

                ! Dopple shift init
                vfac = 1d0

                ! If velocity
                if (lvel) then

                  ! If emergent
                  if (LOS) then

                    ct = Geom%L_mu(ith)
                    st = sqrt(1d0 - ct*ct)
                    cc = cos(Geom%L_phi(iph))
                    sc = sin(Geom%L_phi(iph))

                  else

                    ! Recover the indexes
                    ith1 = Geom%ithv(jdir)
                    iph1 = Geom%iphv(jdir)
                    ct = Geom%V_mu(ith1)
                    st = sqrt(1d0 - ct*ct)
                    cc = Geom%v_mux(iph1)
                    sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

                  end if ! LOS or quadrature

                  ! Get Doppler shift
                  vfac = 1d0 - Atmo%vx(iz)*st*cc - &
                               Atmo%vy(iz)*st*sc - &
                               Atmo%vz(iz)*ct

                end if ! Velocity

                ! Get index
                indx = Red%idzao(nxtran+ia,iz,jdir)

                ! If not storing, why bother
                if (.not.Red%dzao(indx)%VRAM) cycle

                ! Output Doppler width
                Dw = LTElines(ia)%Dfreq*DwT
                iDw = 1d0/Dw

                ! Get indexes
                if0 = LTElines(ia)%if0
                if1 = LTElines(ia)%if1

                ! Common quantities
                Dfreqw = (LTElines(ia)%eu - LTElines(ia)%el)*iDw
                vfacw = vfac*iDw
                auf = LTElines(ia)%damp(iz)*iDw

                !
                ! No magnetic field
                if (.not.field) then

                  !
                  ! Calculate profile
                  !

                  ! For each frequency
                  do ifreq=if0,if1

                    ! Voigt
                    call voigt(Dfreqw - Frec%omega(ifreq)*vfacw, &
                               auf,prof)

                    ! Save profile
                    Red%dzao(indx)%p(ifreq) = dble(prof)

                  end do ! frequencies

                ! Yes magnetic field
                else

                  ! Initialize
                  jj = 0

                  ! Run over Mu
                  do iMu=1,LTElines(ia)%nMu

                    ! Magnetic number
                    rMu = -LTElines(ia)%Ju + dble(iMu-1)

                    ! Run over Ml
                    do iMf=1,LTElines(ia)%nMl

                      ! Magnetic number
                      rMf = -LTElines(ia)%Jl + dble(iMf-1)

                      ! Selection rules
                      if (nint(abs(rMu-rMf)).gt.1) cycle

                      ! Advance index
                      jj = jj + 1

                      ! Get magnetic shift
                      dnubw = B2LK*Bstrength(iz)* &
                              (LTElines(ia)%gu*rMu - &
                               LTElines(ia)%gl*rMf)*iDw

                      !
                      ! Calculate profile
                      !

                      ! For each frequency
                      do ifreq=if0,if1

                        ! Voigt
                        call voigt(Dfreqw + dnubw - &
                                   Frec%omega(ifreq)*vfacw, &
                                   auf,prof)

                        ! Save profile
                        Red%dzao(indx)%cp(ifreq,jj) = prof

                      end do ! frequencies
                    end do ! Ml
                  end do ! Mu

                end if ! Magnetic field presence

              end do ! output direction
            end do ! height

          end if ! MPI

        end do ! LTE lines

      end if ! LTE lines present

      ! And now update the memory in Normp
      call cram_red_norm(Red,num)
      VRAMc = VRAMc + num
      DRAMc = 0d0

      ! Control
      call control

      return

      end subroutine normalize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize the absorption profiles for the first order
      !! profiles in the PRD emissivity for the polarization problem\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Bstrength(double(:)): Magnetic field strength\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!         ofram(logical): If reached the RAM limit
      subroutine normalize_PRD(Atom,Atmo,Bstrength,Frec,Red,ofram)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      logical, intent(out):: ofram
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      logical:: LVRAM,nfield

      integer:: ia,jtran,itermf,itermu,nMu,nMf,iMu,iMf,iU,mF
      integer:: iz,indx,ifreq,if0,if1,nf,indU,indF,indK

      double precision:: RAM,d1,rMu,rMf,el,eu
      double precision:: rLu,rLf,S,rJumax,rJfmax,rJu,rJf
      double precision:: au,af,auf,atuf,Dfreq,DwT,Dw,iDw
      double precision, dimension(:), allocatable:: W0,W1

      complex(kind=8):: prof


      ! Routine name
      urou = 'normalize_PRD'

      ! Initialize
      ofram = .False.

      ! Allocate space for data
      allocate(Red%pzao(Red%nzao))

      ! Initialize RAM numbers
      RAM = cram_add(1)

      ! For each atom
      do ia=1,nA

        ! Allocate extremal weights
        allocate(W0(Atom(ia)%ntran),W1(Atom(ia)%ntran))

        ! Master get stored weights
        if (pid.eq.0) then
          W0 = Atom(ia)%W0
          W1 = Atom(ia)%W1
        end if

        ! Share limit weights with slaves
        call MPI_BCAST(W0,Atom(ia)%ntran,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_RT,ierr)
        call MPI_BCAST(W1,Atom(ia)%ntran,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_RT,ierr)


        !
        ! Allocate the norm array
        !

        ! For each height
        do iz=Rz0,Rz1_PRD

          ! Thermal part of the Doppler width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

          ! Check if no field
          nfield = Bstrength(iz).le.TINYB

          ! For each transition
          do jtran=1,Atom(ia)%ntran

            ! Skip if no PRD
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! Index
            indx = Red%izao(jtran,ia,iz)

            ! Limits
            if0 = Red%zao(indx)%Igf0
            if1 = Red%zao(indx)%Igf1
            nf = if1 - if0 + 1

            ! Correct weights if not extreme wavelength
            if (if0.ne.Atom(ia)%tif0(jtran)) &
              W0(jtran) = Frec%W_freq(if0)
            if (if1.ne.Atom(ia)%tif1(jtran)) &
              W1(jtran) = Frec%W_freq(if1)

            ! Identify the terms
            itermf = Atom(ia)%fst(jtran)%iterml
            itermu = Atom(ia)%fst(jtran)%itermu

            ! Get atomic quantities
            S = Atom(ia)%Sval(itermu)
            rLu = Atom(ia)%rLval(itermu)
            rLf = Atom(ia)%rLval(itermf)
            rJumax = rLu+S
            nMu = nint(2d0*rJumax+1d0)
            rJfmax = rLf + S
            nMf = nint(2d0*rJfmax+1d0)

            ! Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)
            iDw = 1d0/Dw

            ! Get contributions to damping parameter
            au = Atom(ia)%damp(itermu,iz)
            af = Atom(ia)%damp(itermf,iz)
            auf = Atom(ia)%ldamp(jtran,iz)

            ! No magnetic field
            if (nfield) then

              ! Allocate and initialize
              allocate(Red%pzao(indx)% &
                           Norm(Atom(ia)%trano(jtran)%ncomNB))
              Red%pzao(indx)%Norm = 0d0

              ! Allocate profile itself if storing and it is present
              if (VPRAM.and.Red%zao(indx)%nran.gt.0) then

                ! Prediction
                ! Subtract norm already accounted for
                d1 = 16d-6*dble(Atom(ia)%trano(jtran)%ncomNB*nf) - &
                     8d-6*dble(Atom(ia)%trano(jtran)%ncomNB)

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored
                  Red%pzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%pzao(indx)%VRAM = .True.

                  ! Allocate
                  allocate(Red%pzao(indx)% &
                               cp(if0:if1,Atom(ia)%trano(jtran)% &
                                                   ncomNB))
                  ! Update RAM
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%pzao(indx)%VRAM = .False.

              end if ! Storing

              !
              ! Proper normalization
              !

              ! Common quantities
              d1 = 1d-5*iDw/sqrt(PI)
              atuf = (au+af+auf)*iDw

              ! sum over Ju
              do iU=1,Atom(ia)%nJ(itermu)

                ! Get energy
                eu = Atom(ia)%FSfreq(iU,itermu)

                ! Get indexes
                indU = Atom(ia)%irho(itermu)%irho_ij(iU)

                ! Get Ju
                rJu = Atom(ia)%rJval(iU,itermu)

                ! sum over Jf
                do mF=1,Atom(ia)%nJ(itermf)

                  ! Get Jl
                  rJf = Atom(ia)%rJval(mF,itermf)

                  ! Get energy
                  el = Atom(ia)%FSfreq(mF,itermf)

                  ! Get indexes
                  indF = Atom(ia)%irho(itermf)%irho_ij(mF)
                  indK = Atom(ia)%trano(jtran)%indNB(indF,indU)

                  ! Skip 0
                  if (indK.lt.1) cycle

                  ! Common quantities
                  Dfreq = eu - el

                  ! Boundaries

                  ! Lower

                  ! Voigt
                  call voigt((Dfreq - Frec%omega(if0))*iDw,atuf,prof)

                  ! Add to the integral
                  Red%pzao(indx)%Norm(indK) = dble(prof)*W0(jtran)*d1

                  ! Save profile
                  if (Red%pzao(indx)%VRAM) &
                    Red%pzao(indx)%cp(if0,indK) = prof

                  ! For each frequency
                  do ifreq=if0+1,if1-1

                    ! Voigt
                    call voigt((Dfreq - Frec%omega(ifreq))*iDw, &
                               atuf,prof)

                    ! Add to the integral
                    Red%pzao(indx)%Norm(indK) = &
                                   Red%pzao(indx)%Norm(indK) + &
                                   dble(prof)*(Frec%W_freq(ifreq)*d1)

                    ! Save profile
                    if (Red%pzao(indx)%VRAM) &
                      Red%pzao(indx)%cp(ifreq,indK) = prof

                  end do ! frequencies

                  ! Upper

                  ! Voigt
                  call voigt((Dfreq - Frec%omega(if1))*iDw,atuf,prof)

                  ! Add to the integral
                  Red%pzao(indx)%Norm(indK) = &
                                         Red%pzao(indx)%Norm(indK) + &
                                         dble(prof)*W1(jtran)*d1
                                           
                  ! Save profile
                  if (Red%pzao(indx)%VRAM) &
                    Red%pzao(indx)%cp(if1,indK) = prof

                end do ! Jf
              end do ! Ju

              ! If MPI
              if (nproc.gt.1) then

                ! Compute the norm
                call MPI_ALLREDUCE(MPI_IN_PLACE,Red%pzao(indx)%Norm, &
                                   Atom(ia)%trano(jtran)%ncomNB, &
                                   MPI_DOUBLE_PRECISION,MPI_SUM, &
                                   MPI_COMM_RT,ierr)

                !
                ! Check everyone is saving
                !

                ! Valid range
                if (nf.gt.0) then

                  ! Actual variable
                  LVRAM = Red%pzao(indx)%VRAM

                ! No frequencies
                else

                  ! Send a true
                  LVRAM = .True.

                end if ! Valid range

                ! Check
                call MPI_ALLREDUCE(MPI_IN_PLACE,LVRAM, &
                                   1,MPI_LOGICAL,MPI_LAND, &
                                   MPI_COMM_RT,ierr)

                ! Valid range
                if (nf.gt.0) then

                  ! If was storing, but need to free
                  if (Red%pzao(indx)%VRAM.and..not.LVRAM) then

                    ! Free
                    RAM = RAM - 1d-6*sizeof(Red%pzao(indx)%cp) + &
                          8d-6*dble(Atom(ia)%trano(jtran)%ncomNB)
                    deallocate(Red%pzao(indx)%cp)

                  end if ! Was storing but cannot anymore
                end if ! Valid range

              end if ! MPI

              ! Larger than zero norm
              do indK=1,Atom(ia)%trano(jtran)%ncomNB

                ! If non-zero norm
                if (Red%pzao(indx)%Norm(indK).gt.0d0) &
                  Red%pzao(indx)%Norm(indK) = &
                                      1d0/Red%pzao(indx)%Norm(indK)

                ! If valid and storing
                if (nf.gt.0.and.Red%pzao(indx)%VRAM) then

                  ! Larger than 0 norm, normalize profile
                  if (Red%pzao(indx)%Norm(indK).gt.0d0) &
                    Red%pzao(indx)%cp(:,indK) = &
                             dcmplx(dble(Red%pzao(indx)%cp(:,indK)* &
                                         Red%pzao(indx)%Norm(indK)), &
                                    dimag(Red%pzao(indx)%cp(:,indK)))

                end if

              end do ! Components

            ! Yes magnetic field
            else

              ! Allocate and initialize
              allocate(Red%pzao(indx)%Norm(Atom(ia)%trano(jtran)% &
                                                    ncomB))
              Red%pzao(indx)%Norm = 0d0

              ! Allocate profile itself if storing and it is present
              if (VPRAM.and.Red%zao(indx)%nran.gt.0) then

                ! Prediction
                ! Subtract norm already accounted for
                d1 = 16d-6*dble(Atom(ia)%trano(jtran)%ncomB*nf) - &
                     8d-6*dble(Atom(ia)%trano(jtran)%ncomB)

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored
                  Red%pzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%pzao(indx)%VRAM = .True.

                  ! Allocate
                  allocate(Red%pzao(indx)% &
                               cp(if0:if1,Atom(ia)%trano(jtran)% &
                                                   ncomB))
                  ! Update RAM
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%pzao(indx)%VRAM = .False.

              end if ! Storing

              ! Common quantities
              d1 = 1d-5*iDw/sqrt(PI)
              atuf = (au+af+auf)*iDw

              ! Run over Mu
              do iMu=1,nMu

                ! Get M
                rMu = -rJumax + dble(iMu-1)

                ! Run over mu_u
                do iU=1,Atom(ia)%nblk(iMu,itermu)

                  ! Get energy
                  eu = Atom(ia)%eval(iU,iMu,itermu,iz)

                  ! Get index
                  indU = Atom(ia)%irho(itermu)%jM(iU,iMu)

                  ! Run over Mf
                  do iMf=1,nMf

                    ! Get M
                    rMf = -rJfmax + dble(iMf-1)

                    ! Selection rules
                    if (nint(abs(rMu-rMf)).gt.1) cycle

                    ! Run over mu_f
                    do mF=1,Atom(ia)%nblk(iMf,itermf)

                      ! Get energy
                      el = Atom(ia)%eval(mF,iMf,itermf,iz)

                      ! Get indexes
                      indF = Atom(ia)%irho(itermf)%jM(mF,iMf)
                      indK = Atom(ia)%trano(jtran)%indB(indF,indU)

                      ! Skip 0
                      if (indK.lt.1) cycle

                      ! Common quantities
                      Dfreq = eu - el + Atom(ia)%Dfreq(jtran)

                      ! Boundaries

                      ! Lower

                      ! Voigt
                      call voigt((Dfreq - Frec%omega(if0))*iDw, &
                                 atuf,prof)

                      ! Add to the integral
                      Red%pzao(indx)%Norm(indK) = &
                                              dble(prof)*W0(jtran)*d1

                      ! Save profile
                      if (Red%pzao(indx)%VRAM) &
                        Red%pzao(indx)%cp(if0,indK) = prof

                      ! For each frequency
                      do ifreq=if0+1,if1-1

                        ! Voigt
                        call voigt((Dfreq - Frec%omega(ifreq))*iDw, &
                                   atuf,prof)

                        ! Add to the integral
                        Red%pzao(indx)%Norm(indK) = &
                                      Red%pzao(indx)%Norm(indK) + &
                                      dble(prof)*Frec%W_freq(ifreq)*d1

                        ! Save profile
                        if (Red%pzao(indx)%VRAM) &
                          Red%pzao(indx)%cp(ifreq,indK) = prof

                      end do ! frequencies

                      ! Upper

                      ! Voigt
                      call voigt((Dfreq - Frec%omega(if1))*iDw, &
                                 atuf,prof)

                      ! Add to the integral
                      Red%pzao(indx)%Norm(indK) = &
                                     Red%pzao(indx)%Norm(indK) + &
                                     dble(prof)*W1(jtran)*d1 

                      ! Save profile
                      if (Red%pzao(indx)%VRAM) &
                        Red%pzao(indx)%cp(if1,indK) = prof

                    end do ! if
                  end do ! Mf
                end do ! iU
              end do ! Mu

              ! If MPI
              if (nproc.gt.1) then

                ! Compute the norm
                call MPI_ALLREDUCE(MPI_IN_PLACE,Red%pzao(indx)%Norm, &
                                   Atom(ia)%trano(jtran)%ncomB, &
                                   MPI_DOUBLE_PRECISION,MPI_SUM, &
                                   MPI_COMM_RT,ierr)

                !
                ! Check everyone is saving
                !

                ! Valid range
                if (nf.gt.0) then

                  ! Actual variable
                  LVRAM = Red%pzao(indx)%VRAM

                ! No frequencies
                else

                  ! Send a true
                  LVRAM = .True.

                end if ! Valid range

                ! Check
                call MPI_ALLREDUCE(MPI_IN_PLACE,LVRAM, &
                                   1,MPI_LOGICAL,MPI_LAND, &
                                   MPI_COMM_RT,ierr)
                ! Valid range
                if (nf.gt.0) then

                  ! If was storing, but need to free
                  if (Red%pzao(indx)%VRAM.and..not.LVRAM) then

                    ! Free
                    RAM = RAM - 1d-6*sizeof(Red%pzao(indx)%cp) + &
                          8d-6*dble(Atom(ia)%trano(jtran)%ncomB)
                    deallocate(Red%pzao(indx)%cp)

                  end if ! Was storing but cannot anymore
                end if ! Valid range

              end if ! MPI

              ! For each component
              do indK=1,Atom(ia)%trano(jtran)%ncomB

                ! If non-zero norm, get inverse
                if (Red%pzao(indx)%Norm(indK).gt.0d0) &
                  Red%pzao(indx)%Norm(indK) = &
                                      1d0/Red%pzao(indx)%Norm(indK)

                ! If valid and storing
                if (nf.gt.0.and.Red%pzao(indx)%VRAM) then

                  ! Larger than 0 norm, apply norm
                  if (Red%pzao(indx)%Norm(indK).gt.0d0) &
                    Red%pzao(indx)%cp(:,indK) = &
                             dcmplx(dble(Red%pzao(indx)%cp(:,indK)* &
                                         Red%pzao(indx)%Norm(indK)), &
                                    dimag(Red%pzao(indx)%cp(:,indK)))
                end if

              end do ! Components

            end if ! No magnetic field

          end do ! transitions
        end do ! heights

        ! Free
        deallocate(W0,W1)

      end do ! Atoms

      ! Count RAM
      call cram_red_1stord(Red,RAM)
      ORAMc = RAM
      DRAM2c = 0d0

      ! Control
      call control

      return

      end subroutine normalize_PRD

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize the absorption profiles for the intensity problem\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!        folder(character(:)): Output folder path\n
      !!              ofram(logical): If reached the RAM limit\n
      !!                ith(integer): Polar index if LOS\n
      !!                iph(integer): Azimuth index if LOS\n
      !!                LOS(logical): If normalizing for a LOS
      !!                              direction
      subroutine normalizeI(Atom,LTElines,Atmo,Geom,MPID,Frec,Red, &
                            folder,ofram,ith,iph,LOS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(inout):: MPID
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      character(len=500), intent(in):: folder
      logical, intent(in):: LOS
      integer, intent(in):: ith,iph
      logical, intent(out):: ofram

      ! Local

      logical:: extracomm,lVIRAM,lvel

      integer:: ia,jtran,fjtran,ffjtran,ffktran,indx
      integer:: itermf,itermu,iJf,iJu,jdir,iz,njdir
      integer:: ith1,iph1,ifreq,if0,if1,lcheckram,finished
      integer, dimension(5):: id, id1_b
      integer, dimension(:), allocatable:: outofbound
      integer, dimension(:,:,:), allocatable:: checkram

      double precision:: num,RAM,lRAM,prof,d1,W0,W1,el,eu
      double precision:: au,af,auf,atuf,Dfreqw,buff1
      double precision:: DwT,Dw,iDw,vfac,vfacw,ct,st,cc,sc,vel
      double precision, dimension(1):: ad1


      ! Routine name
      urou = 'normalizeI'

      ! Initialize
      ofram = .False.

      ! Initialize comm flag
      if (nproc.gt.cpulimit) then
        extracomm = .True.
      else
        extracomm = .False.
      end if

      ! Get real size of the direction dimension
      if (dyn) then
        njdir = Geom%njdir
      else
        njdir = 1
      end if

      ! Estimate memory in just normalization constants
      call cram_estimate_norm(Atom,LTElines,Atmo,ad1, &
                              njdir,.False.,lRAM)

      ! Get current size
      RAM = cram_add(1)

      ! Check if going over limits already
      if ((RAM+lRAM).gt.RLIM) then

        ! We cannot store profiles at all
        lVIRAM = .False.

      ! If there is space
      else

        ! Store if allowed by user
        lVIRAM = VIRAM

      end if

      ! Initialize velocity flag
      lvel = .False.

      !
      ! Normalize active atoms
      !

      ! For each atom
      do ia=1,nA

        ! If Master doing MPI
        if (pid.eq.0.and.MPID%mpi) then

          ! Allocate the check for the master
          allocate(checkram(Atom(ia)%ntran,Rz0:Rz1,njdir))

          ! Initialize if storing profiles
          if (VIRAM) then
            checkram = 1
          else
            checkram = -1
          end if

        end if ! Master doing MPI

        !
        ! Allocate the norm array
        !

        ! For each height
        do iz=Rz0,Rz1

          ! Check velocity if dynamic
          if (dyn) then
            vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))
            lvel = vel.gt.TINYVEL
          end if

          ! For each direction
          do jdir=1,njdir

            ! Skip static multiple directions
            if (.not.lvel.and.jdir.gt.1) cycle

            ! For each transition (term)
            do jtran=1,Atom(ia)%ntran

              ! If absent, skip
              if (Atom(ia)%fflag(jtran)%absent) cycle

              ! For each transition (FS)
              do fjtran=1,Atom(ia)%fst(jtran)%nt

                ! Rolling index
                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                ffktran = ffjtran + Atom(ia)%tfshift

                ! Get index
                indx = Red%idzao(ffktran,iz,jdir)

                ! Allocate and initialize
                allocate(Red%dzao(indx)%Norm(1))
                Red%dzao(indx)%Norm = 0d0

                ! Master in MPI skips
                if (MPID%mpi.and.pid.eq.0) cycle

                ! Allocate profile itself if storing
                ! and it is present
                if (lVIRAM) then

                  ! Predict
                  ! Subtract already accounter for norm
                  d1 = 8d-6*dble(Atom(ia)%if1(jtran) - &
                                 Atom(ia)%if0(jtran) + 1) - 8d-6

                  ! If no more space
                  if (floor(RAM+d1).gt.RLIM) then

                    ! No stored
                    Red%dzao(indx)%VRAM = .False.
                    ofram = .True.

                  ! If there is space
                  else

                    ! Storing
                    Red%dzao(indx)%VRAM = .True.

                    ! Allocate
                    allocate(Red%dzao(indx)%p(Atom(ia)%if0(jtran): &
                                              Atom(ia)%if1(jtran)))

                    ! Add normalization to RAM
                    RAM = RAM + d1

                  end if ! Space to store

                ! Not storing Voigt
                else

                  ! No stored
                  Red%dzao(indx)%VRAM = .False.

                end if

              end do ! For each transition (FS)
            end do ! For each transition
          end do ! For each direction
        end do ! For each height

        !
        ! MASTER with MPI
        !
        if (MPID%mpi.and.pid.eq.0) then

          ! Initialize finished
          finished = 1

          !
          ! Calculate normalization
          !

          ! While there is work to do
          do while (finished.lt.nproc)

            ! Receive indexes metadata
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

              ! Add to finished CPUs and restart
              finished = finished + 1
              cycle

            end if

            ! Only in extra moode receive another ping
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
                call MPI_recv(buff1, 1, &
                              MPI_DOUBLE_PRECISION, id(1), &
                              id(1), MPI_COMM_RT, &
                              MPI_STATUS_IGNORE, ierr)
                if (ierr.eq.0) exit
              end do
            else
              call MPI_recv(buff1, 1, &
                            MPI_DOUBLE_PRECISION, id(1), &
                            id(1), MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
            end if

            ! If storing but this CPU could not keep this profile
            if (lVIRAM.and.id(5).lt.0) then

              ! Get index and flag checkram
              jtran = Atom(ia)%ifst(id(2))
              checkram(jtran,id(3),id(4)) = -1

            end if

            ! Rolling index
            ffktran = id(2) + Atom(ia)%tfshift

            ! Get norm index
            indx = Red%idzao(ffktran,id(3),id(4))

            ! Accumulate the sub-integrals
            Red%dzao(indx)%Norm = Red%dzao(indx)%Norm + buff1

          end do ! Communications to do

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
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! Dopple shift init
              vfac = 1d0

              ! If velocity
              if (lvel) then

                ! If emergent
                if (LOS) then

                  ! Trigonometry
                  ct = Geom%L_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph))
                  sc = sin(Geom%L_phi(iph))

                ! If quadrature
                else

                  ! Recover the indexes
                  ith1 = Geom%ithv(jdir)
                  iph1 = Geom%iphv(jdir)

                  ! Trigonometry
                  ct = Geom%V_mu(ith1)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph1)
                  sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

                end if ! Los or quadrature

                ! Get Doppler shift
                vfac = 1d0 - Atmo%vx(iz)*st*cc - &
                             Atmo%vy(iz)*st*sc - &
                             Atmo%vz(iz)*ct

              end if ! Velocity

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! If the line is in this process, skip
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! Get term indexes
                itermu = Atom(ia)%fst(jtran)%itermu
                itermf = Atom(ia)%fst(jtran)%iterml

                ! For each transition (FS)
                do fjtran=1,Atom(ia)%fst(jtran)%nt

                  ! Rolling index
                  ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                  ffktran = ffjtran + Atom(ia)%tfshift

                  ! Get index
                  indx = Red%idzao(ffktran,iz,jdir)

                  ! Get iJ indexes
                  iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
                  iJf = Atom(ia)%fst(jtran)%ilevell(fjtran)

                  ! Get contributions to damping parameter
                  au = Atom(ia)%damp(itermu,iz)
                  af = Atom(ia)%damp(itermf,iz)
                  auf = Atom(ia)%ldamp(jtran,iz)

                  ! Get indexes
                  if0 = Atom(ia)%if0(jtran)
                  if1 = Atom(ia)%if1(jtran)

                  ! Get Weights
                  W0 = Atom(ia)%W0(jtran)
                  W1 = Atom(ia)%W1(jtran)

                  !
                  ! Proper normalization
                  !

                  ! Get the FS energies
                  eu = Atom(ia)%FSfreq(iJu,itermu)
                  el = Atom(ia)%FSfreq(iJf,itermf)

                  ! Output Doppler width
                  Dw = (eu - el)*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)
                  iDw = 1d0/Dw

                  ! Scaled Doppler width and normalization
                  vfacw = vfac*iDw
                  d1 = 1d-5*iDw/sqrt(PI)

                  ! Total damping
                  atuf = (au + af + auf)*iDw

                  !
                  ! Calculate profile
                  !

                  ! Common quantities
                  Dfreqw = (eu - el)/Dw

                  ! Lower boundary

                  ! Voigt
                  call voigtI(Dfreqw - Frec%omega(if0)*vfacw, &
                             atuf,prof)

                  ! Add to norm
                  Red%dzao(indx)%Norm = prof*W0*d1

                  ! Store profile
                  if (Red%dzao(indx)%VRAM) &
                    Red%dzao(indx)%p(if0) = prof

                  ! For each frequency
                  do ifreq=if0+1,if1-1

                    ! Voigt
                    call voigtI(Dfreqw - Frec%omega(ifreq)*vfacw, &
                               atuf,prof)

                    ! Add to the integral
                    Red%dzao(indx)%Norm = Red%dzao(indx)%Norm + &
                                          prof*Frec%W_freq(ifreq)*d1

                    ! Store profile
                    if (Red%dzao(indx)%VRAM) &
                      Red%dzao(indx)%p(ifreq) = prof

                  end do ! Frequencies

                  ! Upper boundary

                  ! Voigt
                  call voigtI(Dfreqw - Frec%omega(if1)*vfacw, &
                              atuf,prof)

                  ! Add to the integral
                  Red%dzao(indx)%Norm = Red%dzao(indx)%Norm + &
                                        prof*W1*d1

                  ! Store profile
                  if (Red%dzao(indx)%VRAM) &
                    Red%dzao(indx)%p(if1) = prof

                  ! If MPI
                  if (MPID%mpi) then

                    ! Check last send was received
                    if (.not.extracomm) then
                      call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE, &
                                    ierr)
                      call MPI_WAIT(MPID%request2,MPI_STATUS_IGNORE, &
                                    ierr)
                    end if

                    ! Prepare metadata
                    if (Red%dzao(indx)%VRAM) then
                      id1_b = (/ pid, ffjtran, iz, jdir, 1 /)
                    else
                      id1_b = (/ pid, ffjtran, iz, jdir,-1 /)
                    end if

                    ! Send metadata
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
                    buff1 = Red%dzao(indx)%Norm(1)

                    ! If extra communication, send ping
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
                        call MPI_SEND(buff1,1, &
                                      MPI_DOUBLE_PRECISION, 0, pid, &
                                      MPI_COMM_RT, ierr)
                        if (ierr.eq.0) exit
                      end do
                    else
                      call MPI_ISEND(buff1,1, &
                                     MPI_DOUBLE_PRECISION, 0, pid, &
                                     MPI_COMM_RT, &
                                     MPID%request2, ierr)
                    end if
                  end if ! MPI

                end do ! Fine transition
              end do ! output transition
            end do ! output direction
          end do ! height

          ! If MPI send finished signal
          if (MPID%mpi) then

            ! Check last send was received
            if (.not.extracomm) &
              call MPI_WAIT(MPID%request1,MPI_STATUS_IGNORE,ierr)

            ! Prepare finished metadata
            id = (/ -1, -1, -1, -1, -1 /)

            ! Send metadata
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

          ! Slaves allocating profiles
          if (pid.gt.0.and.VIRAM) then

            ! Allocate checkram
            allocate(checkram(Atom(ia)%ntran,Rz0:Rz1,njdir))

          end if

          ! For each height
          do iz=Rz0,Rz1

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! For each transition (FS)
                do fjtran=1,Atom(ia)%fst(jtran)%nt

                  ! Rolling index
                  ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

                  ! Slave without line
                  if (pid.gt.0.and.Atom(ia)%fflag(jtran)%absent) then

                    ! Get dummy
                    call MPI_BCAST(buff1, 1, &
                                   MPI_DOUBLE_PRECISION, 0, &
                                   MPI_COMM_RT, ierr)

                  ! Master or slave with line
                  else

                    ! Rolling index
                    ffktran = ffjtran + Atom(ia)%tfshift

                    ! Get index
                    indx = Red%idzao(ffktran,iz,jdir)

                    ! Share Norm
                    call MPI_BCAST(Red%dzao(indx)%Norm, 1, &
                                   MPI_DOUBLE_PRECISION, 0, &
                                   MPI_COMM_RT, ierr)

                  end if

                end do ! transitions (FS)
              end do ! transitions
            end do ! directions
          end do ! heights

          ! If storing in RAM
          if (lVIRAM) then

            ! Send checkram
            lcheckram = Atom(ia)%ntran*Rnz*njdir

            ! Share Norm
            call MPI_BCAST(checkram(1,Rz0,1),lcheckram, &
                           MPI_INTEGER, 0, &
                           MPI_COMM_RT, ierr)

            ! Slaves deal with it
            if (pid.gt.0) then

              ! For each height
              do iz=Rz0,Rz1

                ! If dynamic
                if (dyn) then

                  ! Check local velocity
                  vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                             Atmo%vy(iz)*Atmo%vy(iz) + &
                             Atmo%vz(iz)*Atmo%vz(iz))
                  lvel = vel.gt.TINYVEL

                end if

                ! For each direction
                do jdir=1,njdir

                  ! Skip static multiple directions
                  if (.not.lvel.and.jdir.gt.1) cycle

                  ! For each transition
                  do jtran=1,Atom(ia)%ntran

                    ! If absent, skip
                    if (Atom(ia)%fflag(jtran)%absent) cycle

                    ! For each transition (FS)
                    do fjtran=1,Atom(ia)%fst(jtran)%nt

                      ! Rolling index
                      ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                      ffktran = ffjtran + Atom(ia)%tfshift

                      ! Get index
                      indx = Red%idzao(ffktran,iz,jdir)

                      ! If was saving but cannot
                      if (Red%dzao(indx)%VRAM.and. &
                          checkram(jtran,iz,jdir).lt.0) then

                        ! Not storing now
                        Red%dzao(indx)%VRAM = .False.

                        ! Deallocate
                        RAM = RAM - 1d-6*sizeof(Red%dzao(indx)%p) + &
                              8d-6
                        deallocate(Red%dzao(indx)%p)

                      end if ! Was storing and not now

                    end do ! transition (FS)
                  end do ! transition
                end do ! height
              end do ! direction

            end if ! Slaves
          end if ! If saving in RAM
        end if ! MPI

        !
        ! Calculate multiplicative factor (instead of division factor)
        !

        ! Allocate and initialize out of bound counter
        allocate(outofbound(Atom(ia)%nftran))
        outofbound = 0

        ! If MPI, master does not need this
        if (MPID%mpi.and.pid.eq.0) then

          ! Run over indexes
          do indx=1,Red%ndzao

            ! Deallocate
            deallocate(Red%dzao(indx)%Norm)

          end do

          ! Free
          deallocate(checkram)

        ! Serial or slave
        else

          ! For each height
          do iz=Rz0,Rz1

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! Skip absent
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! For each transition (FS)
                do fjtran=1,Atom(ia)%fst(jtran)%nt

                  ! Rolling index
                  ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                  ffktran = ffjtran + Atom(ia)%tfshift

                  ! Get index
                  indx = Red%idzao(ffktran,iz,jdir)

                  ! Easier to write variable
                  d1 = Red%dzao(indx)%Norm(1)

                  ! If the norm is not zero
                  if (d1.gt.TINYN) then

                    ! Check close to 1
                    if (d1.lt.BADNORM.or.d1.gt.2d0-BADNORM) then

                      ! Add count of bad norm
                      outofbound(ffjtran) = outofbound(ffjtran) + 1

                      ! If can write, do it
                      if (obadnorm) &
                        call writebadbound(folder, &
                                           Atom(ia)%Element,iz,jdir, &
                                           .False.,jtran,fjtran,d1)
                    end if ! Bad norm
                  end if ! Norm not zero

                  ! If storing
                  if (Red%dzao(indx)%VRAM) then

                    ! If the norm is not zero
                    if (d1.gt.TINYN) then

                      ! Apply norm
                      Red%dzao(indx)%p = Red%dzao(indx)%p/d1

                      ! And deallocate the norm
                      deallocate(Red%dzao(indx)%Norm)

                    end if ! Norm not zero

                  ! Not storing
                  else

                    ! If the norm is not zero, store inverse
                    if (d1.gt.TINYN) Red%dzao(indx)%Norm = 1d0/d1

                  end if ! Storing

                end do ! FS transitions
              end do ! transitions
            end do ! Directions
          end do ! heights

        end if ! Master and MPI

        ! If MPI
        if (MPID%mpi) then

          ! Add together all bad normalization data

          ! Master
          if (pid.eq.0) then

            ! Master in place
            call MPI_REDUCE(MPI_IN_PLACE,outofbound, &
                            Atom(ia)%nftran,MPI_INTEGER, &
                            MPI_SUM,0,MPI_COMM_RT,ierr)
          ! Slave
          else

            ! Slave just send
            call MPI_REDUCE(outofbound,outofbound, &
                            Atom(ia)%nftran,MPI_INTEGER, &
                            MPI_SUM,0,MPI_COMM_RT,ierr)

          end if ! Master/slave
        end if ! MPI

        ! Master
        if (pid.eq.0) then

          ! Check bad limits
          if (maxval(outofbound).gt.0) then

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! For each FS transition
              do fjtran=1,Atom(ia)%fst(jtran)%nt

                ! Rolling index
                ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

                ! Check line is affected
                if (outofbound(ffjtran).gt.0) then

                  ! Write warning
                  write(umsg,'(A,i4,",",i4,4A,i4,A)') &
                    ' # Warning: transition ',jtran,fjtran,' in ', &
                    Atom(ia)%Element,' has bad normalization', &
                    ' for the chosen width in ', &
                    outofbound(ffjtran),' heights X directions.'
                  call verbose

                end if

              end do ! FS transition
            end do ! Transition

          end if ! Any transition had bad normalization
        end if ! Master

        ! Deallocate
        deallocate(outofbound)

      end do ! Atoms

      !
      ! LTE lines
      !

      ! LTE lines present
      if (allocated(LTElines)) then

        ! For each LTE line
        do ia=1,size(LTElines)

          ! Skip absent
          if (LTElines(ia)%absent) cycle

          !
          ! Allocate the norm array
          !

          ! For each height
          do iz=LTElines(ia)%Rz0,Rz1

            ! If dynamic
            if (dyn) then

              ! Check local velocity
              vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                         Atmo%vy(iz)*Atmo%vy(iz) + &
                         Atmo%vz(iz)*Atmo%vz(iz))
              lvel = vel.gt.TINYVEL

            end if

            ! For each direction
            do jdir=1,njdir

              ! Skip static multiple directions
              if (.not.lvel.and.jdir.gt.1) cycle

              ! Skip slave
              if (MPID%mpi.and.pid.eq.0) cycle

              ! Get index
              indx = Red%idzao(nxt+ia,iz,jdir)

              ! Allocate profile itself if storing and it is present
              if (lVIRAM) then

                ! Prediction
                ! Do not subtract because LTE lines never store
                ! the normalization constant
                d1 = 16d-6*dble(LTElines(ia)%if1 - &
                                LTElines(ia)%if0 + 1)

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored
                  Red%dzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%dzao(indx)%VRAM = .True.

                  ! Allocate
                  allocate(Red%dzao(indx)%p(LTElines(ia)%if0: &
                                            LTElines(ia)%if1))

                  ! Update RAM
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%dzao(indx)%VRAM = .False.

              end if ! Storing

            end do ! directions
          end do ! heights

          !
          ! SLAVE OR SINGLE PROCESSOR
          !
          if (.not.MPID%mpi.or.pid.gt.0) then

            !
            ! Calculate normalization
            !

            ! For each height
            do iz=LTElines(ia)%Rz0,Rz1

              ! Thermal part of the Doppler width
              DwT = sqrt(LTElines(ia)%cDopp*LTElines(ia)%cDopp* &
                         Atmo%T(iz) + &
                         Atmo%vmi(iz)*Atmo%vmi(iz))

              ! If dynamic
              if (dyn) then

                ! Check local velocity
                vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                           Atmo%vy(iz)*Atmo%vy(iz) + &
                           Atmo%vz(iz)*Atmo%vz(iz))
                lvel = vel.gt.TINYVEL

              end if

              ! For each direction
              do jdir=1,njdir

                ! Skip static multiple directions
                if (.not.lvel.and.jdir.gt.1) cycle

                ! Get index
                indx = Red%idzao(nxt+ia,iz,jdir)

                ! If not storing, why bother
                if (.not.Red%dzao(indx)%VRAM) cycle

                ! Dopple shift init
                vfac = 1d0

                ! If velocity
                if (lvel) then

                  ! If emergent
                  if (LOS) then

                    ! Trigonometry
                    ct = Geom%L_mu(ith)
                    st = sqrt(1d0 - ct*ct)
                    cc = cos(Geom%L_phi(iph))
                    sc = sin(Geom%L_phi(iph))

                  ! If quadrature
                  else

                    ! Recover the indexes
                    ith1 = Geom%ithv(jdir)
                    iph1 = Geom%iphv(jdir)

                    ! Trigonometry
                    ct = Geom%V_mu(ith1)
                    st = sqrt(1d0 - ct*ct)
                    cc = Geom%v_mux(iph1)
                    sc = Geom%v_muy(iph1)*sqrt(1d0 - cc*cc)

                  end if ! LOS or quadrature

                  ! Get Doppler shift
                  vfac = 1d0 - Atmo%vx(iz)*st*cc - &
                               Atmo%vy(iz)*st*sc - &
                               Atmo%vz(iz)*ct

                end if ! Velocity

                ! Output Doppler width
                Dw = LTElines(ia)%Dfreq*DwT
                iDw = 1d0/Dw

                ! Get indexes
                if0 = LTElines(ia)%if0
                if1 = LTElines(ia)%if1

                ! Common quantities
                Dfreqw = (LTElines(ia)%eu - LTElines(ia)%el)*iDw
                vfacw = vfac*iDw
                auf = LTElines(ia)%damp(iz)*iDw

                !
                ! Calculate profile
                !

                ! For each frequency
                do ifreq=if0,if1

                  ! Voigt
                  call voigtI(Dfreqw - Frec%omega(ifreq)*vfacw, &
                              auf,prof)

                  Red%dzao(indx)%p(ifreq) = prof

                end do ! frequencies
              end do ! output direction
            end do ! height

          end if ! MPI

        end do ! LTE lines

      end if ! LTE lines present

      ! And now update the memory in Normp
      call cram_red_norm(Red,num)
      VRAMc = VRAMc + num
      DRAMc = 0d0

      ! Control
      call control

      return

      end subroutine normalizeI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Normalize the absorption profiles for the first order
      !! profiles in the PRD emissivity for the intensity problem\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!         ofram(logical): If reached the RAM limit
      subroutine normalizeI_PRD(Atom,Atmo,Frec,Red,ofram)

      ! I/O

      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      logical, intent(out):: ofram

      ! Local

      logical:: LVRAM

      integer:: ia,jtran,fjtran,ffjtran,itermf,itermu,iJf,iJu
      integer:: iz,indx,ifreq,if0,if1,nf

      double precision:: RAM,prof,d1,W0,W1
      double precision:: el,eu,au,af,auf,atuf,Dfreq,DwT,Dw,iDw
      double precision, dimension(:), allocatable:: W0_0,W1_0


      ! Routine name
      urou = 'normalizeI_PRD'

      ! Initialize
      ofram = .False.

      ! Allocate space for data
      allocate(Red%pzao(Red%nzao))

      ! Initialize RAM numbers
      RAM = cram_add(1)

      ! For each atom
      do ia=1,nA

        !
        ! Share limit weights from master
        !

        ! Allocate extremal weights
        allocate(W0_0(Atom(ia)%ntran),W1_0(Atom(ia)%ntran))

        ! Master gets them from atom structure
        if (pid.eq.0) then
          W0_0 = Atom(ia)%W0
          W1_0 = Atom(ia)%W1
        end if

        ! Share with the rest
        call MPI_BCAST(W0_0,Atom(ia)%ntran,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_RT,ierr)
        call MPI_BCAST(W1_0,Atom(ia)%ntran,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_RT,ierr)

        !
        ! Allocate the norm array
        !

        ! For each height
        do iz=Rz0,Rz1_PRD

          ! Thermal part of the Doppler width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

          ! For each transition
          do jtran=1,Atom(ia)%ntran

            ! Skip if no PRD
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! Find the term indexes for this transition
            itermf = Atom(ia)%fst(jtran)%iterml
            itermu = Atom(ia)%fst(jtran)%itermu

            ! Get contributions to damping parameter
            au = Atom(ia)%damp(itermu,iz)
            af = Atom(ia)%damp(itermf,iz)
            auf = Atom(ia)%ldamp(jtran,iz)

            ! For each FS transition
            do fjtran=1,Atom(ia)%fst(jtran)%nt

              ! Reset limits
              W0 = W0_0(jtran)
              W1 = W1_0(jtran)

              ! Global index
              ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

              ! Index
              indx = Red%izao(ffjtran,ia,iz)

              ! Limits
              if0 = Red%zao(indx)%Igf0
              if1 = Red%zao(indx)%Igf1
              nf = if1 - if0 + 1

              ! Initialize norm
              allocate(Red%pzao(indx)%Norm(1))
              Red%pzao(indx)%Norm = 0d0

              ! Correct weights
              if (if0.ne.Atom(ia)%tif0(jtran)) &
                W0 = Frec%W_freq(if0)
              if (if1.ne.Atom(ia)%tif1(jtran)) &
                W1 = Frec%W_freq(if1)

              ! Allocate profile itself if storing and it is present
              if (VIRAM.and.nf.gt.0) then

                ! Predict
                ! Subtract already counter for norm
                d1 = 8d-6*dble(nf) - 8d-6

                ! If no more space
                if (floor(RAM+d1).gt.RLIM) then

                  ! No stored
                  Red%pzao(indx)%VRAM = .False.
                  ofram = .True.

                ! If there is space
                else

                  ! Storing
                  Red%pzao(indx)%VRAM = .True.

                  ! Allocate
                  allocate(Red%pzao(indx)%p(if0:if1))
                  RAM = RAM + d1

                end if ! Space to store

              ! Not storing Voigt
              else

                ! No stored
                Red%pzao(indx)%VRAM = .False.

              end if

              ! Find the level indexes for this transition
              iJf = Atom(ia)%fst(jtran)%ilevell(fjtran)
              iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)

              ! Get the FS energies
              eu = Atom(ia)%FSfreq(iJu,itermu)
              el = Atom(ia)%FSfreq(iJf,itermf)

              ! Output Doppler width
              Dw = (eu - el)*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)
              iDw = 1d0/Dw

              ! Common quantities
              Dfreq = eu - el
              d1 = 1d-5*iDw/sqrt(PI)
              atuf = (au+af+auf)*iDw

              ! If valid
              if (nf.gt.0) then

                ! Lower boundary
                call voigtI((Dfreq - Frec%omega(if0))*iDw,atuf,prof)
                Red%pzao(indx)%Norm = Red%pzao(indx)%Norm + &
                                      prof*W0*d1

                ! Store
                if (Red%pzao(indx)%VRAM) &
                  Red%pzao(indx)%p(if0) = prof

                ! For each frequency
                do ifreq=if0+1,if1-1

                  ! Add to the integral
                  call voigtI((Dfreq - Frec%omega(ifreq))*iDw, &
                              atuf,prof)
                  Red%pzao(indx)%Norm = Red%pzao(indx)%Norm + &
                                        prof*Frec%W_freq(ifreq)*d1

                  ! Store
                  if (Red%pzao(indx)%VRAM) &
                    Red%pzao(indx)%p(ifreq) = prof

                end do

                ! Upper boundary
                call voigtI((Dfreq - Frec%omega(if1))*iDw,atuf,prof)
                Red%pzao(indx)%Norm = Red%pzao(indx)%Norm + &
                                      prof*W1*d1

                ! Store
                if (Red%pzao(indx)%VRAM) &
                  Red%pzao(indx)%p(if1) = prof

              end if ! Valid frequencies

              ! If MPI
              if (nproc.gt.1) then

                ! Compute the norm
                call MPI_ALLREDUCE(MPI_IN_PLACE,Red%pzao(indx)%Norm, &
                                   1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                   MPI_COMM_RT,ierr)

                !
                ! Check everyone is saving
                !

                ! Valid range
                if (nf.gt.0) then

                  ! Actual variable
                  LVRAM = Red%pzao(indx)%VRAM

                ! No frequencies
                else

                  ! Send a true
                  LVRAM = .True.

                end if ! Valid range

                ! Check
                call MPI_ALLREDUCE(MPI_IN_PLACE,LVRAM, &
                                   1,MPI_LOGICAL,MPI_LAND, &
                                   MPI_COMM_RT,ierr)
                ! Valid range
                if (nf.gt.0) then

                  ! If was storing, but need to free
                  if (Red%pzao(indx)%VRAM.and..not.LVRAM) then

                    ! Free
                    ! Put back Norm memory
                    RAM = RAM - 1d-6*sizeof(Red%pzao(indx)%p) + 8d-6
                    deallocate(Red%pzao(indx)%p)

                  end if ! Was storing but cannot anymore
                end if ! Valid range
              end if ! MPI

              ! Larger than zero norm
              if (Red%pzao(indx)%Norm(1).gt.0d0) &
                Red%pzao(indx)%Norm(1) = 1d0/Red%pzao(indx)%Norm(1)

              ! If valid and storing
              if (nf.gt.0.and.Red%pzao(indx)%VRAM) then

                ! Larger than 0 norm
                if (Red%pzao(indx)%Norm(1).gt.0d0) &
                  Red%pzao(indx)%p = Red%pzao(indx)%p* &
                                     Red%pzao(indx)%Norm(1)
              end if

            end do ! For each FS transition
          end do ! For each transition
        end do ! For each height

        ! Free
        deallocate(W0_0,W1_0)

      end do ! Atoms

      ! Count RAM
      call cram_red_1stord(Red,RAM)
      ORAMc = RAM
      DRAM2c = 0d0

      ! Control
      call control

      return

      end subroutine normalizeI_PRD

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write warning about bad normalization of Voigt profiles in a
      !! file\n
      !!   folder(character(:)): Output folder path\n
      !!  element(character(:)): Name of the atomic element\n
      !!            iz(integer): Height index\n
      !!          jdir(integer): Direction index\n
      !!           pol(logical): If polarization profiles\n
      !!         field(logical): If there is magnetic field\m
      !!         jtran(integer): Transition index\n
      !!            jj(integer): Line component\n
      !!                         index within the term\n
      !!             d1(double): Normalization value
      subroutine writebadbound(folder,element,iz,jdir,pol,jtran,jj,d1)

      ! I/O
      character(len=2), intent(in):: element
      character(len=500), intent(in):: folder
      logical, intent(in):: pol
      integer, intent(in):: iz,jdir,jtran,jj
      double precision, intent(in):: d1

      ! Local

      logical:: exists

      character(LEN=5):: CPUC

      ! Get ID in character
      write(CPUC,'(I0.5)') gpid

      !
      ! Check if file exists
      inquire(file=trim(folder)//'/badnorm'//CPUC, exist=exists)

      ! If does not exist
      if(.not.exists)then

        ! Create new
        open(800,file=trim(folder)//'/badnorm'//CPUC)

      ! If it exists
      else

        ! Append to file
        open(800,file=trim(folder)//'/badnorm'//CPUC, &
             position='append')

      endif

      ! 1D case
      if (run_mode.eq.0) then

        ! If polarized case
        if (pol) then

          ! Write
          write(800,'("Atom",1x,A2,1x,'// &
                    '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                    '"Transition",1x,i4,1x,"Component",'// &
                    '1x,i3,1x,"Pol Norm",1x,f18.16)') &
              element,iz,jdir,jtran,jj,d1

        ! Intensity case
        else

          ! Write
          write(800,'("Atom",1x,A2,1x,'// &
                    '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                    '"Transition",1x,i4,1x,"Component",'// &
                    '1x,i3,1x,"Int Norm",1x,f18.16)') &
              element,iz,jdir,jtran,jj,d1

        end if ! Polarization or intensity

      ! Any other case
      else

        ! If polarized case
        if (pol) then

          ! Write
          write(800,'("Atom",1x,A2,1x,'// &
                      '"LOS",1x,"(",i4,",",i4,")",1x,'// &
                      '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                      '"Transition",1x,i4,1x,"Components",'// &
                      '1x,i3,1x,"Pol Norm",1x,f18.16)') &
              element,icoords(1:2),iz,jdir,jtran,jj,d1

        ! Intensity case
        else

          ! Write
          write(800,'("Atom",1x,A2,1x,'// &
                    '"LOS",1x,"(",i4,",",i4,")",1x,'// &
                    '"Height",1x,i3,1x,"Direction",1x,i3,1x,'// &
                    '"Transition",1x,i4,1x,"Component",'// &
                    '1x,i3,1x,"Int Norm",1x,f18.16)') &
              element,icoords(1:2),iz,jdir,jtran,jj,d1

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

      !> Dummy initialization of the normalization structure for the
      !! CLE synthesis mode\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!       Red(Red_class): Structure with redistribution input
      !!                       frequency data, redistribution function
      !!                       data, and profile or normalization data
      subroutine normalize_cle(Atom,Red)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Red_class), intent(inout):: Red

      ! Local

      integer:: ia,jtran,ncom


      ! Routine name
      urou = 'normalize_cle'

      ! Initialize largest number of components
      ncom = 1

      ! For each atom
      do ia=1,nA

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! Check larger
          if (Atom(ia)%trano(jtran)%ncomB.gt.ncom) &
            ncom = Atom(ia)%trano(jtran)%ncomB

        end do ! Transitions
      end do ! Atoms

      ! Only one point
      Red%ndzao = 1

      !
      ! Allocate norms
      allocate(Red%dzao(Red%ndzao))
      allocate(Red%dzao(1)%Norm(ncom))

      ! Dummy initialization
      Red%dzao(1)%VRAM = .False.
      Red%dzao(1)%Norm = 1d0
      VRAMc = 12d-6

      return

      end subroutine normalize_cle

!#####################################################################
!#####################################################################
!#####################################################################

      end module normalizer_mod
