      !> Background opacities
      module background_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Revised headers (TdPA)
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
!     Manage the need for directional dependency in a more clever way
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  background
!    Computes background absortivity, scattering coefficient, and
!  emissivity
!
!  chi_freq
!    Computes background absortivity, scattering coefficient, and
!  emissivity at a single frequency neglecting b-b transitions
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use backgroundaux_mod
      use commons_mod
      use inter_mod
      use kurucz_mod
      use parameters_mod , ONLY : vacuum , TINYVEL
      use planck_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes background absortivity, scattering coefficient, and
      !! emissivity\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!   Atomb(Atom_class(:)): Structures with atomic data for
      !!                         background atoms\n
      !!      Mol(Mol_class(:)): Structures with molecular data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!     fudge(fudge_class): Structure with fudge data\n
      !!   kurucz(kurucz_class): Structure with Kurucz line data\n
      !!     Input(Input_class): Structure with configuration data\n
      !!       omega(double(:)): Frequency array\n
      !!  Cont(Continuum_class): Structure with background opacity
      !!                         data\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      subroutine background(Atom,Atomb,Mol,Atmo,fudge,kurucz,Input, &
                            omega,Cont,Geom,MPID,Flgsg)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atom_class), dimension(:), intent(in):: Atomb
      type(Mol_class), dimension(:), intent(in):: Mol
      type(Atmo_class), intent(in):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class), intent(in):: Input
      type(kurucz_class), intent(in):: kurucz
      type(Continuum_class), intent(inout):: Cont
      type(Geometry_class), intent(in):: Geom
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(in):: Flgsg
      double precision, dimension(:), intent(in):: omega

      ! Local

      type(Continuum_class):: Cont_aux

      logical:: lH,lHe,lOH,lCH,lH2,nfline,fline,fkline

      integer:: if0,if1,ith,iph
      integer:: iaH,iabH,iaHe,iabHe,imOH,imCH, imH2
      integer:: ia,imol,ifreq,iz,idir,jdir,ndir
      integer, dimension(:), allocatable:: ithV,iphV,ithLV,iphLV

      double precision:: freq,DwT,ct,st,cc,sc
      double precision, dimension(3):: fudge_fact
      double precision, dimension(:), allocatable:: eta,sig,eps, &
                                                    bplanck,vfac


      ! Routine name
      urou = 'background'

      ! The master does not care about the continuum
      ! unless it is in CLE mode
      if (MPID%mpi.and.pid.eq.0.and.run_mode.ne.2) then
        call control
        return
      end if

      ! Store the MPI limits into easy to handle variables
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)

      ! Allocate
      allocate(eta(nz),sig(nz),eps(nz),bplanck(nz),vfac(nz))

      ! Allocate the continuum structure, we assume first that
      ! we will not need the extra directions
      allocate(Cont%c(if0:if1,3,1,nz))

      ! Initialize
      Cont%c = 0d0

      ! Check if H is not being calculated (if it is not active)
      lH = .True.
      do ia=1,NA
        if (Atom(ia)%Element.eq.' H') then
          iaH = ia
          lH = .False.
          iabH = -1
          exit
        end if
      end do

      ! If HI not being calculated, find the index in the list of
      ! backgrounds, it should be the first one anyways
      if (lH) then
        do ia=1,NAb
          if (Atomb(ia)%element.eq.' H') then
            iabH = ia
            iaH = -1
            exit
          end if
        end do
      end if

      ! Check if He is being calculated
      lHe = .True.
      do ia=1,NA
        if (Atom(ia)%Element.eq.'HE') then
          iaHe = ia
          lHe = .False.
          exit
        end if
      end do

      ! If He not being calculated, find the index in the list of
      ! backgrounds
      if (lHe) then
        iabHe = -1
        do ia=1,NAb
          if (Atomb(ia)%element.eq.'HE') then
            iabHe = ia
            exit
          end if
        end do
      end if

      ! Check if there is OH
      lOH = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'OH') then
            lOH = .True.
            imOH = imol
            exit
          end if
        end do
      end if

      ! Check if there is CH
      lCH = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'CH') then
            lCH = .True.
            imCH = imol
            exit
          end if
        end do
      end if

      ! Check if there is H2
      lH2 = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'H2') then
            lH2 = .True.
            imH2 = imol
            exit
          end if
        end do
      end if

      !
      ! Allocate the continuum quantities
      !

      ! If there are velocities
      if (dyn) then

        ! We may need to specify the quantities for each direction
        ndir = Geom%nPh*Geom%nTh + Geom%nPhLOS*Geom%nThLOS

        !
        ! De-index the directions
        !
        allocate(ithv(ndir))
        allocate(iphv(ndir))
        allocate(ithlv(ndir))
        allocate(iphlv(ndir))
        idir = 0
        do ith=1,Geom%nTh
          do iph=1,Geom%nPh
            idir = idir + 1
            ithv(idir) = ith
            iphv(idir) = iph
            ithlv(idir) = -1
            iphlv(idir) = -1
          end do
        end do
        do ith=1,Geom%nThLOS
          do iph=1,Geom%nPhLOS
            idir = idir + 1
            ithv(idir) = -1
            iphv(idir) = -1
            ithlv(idir) = ith
            iphlv(idir) = iph
          end do
        end do

      ! If there are no velocities
      else

        ! Isotropy
        ndir = 1
        vfac = 1d0

      end if

      ! Store the number of directions we may have to take into
      ! account
      Cont%ndir = ndir

      !
      ! Add Thomson scattering
      !

      ! Get scattering coefficient
      call Thomson(Atmo%ne,1,nz,sig)
      do iz=1,nz
        Cont%c(:,2,1,iz) = sig(iz)
      end do


      !
      ! Initialize nfline, flag to indicates that we haven't find any
      ! spectral line in the background
      !
      nfline = .True.


      !
      ! The rest of contributions depend on frequencies
      !

      ! For each output frequency
      do ifreq=if0,if1

        ! Current frequency
        freq = omega(ifreq)

        ! If there are fudge factors
        if (fudge%nfreq_f.gt.0) then

          ! Interpolate them linearly
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,2), &
                      1d2/freq,fudge_fact(1))
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,3), &
                      1d2/freq,fudge_fact(2))
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,4), &
                      1d2/freq,fudge_fact(3))

        ! No fudge factors
        else

          ! Set factors to 1
          fudge_fact = 1d0

        end if ! There are fudge factors


        ! Calculate the Planck function for each height at this
        ! frequency
        do iz=1,nz
          bplanck(iz) = planck(freq,Atmo%T(iz))
        end do


        !
        ! Compute H- b-f contribution
        !
        call Hminus_bf(freq,Atmo%nHm,Atmo%T,1,nz,eta,eps)
        do iz=1,nz
          Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
          Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + eps(iz)
        end do


        !
        ! Compute H- f-f contritubion
        !
        call Hminus_ff(freq,Atmo%nh(1:nz,1),Atmo%ne, &
                       Atmo%T,1,nz,eta)
        do iz=1,nz
          Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
          Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                 eta(iz)*bplanck(iz)
        end do


        !
        ! Apply H- opacity fudge
        !
        Cont%c(ifreq,1,1,:) = Cont%c(ifreq,1,1,:)*fudge_fact(1)
        Cont%c(ifreq,3,1,:) = Cont%c(ifreq,3,1,:)*fudge_fact(1)


        !
        ! Compute OH b-f contribution if it is present
        !
        if (lOH) then
          call OH_bf(freq,Mol(imOH)%n,Atmo%T,1,nz,eta,eps)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do
        end if


        !
        ! Compute CH b-f contribution if it is present
        !
        if (lCH) then
          call CH_bf(freq,Mol(imCH)%n,Atmo%T,1,nz,eta,eps)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do
        end if


        !
        ! Contributions of the HI atom
        !

        ! If H is passive
        if (lH) then

          !
          ! Compute HI b-f contribution
          !
          call HI_bf(freq,Atomb(iabH),Atmo%T,1,nz,eta,eps)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + eps(iz)
          end do


          !
          ! Compute HI f-f contribution
          !
          call HI_ff(freq,Atomb(iabH),Atmo%T,Atmo%ne,1,nz,eta)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do


          !
          ! Compute H + p+ f-f contribution
          !
          call HHp_ff(freq,Atomb(iabH),Atmo%T,1,nz,eta)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do


          !
          ! Compute HI Rayleigh contribution
          !

          ! Compute Doppler width
          if(Input%dws.eq.'MAX')then
            DwT = Atomb(iabH)%cDopp*sqrt(maxval(Atmo%T))
          else if(Input%dws.eq.'MIN')then
            DwT = Atomb(iabH)%cDopp*sqrt(minval(Atmo%T))
          else if(Input%dws.eq.'NUM')then
            DwT = Input%dw*1d-9/c
          end if

          ! Compute proper scattering coefficient
          call rayleigh(freq,ifreq,Atomb(iabH),DwT,1,nz,0,sig)
          do iz=1,nz
            Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
          end do


        ! If H is active
        else


          !
          ! Compute HI f-f contribution
          !
          call HI_ff(freq,Atom(iaH),Atmo%T,Atmo%ne,1,nz,eta)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do


          !
          ! Compute H + p+ f-f contribution
          !
          call HHp_ff(freq,Atom(iaH),Atmo%T,1,nz,eta)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do


          !
          ! Compute HI Rayleigh contribution
          !

          ! Compute Doppler width
          if(Input%dws.eq.'MAX')then
            DwT = Atom(iabH)%cDopp*sqrt(maxval(Atmo%T))
          else if(Input%dws.eq.'MIN')then
            DwT = Atom(iabH)%cDopp*sqrt(minval(Atmo%T))
          else if(Input%dws.eq.'NUM')then
            DwT = Input%dw*1d-9/c
          end if

          ! Compute proper scattering coefficient
          call rayleigh(freq,ifreq,Atom(iaH),DwT,1,nz,1,sig)
          do iz=1,nz
            Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
          end do

        end if ! H active or not


        !
        ! Compute He I Rayleigh contribution
        !

        ! If He is not active
        if (lHe) then

          ! There is actual helium
          if (iabHe.gt.0) then

            ! Compute Doppler width
            if(Input%dws.eq.'MAX')then
              DwT = Atomb(iabHe)%cDopp*sqrt(maxval(Atmo%T))
            else if(Input%dws.eq.'MIN')then
              DwT = Atomb(iabHe)%cDopp*sqrt(minval(Atmo%T))
            else if(Input%dws.eq.'NUM')then
              DwT = Input%dw*1d-9/c
            end if

            ! Compute proper scattering coefficient
            call rayleigh(freq,ifreq,Atomb(iabHe),DwT,1,nz,0,sig)
            do iz=1,nz
              Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
            end do

          end if ! There is actual helium

        ! If He is active
        else

          ! Compute Doppler width
          if(Input%dws.eq.'MAX')then
            DwT = Atom(iabHe)%cDopp*sqrt(maxval(Atmo%T))
          else if(Input%dws.eq.'MIN')then
            DwT = Atom(iabHe)%cDopp*sqrt(minval(Atmo%T))
          else if(Input%dws.eq.'NUM')then
            DwT = Input%dw*1d-9/c
          end if

          ! Compute proper scattering coefficient
          call rayleigh(freq,ifreq,Atom(iaHe),DwT,1,nz,1,sig)
          do iz=1,nz
            Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
          end do

        end if ! Active He


        !
        ! Contributions of the H2 molecule, if present
        !
        if (lH2) then


          !
          ! Compute H2 f-f contribution
          !
          call H2m_ff(freq,Mol(imH2)%n,Atmo%T,Atmo%ne,1,nz,eta)
          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + eta(iz)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eta(iz)*bplanck(iz)
          end do


          !
          ! Compute H2 Rayleigh contribution
          !
          call rayleigh_H2(freq,Mol(imH2)%n,1,nz,sig)
          do iz=1,nz
            Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
          end do

        end if ! H2 is present


        !
        ! Compute Background atoms bound-free contribution
        !
        do ia=1,nAb

          ! Skip the HI
          if (ia.eq.iabH) cycle

          call back_bf(freq,Atomb(ia),Atmo%T,1,nz,eta,eps)

          do iz=1,nz
            Cont%c(ifreq,1,1,iz) = Cont%c(ifreq,1,1,iz) + &
                                   eta(iz)*fudge_fact(3)
            Cont%c(ifreq,3,1,iz) = Cont%c(ifreq,3,1,iz) + &
                                   eps(iz)*fudge_fact(3)
          end do
        end do


        !
        ! Apply the scattering fudge and add the scattering
        ! coefficient to the total opacity
        !
        Cont%c(ifreq,2,1,:) = Cont%c(ifreq,2,1,:)*fudge_fact(2)
        Cont%c(ifreq,1,1,:) = Cont%c(ifreq,1,1,:) + &
                              Cont%c(ifreq,2,1,:)


        !
        ! Contribution of spectral lines from background atoms
        !

        ! If including bound-bound transitions
        if (Input%addbb) then

          ! For each direction in the problem
          do idir=1,ndir

            ! For each background atom
            do ia=1,nAb

              ! Get the Doppler width
              if(Input%dws.eq.'MAX')then
                DwT = Atom(ia)%cDopp*sqrt(maxval(Atmo%T))
              else if(Input%dws.eq.'MIN')then
                DwT = Atom(ia)%cDopp*sqrt(minval(Atmo%T))
              else if(Input%dws.eq.'NUM')then
                DwT = Input%dw*1d-9/c
              end if

              ! If there are velocities, compute the frequency Doppler
              ! factor
              if (dyn) then

                ! Get polar index
                ith = ithv(idir)

                ! If LOS index
                if (ith.lt.0) then

                  ! Get real indexes and trigonometry
                  ith = ithlv(idir)
                  iph = iphlv(idir)
                  ct = Geom%L_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph))
                  sc = sin(Geom%L_phi(iph))

                ! If quadrature index
                else

                  ! Complete the indexes and get trigonometry
                  iph = iphv(idir)
                  ct = Geom%V_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph)
                  sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

                end if ! LOS/quadrature index

                ! Check velocity amplitude
                vfac = sqrt(Atmo%vx*Atmo%vx + Atmo%vy*Atmo%vy + &
                            Atmo%vz*Atmo%vz)

                ! For every height
                do iz=1,nz

                  ! If velocity
                  if (vfac(iz).gt.TINYVEL) then

                    vfac(iz) = 1d0 - atmo%vx(iz)*st*cc - &
                                     atmo%vy(iz)*st*sc - &
                                     atmo%vz(iz)*ct

                  ! No velocity
                  else

                    vfac(iz) = 1d0

                  end if ! Velocity at this height

                end do ! Heights

              end if ! There are velocities

              ! If it is the hard.coded HI model
              if (Atomb(ia)%cust) then

                ! Get bound-bound contribution
                call backH_bb(freq,Atomb(ia),Atmo%T, &
                              Atmo%vmi,DwT,vfac,1,nz, &
                              fline,eta,eps)

              ! It is a read model atom
              else

                ! Get bound-bound contribution
                call back_bb(freq,Atomb(ia),Flgsg,Atmo%T, &
                             Atmo%vmi,DwT,vfac,1,nz, &
                             fline,eta,eps)

              end if

              ! If we found a line for the first time and there are
              ! velocities
              if (fline.and.nfline.and.dyn) then

                ! Allocate an auxiliar and store the known data there
                allocate(Cont_aux%c(if0:if1,3,1,nz))
                Cont_aux%c = Cont%c

                ! Reallocate the true continuum structure
                deallocate(Cont%c)
                allocate(Cont%c(if0:if1,3,ndir,nz))

                ! Recover the data we already had
                do jdir=1,ndir
                  Cont%c(if0:if1,:,jdir,:) = Cont_aux%c(if0:if1,:,1,:)
                end do

                ! Deallocate the auxiliar
                deallocate(Cont_aux%c)

                ! Flag that we already found a line
                nfline = .False.

              end if

              ! If there was a line in the background atom, add its
              ! contributions
              if (fline) then
                do iz=1,nz
                  Cont%c(ifreq,1,idir,iz) = eta(iz) + &
                                            Cont%c(ifreq,1,idir,iz)
                  Cont%c(ifreq,3,idir,iz) = eps(iz) + &
                                            Cont%c(ifreq,3,idir,iz)
                end do
              end if

            end do ! Background atoms

            ! If we read Kurucz lines
            if (Kurucz%ntran.ge.1) then

              ! Get bound-bound contributions
              call kurucz_bb(freq,Kurucz,Atmo,vfac, &
                             1,nz,fkline,eta,eps)

            ! No Kurucz lines
            else

              ! Flag as not lines found
              fkline = .False.

            end if ! Kurucz data read

            ! If we found a line for the first time and there are
            ! velocities
            if (fkline.and.nfline.and.dyn) then

              ! Allocate an auxiliar and store the known data there
              allocate(Cont_aux%c(if0:if1,3,1,nz))
              Cont_aux%c = Cont%c

              ! Reallocate the true continuum structure
              deallocate(Cont%c)
              allocate(Cont%c(if0:if1,3,ndir,nz))

              ! Recover the data we already had
              do jdir=1,ndir
                Cont%c(if0:if1,:,jdir,:) = Cont_aux%c(if0:if1,:,1,:)
              end do

              ! Deallocate the auxiliar
              deallocate(Cont_aux%c)

              ! Flag that we already found a line
              nfline = .False.

            end if ! First time we find a line

            ! If there was a line in the kurucz list, add its
            ! contributions
            if (fkline) then
              do iz=1,nz
                Cont%c(ifreq,1,idir,iz) = eta(iz) + &
                                          Cont%c(ifreq,1,idir,iz)
                Cont%c(ifreq,3,idir,iz) = eps(iz) + &
                                          Cont%c(ifreq,3,idir,iz)
              end do
            end if

            ! If there are no lines, but there are dynamics, it is
            ! still isotropic
            if (nfline.and.dyn) then
              Cont%ndir = 1
              exit
            end if

          end do ! Directions

        ! Not including bound-bound transitions
        else

          ! No b-b lines
          Cont%ndir = 1

        end if ! Including b-b transitions

      end do ! Frequencies

      ! Count memory in background quantities
      BRAMc = 1d-6*sizeof(Cont%c)

      ! Free
      deallocate(eta,sig,eps,bplanck,vfac)
      if (dyn) deallocate(ithv,iphv,ithlv,iphlv)

      ! Control
      call control

      end subroutine background

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes background absortivity, scattering coefficient, and
      !! emissivity at a single frequency neglecting b-b transitions\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms\n
      !!     Mol(Mol_class(:)): Structures with molecular data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!    fudge(fudge_class): Structure with fudge data\n
      !!    Input(Input_class): Structure with configuration data\n
      !!          freq(double): Frequency\n
      !!        chi(double(:)): Background opacity\n
      !!          iz0(integer): First height index to consider\n
      !!          iz1(integer): Last height index to consider\n
      !!  skip_master(logical): If master can skip this
      subroutine chi_freq(Atom,Atomb,Mol,Atmo,fudge, &
                          Input,freq,chi,iz0,iz1,skip_master)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atom_class), dimension(:), intent(in):: Atomb
      type(Mol_class), dimension(:), intent(in):: Mol
      type(Atmo_class), intent(in):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class), intent(in):: Input
      logical, intent(in):: skip_master
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(:), intent(out):: chi

      ! Local

      logical:: lH,lHe,lOH,lCH,lH2

      integer:: iaH,iabH,iaHe,iabHe,imOH,imCH,imH2
      integer:: ia,imol,ifreq,lnz

      double precision:: DwT
      double precision, dimension(3):: fudge_fact
      double precision, dimension(:), allocatable:: eta,sig,eps


      ! The master does not care about the continuum
      ! unless is in CLE mode
      if (skip_master.and.pid.eq.0.and.run_mode.ne.2) then
        call control
        return
      end if

      ! Local size of chi
      lnz = iz1 - iz0 + 1

      ! Allocate eta, sig, eps
      allocate(eta(nz),sig(nz),eps(nz))

      ! Check if H is not being calculated (if it is not active)
      lH = .True.
      do ia=1,NA
        if (Atom(ia)%Element.eq.' H') then
          iaH = ia
          lH = .False.
          iabH = -1
          exit
        end if
      end do

      ! If HI not being calculated, find the index in the list of
      ! backgrounds, it should be the first one anyways
      if (lH) then
        do ia=1,NAb
          if (Atomb(ia)%element.eq.' H') then
            iabH = ia
            iaH = -1
            exit
          end if
        end do
      end if

      ! Check if He is being calculated
      lHe = .True.
      do ia=1,NA
        if (Atom(ia)%Element.eq.'HE') then
          iaHe = ia
          lHe = .False.
          exit
        end if
      end do

      ! If He not being calculated, find the index in the list of
      ! backgrounds
      if (lHe) then
        iabHe = -1
        do ia=1,NAb
          if (Atomb(ia)%element.eq.'HE') then
            iabHe = ia
            exit
          end if
        end do
      end if

      ! Check if there is OH
      lOH = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'OH') then
            lOH = .True.
            imOH = imol
            exit
          end if
        end do
      end if

      ! Check if there is CH
      lCH = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'CH') then
            lCH = .True.
            imCH = imol
            exit
          end if
        end do
      end if

      ! Check if there is H2
      lH2 = .False.
      if (nM.gt.0) then
        do imol=1,nM
          if (Mol(imol)%Molecule.eq.'H2') then
            lH2 = .True.
            imH2 = imol
            exit
          end if
        end do
      end if

      ! If there are fudge factors
      if (fudge%nfreq_f.gt.0) then

        ! Interpolate them linearly
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,2), &
                    1d2/freq,fudge_fact(1))
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,3), &
                    1d2/freq,fudge_fact(2))
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,4), &
                    1d2/freq,fudge_fact(3))

      ! No fudge factors
      else

        ! Set factors to 1
        fudge_fact = 1d0

      end if ! There are fudge factors


      ! Add Thomson scattering
      call Thomson(Atmo%ne,iz0,iz1,sig)
      chi(1:lnz) = sig(iz0:iz1)*fudge_fact(2)

      !
      ! Compute H- b-f contribution
      !
      call Hminus_bf(freq,Atmo%nHm,Atmo%T,iz0,iz1,eta,eps)
      chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

      !
      ! Compute H- f-f contritubion
      !
      call Hminus_ff(freq,Atmo%nh(1:nz,1),Atmo%ne,Atmo%T,iz0,iz1,eta)
      chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

      !
      ! Apply H- opacity fudge
      !
      chi(1:lnz) = chi(1:lnz)*fudge_fact(1)

      !
      ! Compute OH b-f contribution if it is present
      !
      if (lOH) then
        call OH_bf(freq,Mol(imOH)%n,Atmo%T,iz0,iz1,eta,eps)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)
      end if

      !
      ! Compute CH b-f contribution if it is present
      !
      if (lCH) then
        call CH_bf(freq,Mol(imCH)%n,Atmo%T,iz0,iz1,eta,eps)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)
      end if

      !
      ! Contributions of the HI atom
      !

      ! If H is passive
      if (lH) then

        !
        ! Compute HI b-f contribution
        !
        call HI_bf(freq,Atomb(iabH),Atmo%T,iz0,iz1,eta,eps)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

        !
        ! Compute HI f-f contribution
        !
        call HI_ff(freq,Atomb(iabH),Atmo%T,Atmo%ne,iz0,iz1,eta)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

        !
        ! Compute H + p+ f-f contribution
        !
        call HHp_ff(freq,Atomb(iabH),Atmo%T,iz0,iz1,eta)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

        !
        ! Compute HI Rayleigh contribution
        !

        ! Compute Doppler width
        if(Input%dws.eq.'MAX')then
          DwT = Atomb(iabH)%cDopp*sqrt(maxval(Atmo%T))
        else if(Input%dws.eq.'MIN')then
          DwT = Atomb(iabH)%cDopp*sqrt(minval(Atmo%T))
        else if(Input%dws.eq.'NUM')then
          DwT = Input%dw*1d-9/c
        end if

        ! Compute proper scattering coefficient
        call rayleigh(freq,ifreq,Atomb(iabH),DwT,iz0,iz1,0,sig)
        chi(1:lnz) = chi(1:lnz) + sig(iz0:iz1)*fudge_fact(2)

      ! If H is active
      else

        !
        ! Compute HI f-f contribution
        !
        call HI_ff(freq,Atom(iaH),Atmo%T,Atmo%ne,iz0,iz1,eta)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

        !
        ! Compute H + p+ f-f contribution
        !
        call HHp_ff(freq,Atom(iaH),Atmo%T,iz0,iz1,eta)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)

        !
        ! Compute HI Rayleigh contribution
        !

        ! Compute Doppler width
        if(Input%dws.eq.'MAX')then
          DwT = Atom(iabH)%cDopp*sqrt(maxval(Atmo%T))
        else if(Input%dws.eq.'MIN')then
          DwT = Atom(iabH)%cDopp*sqrt(minval(Atmo%T))
        else if(Input%dws.eq.'NUM')then
          DwT = Input%dw*1d-9/c
        end if

        ! Compute proper scattering coefficient
        call rayleigh(freq,ifreq,Atom(iaH),DwT,iz0,iz1,1,sig)
        chi(1:lnz) = chi(1:lnz) + sig(iz0:iz1)*fudge_fact(2)

      end if ! H active or not


      !
      ! Compute HeI Rayleigh contribution
      !

      ! If He is not active
      if (lHe) then

        ! There is actual helium
        if (iabHe.gt.0) then

          ! Compute Doppler width
          if(Input%dws.eq.'MAX')then
            DwT = Atomb(iabHe)%cDopp*sqrt(maxval(Atmo%T))
          else if(Input%dws.eq.'MIN')then
            DwT = Atomb(iabHe)%cDopp*sqrt(minval(Atmo%T))
          else if(Input%dws.eq.'NUM')then
            DwT = Input%dw*1d-9/c
          end if

          ! Compute proper scattering coefficient
          call rayleigh(freq,ifreq,Atomb(iabHe),DwT,iz0,iz1,0,sig)
          chi(1:lnz) = chi(1:lnz) + sig(iz0:iz1)*fudge_fact(2)

        end if

      ! If He is active
      else

        ! Compute Doppler width
        if(Input%dws.eq.'MAX')then
          DwT = Atom(iabHe)%cDopp*sqrt(maxval(Atmo%T))
        else if(Input%dws.eq.'MIN')then
          DwT = Atom(iabHe)%cDopp*sqrt(minval(Atmo%T))
        else if(Input%dws.eq.'NUM')then
          DwT = Input%dw*1d-9/c
        end if

        ! Compute proper scattering coefficient
        call rayleigh(freq,ifreq,Atom(iaHe),DwT,iz0,iz1,1,sig)
        chi(1:lnz) = chi(1:lnz) + sig(iz0:iz1)*fudge_fact(2)

      end if ! He active or not


      !
      ! Contributions of the H2 molecule, if present
      !
      if (lH2) then


        !
        ! Compute H2 f-f contribution
        !
        call H2m_ff(freq,Mol(imH2)%n,Atmo%T,Atmo%ne,iz0,iz1,eta)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)


        !
        ! Compute H2 Rayleigh contribution
        !
        call rayleigh_H2(freq,Mol(imH2)%n,iz0,iz1,sig)
        chi(1:lnz) = chi(1:lnz) + sig(iz0:iz1)*fudge_fact(2)

      end if ! There is H2


      !
      ! Compute Background atoms bound-free contribution
      !
      do ia=1,nAb

        ! Skip the HI
        if (ia.eq.iabH) cycle

        call back_bf(freq,Atomb(ia),Atmo%T,iz0,iz1,eta,eps)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)*fudge_fact(3)

      end do ! Background atoms

      ! Add vaccum value
      chi(1:lnz) = chi(1:lnz) + vacuum

      ! Free
      deallocate(eta,sig,eps)

      ! Control
      call control

      return

      end subroutine chi_freq

!#####################################################################
!#####################################################################
!#####################################################################

      end module background_mod
