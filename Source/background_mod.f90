      !> Background opacities
      module background_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/19/2017
!  Last version:
!     07/03/2023 V3.1.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     07/03/2023:    V3.1.0 - The calculation of the continuum
!                             opacity at the reference frequency has
!                             its own subroutine now (TdPA)
!
!     11/24/2022:    V3.0.3 - Added an exception to not remove the
!                             continuum data from the master if we
!                             are running CLE (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!
!     07/13/2022:    V3.0.1 - Changed the arguments in kurucz_bb call
!                             to accomodate other changes (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Now fudge factor and Kurucz lines
!                                are inputs to avoid reading the files
!                                repeatedly.
!                              o Consequently, the reading of the
!                                fudge and Kurucz files have been
!                                moved elsewhere.
!                              o The fudge is interpolated every
!                                time, though, as it is likely
!                                relatively lightweight.
!                              o Atmo%v has changed to Atmo%vx,%vy,
!                                and %vz.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!                           - Moved up the allocation of RT
!                             coefficients (TdPA)
!
!     09/24/2020:    V1.3.7 - When defining the vfac array for the
!                             bound-bound background contributions,
!                             the whole velocity vector was casted,
!                             instead of the portion in the relevant
!                             domain. Error found and solution
!                             proposed by David Afonso, IAC (TdPA)
!
!     09/11/2020:    V1.3.6 - Added memory counter for background
!                             quantities (TdPA)
!
!     05/11/2020:    V1.3.5 - Removed indexes in kurucz_bb call (TdPA)
!
!     03/05/2020:    V1.3.4 - kurucz_bb now asks for the total
!                             hydrogen density, also non-atomic (TdPA)
!                           - When the index where the hydrogen is
!                             stored is found in active or passive
!                             atoms, the complementary index is set
!                             to -1 to avoid undefined (TdPA)
!
!     12/17/2019:    V1.3.3 - Bugfix: The Thomson scattering for the
!                             reference frequency was using a
!                             quantity defined after its use (TdPA)
!
!     12/10/2019:    V1.3.2 - Always compute reference continuum
!                             opacity to transform to the other
!                             vertical axis (TdPA)
!
!     11/19/2019:    V1.3.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     09/13/2019:    V1.3.0 - When the height scale is in tau units,
!                             compute the absorptivity of the
!                             reference wavelength, only continuum,
!                             without b-b transitions (TdPA)
!
!     08/09/2019:    V1.2.2 - Option to skip bound-bound transitions
!                             without modifying the models (TdPA)
!
!     07/23/2019:    V1.2.1 - Moved Kurucz lines outside the
!                             background atom loop (TdPA)
!
!     04/08/2019:    V1.2.0 - Added possibility to add opacity from
!                             Kurucz line files (TdPA)
!
!     02/20/2019:    V1.1.0 - Changed verbosity (TdPA)
!                           - fudge read checks for success and unit
!                             is now 100 (TdPA)
!
!     11/28/2018:    V1.0.5 - Fix of doxygen comment (TdPA)
!
!     07/10/2018:    V1.0.4 - Need to pass Flgsg to back_bb (TdPA)
!
!     12/18/2017:    V1.0.3 - Bugfix: Two of the options to compute
!                             Doppler width for background atoms used
!                             the helium index (TdPA)
!
!     09/14/2017:    V1.0.2 - Added a path and ID to the file (TdPA)
!
!     06/12/2017:    V1.0.1 - Removed Cont%l (TdPA)
!
!     04/19/2017:    V1.0.0 - Started coding (TdPA)
!
!#####################################################################
!#####################################################################
!
!  Known bugs:
!
!    Not a bug itself, but never tested it with velocities
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  background:
!    Calculates the continuum opacity, scattering coefficient, and
!  emissivity
!
!  chi_freq:
!    Calculates the continuum opacity at the reference frequency,
!  neglecting b-b transitions
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
      use parameters_mod , ONLY : vacuum
      use planck_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes background quantities\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!        Input(Input_class): Structure with settings data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs
      subroutine background(Atom,Atomb,Mol,Atmo,fudge,kurucz,Input, &
                            omega,Cont,Geom,MPID,Flgsg)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atom_class), dimension(:), intent(in):: Atomb
      type(Mol_class), dimension(:), intent(in):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class):: Input
      type(kurucz_class):: kurucz
      type(Continuum_class):: Cont
      type(Geometry_class):: Geom
      type(MPI_class):: MPID
      type(Fctsg_class), intent(in):: Flgsg
      double precision, dimension(:), intent(in):: omega

      ! Local

      type(Continuum_class):: Cont_aux

      logical:: lH, lHe, lOH, lCH, lH2, nfline, fline, fkline

      integer:: if0,if1,ith,iph
      integer:: iaH,iabH,iaHe,iabHe,imOH,imCH, imH2
      integer:: ia,imol,ifreq,iz,idir,jdir
      integer:: ndir
      integer, dimension(:), allocatable:: ithV,iphV,ithLV,iphLV

      double precision:: freq, DwT
      double precision:: ct, st, cc, sc
      double precision, dimension(3):: fudge_fact
      double precision, dimension(:), allocatable:: eta, sig, eps, &
                                                    bplanck, vfac


      ! Routine name
      urou = 'background'

      ! The master does not care about the continuum
      ! unless is in CLE mode
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

      ! If there are no velocities, we have isotropy
      else

        ndir = 1
        vfac = 1d0

      end if

      ! Store the number of directions we may have to take into
      ! account
      Cont%ndir = ndir

      !
      ! Add Thomson scattering
      !
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
      ! The rest of contributions depends on frequencies
      !

      ! For each output frequency
      do ifreq=if0,if1

        ! Current frequency
        freq = omega(ifreq)

        ! Get the fudge factors for this frequency
        if (fudge%nfreq_f.gt.0) then
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,2), &
                      1d2/freq,fudge_fact(1))
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,3), &
                      1d2/freq,fudge_fact(2))
          call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,4), &
                      1d2/freq,fudge_fact(3))
        else
          fudge_fact = 1d0
        end if


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

        end if


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
            call rayleigh(freq,ifreq,Atomb(iabHe),DwT,1,nz,0,sig)
            do iz=1,nz
              Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
            end do

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
          call rayleigh(freq,ifreq,Atom(iaHe),DwT,1,nz,1,sig)
          do iz=1,nz
            Cont%c(ifreq,2,1,iz) = Cont%c(ifreq,2,1,iz) + sig(iz)
          end do

        end if


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

        end if


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

                ith = ithv(idir)

                if (ith.lt.0) then

                  ith = ithlv(idir)
                  iph = iphlv(idir)
                  ct = Geom%L_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = cos(Geom%L_phi(iph))
                  sc = sin(Geom%L_phi(iph))

                else

                  iph = iphv(idir)
                  ct = Geom%V_mu(ith)
                  st = sqrt(1d0 - ct*ct)
                  cc = Geom%v_mux(iph)
                  sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

                end if

                vfac(:) = 1d0 - atmo%vx*st*cc - &
                                atmo%vy*st*sc - &
                                atmo%vz*ct

              end if ! There are velocities

              ! If it is the hardwired HI model
              if (Atomb(ia)%cust) then

                call backH_bb(freq,Atomb(ia),Atmo%T, &
                              Atmo%vmi,DwT,vfac,1,nz, &
                              fline,eta,eps)

              ! It is a read model atom
              else

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

              call kurucz_bb(freq,Kurucz,Atmo,vfac, &
                             1,nz,fkline,eta,eps)

            else

              fkline = .False.

            end if

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

            end if

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

      ! Control
      call control

      end subroutine background

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes continuum opacity at a single frequency\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!            Mol(Mol_class): Structure with the molecule data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!        Input(Input_class): Structure with settings data\n
      !!              freq(dfloat): Frequency\n
      !!            chi(dfloat(:)): Background opacity\n
      !!              iz0(integer): First height index\n
      !!              iz1(integer): Last height index\n
      !!      skip_master(logical): If master can skip
      subroutine chi_freq(Atom,Atomb,Mol,Atmo,fudge, &
                          Input,freq,chi,iz0,iz1,skip_master)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atom_class), dimension(:), intent(in):: Atomb
      type(Mol_class), dimension(:), intent(in):: Mol
      type(Atmo_class), intent(inout):: Atmo
      type(fudge_class), intent(in):: fudge
      type(Input_class):: Input
      logical, intent(in):: skip_master
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(:), intent(out):: chi

      ! Local
      logical:: lH, lHe, lOH, lCH, lH2

      integer:: iaH,iabH,iaHe,iabHe,imOH,imCH, imH2
      integer:: ia,imol,ifreq,lnz

      double precision:: DwT
      double precision, dimension(3):: fudge_fact
      double precision, dimension(:), allocatable:: eta, sig, eps


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

      ! Get the fudge factors for this frequency
      if (fudge%nfreq_f.gt.0) then
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,2), &
                    1d2/freq,fudge_fact(1))
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,3), &
                    1d2/freq,fudge_fact(2))
        call linear(fudge%fudge_v(:,1),fudge%fudge_v(:,4), &
                    1d2/freq,fudge_fact(3))
      else
        fudge_fact = 1d0
      end if


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

      end if


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

      end if


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

      end if


      !
      ! Compute Background atoms bound-free contribution
      !
      do ia=1,nAb

        ! Skip the HI
        if (ia.eq.iabH) cycle

        call back_bf(freq,Atomb(ia),Atmo%T,iz0,iz1,eta,eps)
        chi(1:lnz) = chi(1:lnz) + eta(iz0:iz1)*fudge_fact(3)

      end do

      ! Add vaccum value
      chi(1:lnz) = chi(1:lnz) + vacuum

      ! Control
      call control

      return

      end subroutine chi_freq

!#####################################################################
!#####################################################################
!#####################################################################

      end module background_mod
