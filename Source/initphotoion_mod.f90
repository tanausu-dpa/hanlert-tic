      !> Initialization of photoionization quantities
      module initphotoion_mod
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
!     10/16/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/16/2023:    V3.0.5 - The decision to free the photoionization
!                             cross-section comes from the argument to
!                             setphotoTEI, and not the mode (TdPA)
!
!     11/24/2022:    V3.0.4 - Added a condition for which the
!                             photoionization cross-section must not
!                             be erased if doing CLE (TdPA)
!
!     10/25/2022:    V3.0.3 - Implemented the limitation of the
!                             height axis (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!                           - Removed MPI%ierr variable (TdPA)
!
!     07/08/2022:    V3.0.1 - Changed a control call into a gcontrol
!                             call (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o setphoto has been split into
!                                setphoto, which computes the cross-
!                                section in the problem grid, and
!                                setphotoTEI, which computes the
!                                thermal ionization part that goes
!                                into the SEE.
!                              o setphoto does not require wfreq,
!                                T, or ne variables anymore.
!                              o The calculation of TEI needs MPI
!                                now because is called after splitting
!                                the frequency axes.
!                              o Frec%exu is nullified when not needed
!                                to avoid trying to deallocate it
!                                when freeing memory.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     09/11/2020:    V1.2.6 - The total RAM counter is not changed in
!                             this module anymore (TdPA)
!
!     07/31/2020:    V1.2.5 - Allocate a dummy Frec%exu(1,1) if not
!                             properly storing in ramphot (TdPA)
!
!     11/19/2019:    V1.2.4 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - It initializes RAM again (TdPA)
!
!     11/12/2019:    V1.2.3 - Moved elsewhere the initialization of
!                             MPI%RAM (TdPA)
!
!     08/08/2019:    V1.2.2 - Cross sections are interpolated in
!                             wavelength and not in frequencies (TdPA)
!
!     03/12/2019:    V1.2.1 - No need to store iexu (TdPA)
!                           - Added information to spline
!                             interpolation message (TdPA)
!                           - Master stores exponentials (TdPA)
!
!     02/20/2019:    V1.2.0 - New verbosity (TdPA)
!                           - Now uses diexp function (TdPA)
!
!     11/28/2018:    V1.1.3 - Removed a non-used parameter (TdPA)
!
!     09/05/2018:    V1.1.2 - Forgot to allocate the new variable
!                             iexu in ramphot (TdPA)
!
!     09/04/2018:    V1.1.1 - Now the inverse of the exponential is
!                             also stored in ramphoto (TdPA)
!
!     08/06/2018:    V1.1.0 - Added ramphoto that stores in RAM some
!                             quantities that are slow to compute in
!                             epsphoto and epsIphoto (TdPA)
!
!     09/08/2017:    V1.0.1 - If splines in gives a negative cross
!                             section, use linear instead (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
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
!  setphoto:
!    Interpolates cross section to the frequency axis
!
!  setphotoTEI:
!    Define the thermal contribution that goes into the SEE
!
!  ramphoto:
!    Allocates and pre-compute frequency quantities used in photoeps
!    and photoepsI
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use backgroundaux_mod
      use commons_mod
      use inter_mod
      use math_mod
      use parameters_mod , only :  c , ryd , c2,  cSaha, fktoJ, &
                                   kb , pi, bigexp , vbigexp , &
                                   vbigexpv
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes some photoionization related
      !! quantities independent of the model atmosphere.\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!     freq(dfloat(:)): Frequency array\n
      !!     MPID(MPI_class): Structure with MPI data
      subroutine setphoto(Atom,freq,MPID)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(MPI_class), intent(in):: MPID
      double precision, dimension(:), intent(in):: freq

      ! Local

      logical:: llinear

      integer:: nfr,ifreq,iphot,if0,if1,iterm,iterm1
      integer:: ilevel,ilevel1,iJ,iJ1

      double precision:: Z,neff,gbf,gbfe,cfreq
      double precision, dimension(:), allocatable:: sb,sc,sd,sx,sy


      ! For every pair of levels
      do ilevel=1,Atom%nlevel-1
        do ilevel1=ilevel+1,Atom%nlevel

          ! Check the indexing of photoionizations
          iphot = Atom%iphot(ilevel,ilevel1)

          ! If there is not, cycle to another pair
          if (iphot.lt.1) cycle

          ! Initialize logical
          llinear = .False.

          ! Get the transition limits into shorter variables for
          ! convenience
          if0 = Atom%phot(iphot)%if0
          if1 = Atom%phot(iphot)%if1

          ! Adjust integration weights in boundaries
          Atom%phot(iphot)%W0 = .5d5*(freq(if0+1) - freq(if0))
          Atom%phot(iphot)%W1 = .5d5*(freq(if1) - freq(if1-1))
          if (pid.eq.0.and.MPID%mpi) then
            Atom%phot(iphot)%MW0 = Atom%phot(iphot)%W0
            Atom%phot(iphot)%MW1 = Atom%phot(iphot)%W1
          end if

          ! Allocations
          ! Cross section
          allocate(Atom%phot(iphot)%alpha(if0:if1))

          ! If the input mode is explicit
          if (Atom%phot(iphot)%mode.eq.0) then

            ! Number of frequencies in the input atomic model
            nfr = Atom%phot(iphot)%nfreq

            ! Check that we have enough space for the spline
            ! coefficients
            if (allocated(sx)) then
              if (size(sx).lt.nfr) then
                deallocate(sx)
                deallocate(sy)
                deallocate(sb)
                deallocate(sc)
                deallocate(sd)
                allocate(sx(nfr))
                allocate(sy(nfr))
                allocate(sb(nfr))
                allocate(sc(nfr))
                allocate(sd(nfr))
              end if
            else
              allocate(sx(nfr))
              allocate(sy(nfr))
              allocate(sb(nfr))
              allocate(sc(nfr))
              allocate(sd(nfr))
            end if

            sx(nfr:1:-1) = 1d2/Atom%phot(iphot)%infreq
            sy(nfr:1:-1) = Atom%phot(iphot)%inalpha

            ! Get the Spline coefficients
            call spline(sx(1:nfr),sy(1:nfr), &
                        sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)

            ! Calculate the cross section for the relevant frequencies
            ! in cm^2
            do ifreq=if0,if1

              Atom%phot(iphot)%alpha(ifreq) = &
                ispline(1d2/freq(ifreq),sx(1:nfr),sy(1:nfr), &
                        sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)*1d4

              ! If something negative, do linear instead
              if (Atom%phot(iphot)%alpha(ifreq).lt.0d0) then
                if (pid.eq.0) then
                  write(umsg,'(A,1x,i4,A,1x,i4,1x,A,A,A)') &
                             ' # Spline interpolation of '// &
                             'ionization cross section of '// &
                             'b-f transition',ilevel,' -->', &
                             ilevel1,' of ', &
                             Atom%Element,' atom gave '// &
                             'negative value, doing linear.'
                  call verbose
                end if
                llinear = .True.
                exit
              end if

            end do

            if (llinear) then

              ! Calculate the cross section for the relevant
              ! frequencies in cm^2
              do ifreq=if0,if1

                call linear(sx(1:nfr),sy(1:nfr), &
                            1d2/freq(ifreq), &
                            Atom%phot(iphot)%alpha(ifreq))

                Atom%phot(iphot)%alpha(ifreq) = &
                                     Atom%phot(iphot)%alpha(ifreq)*1d4

              end do

            end if

          ! If the input mode is Hydrogenic
          else

            ! Find the term and J index of the levels involved and
            ! determine that charge and effective principal quantum
            ! number
            iterm1 = Atom%term(ilevel1)
            iJ1 = Atom%sublevel(ilevel1)
            iterm = Atom%term(ilevel)
            iJ = Atom%sublevel(ilevel)
            Z = dble(Atom%stage(iterm1) - 1)
            neff = Z*sqrt(ryd/ &
                   (Atom%FSfreq(iJ1,iterm1) - Atom%FSfreq(iJ,iterm)))

            ! Gaunt factors for b-f at the edge frequency
            gbfe = gHI_bf(Atom%phot(iphot)%edge,neff,Z)

            ! Calculate the cross section for the relevant frequencies
            do ifreq=if0,if1

              ! Ratio of frequencies (edge and current)
              cfreq = Atom%phot(iphot)%edge/freq(ifreq)

              ! Gaunt factor at this frequency
              gbf  = gHI_bf(freq(ifreq),neff,Z)

              ! Hydrogenic cross section
              Atom%phot(iphot)%alpha(ifreq) = cfreq*cfreq*cfreq* &
                              Atom%phot(iphot)%inalpha(1)*1d4*gbf/gbfe

            end do

          end if ! Photoionization input mode

        end do ! upper levels
      end do ! lower levels

      ! Check if everything is fine
      call gcontrol

      return

      end subroutine setphoto

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates and initializes some photoionization related
      !! quantities.\n
      !!       Atom(Atom_class): Structure with the atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!           T(dfloat(:)): Temperature\n
      !!          ne(dfloat(:)): Electron density\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!          free(logical): If allowed to free cross section
      !!                         data
      subroutine setphotoTEI(Atom,Frec,T,ne,MPID,free)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Frequency_class):: Frec
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: free
      double precision, dimension(:), intent(in):: T,ne

      ! Local

      integer:: ifreq,iphot,if0,if1,iz
      integer:: ilevel,ilevel1

      double precision:: c0,C1,pE,arg,argz
      double precision:: Saha,exu
      double precision, dimension(Rnz):: buffer


      ! 2h/c^2 in k units, convF cancels with the denominator in the
      ! integral and the 1d3 is to take alpha in cm^2
      ! 2d13 = 2d21 (k units) * 1d-8 (s^-1 -> 10^8 s^-1)
      c0 = 2d13*c*4d0*pi
      ! h/K for the argument of the exponential
      arg = c2*1d4

      ! For every pair of levels
      do ilevel=1,Atom%nlevel-1
        do ilevel1=ilevel+1,Atom%nlevel

          ! Check the indexing of photoionizations
          iphot = Atom%iphot(ilevel,ilevel1)

          ! If there is not, cycle to another pair
          if (iphot.lt.1) cycle

          ! Allocations
          ! Integral of spontaneous or LTE recombination
          allocate(Atom%phot(iphot)%TEI(Rz0:Rz1))
          Atom%phot(iphot)%TEI = 0d0


          ! Absent Transition or Master in MPI case running
          ! a non-CLE case
          if (Atom%phot(iphot)%absent.or. &
              (MPID%mpi.and.pid.eq.0.and.run_mode.ne.2)) then

            ! The Master with MPI cannot remove alpha
            if (.not.(pid.eq.0.and.MPID%mpi)) then

              ! Does not care about the value of alpha
              if (free) deallocate(Atom%phot(iphot)%alpha)

            end if

          ! If not absent
          else

            ! Get the transition limits into shorter variables for
            ! convenience
            if0 = Atom%phot(iphot)%if0
            if1 = Atom%phot(iphot)%if1

            !
            ! Calculate the spontaneous or LTE contribution
            !

            ! For each height
            do iz=Rz0,Rz1

              ! Calculate the argument of the exponential
              argz = arg/T(iz)

              ! Calculate the Saha function times electron density
              c1 = fktoJ/kb/T(iz)
              Saha = cSaha*ne(iz)*Atom%phot(iphot)%glu* &
                     exp(Atom%phot(iphot)%edge*c1)/(T(iz)**(1.5d0))

              ! For each non boundary frequency
              do ifreq=if0+1,if1-1

                ! Energy density. One of the freq gets cancel with the
                ! denominator of the integral
                pE = c0*Frec%omega(ifreq)*Frec%omega(ifreq)

                ! Exponential argument
                exu = argz*Frec%omega(ifreq)
                exu = diexp(exu)

                ! Add the contribution to the integral
                Atom%phot(iphot)%TEI(iz) = &
                                     Atom%phot(iphot)%TEI(iz) + &
                                     Frec%W_freq(ifreq)*pE*exu* &
                                     Atom%phot(iphot)%alpha(ifreq)

              end do ! Frequencies

              ! Contribution of the first boundary
              pE = c0*Frec%omega(if0)*Frec%omega(if0)
              exu = argz*Frec%omega(if0)
              exu = diexp(exu)

              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz) + &
                                         Atom%phot(iphot)%W0*pE*exu* &
                                         Atom%phot(iphot)%alpha(if0)

              ! Contribution of the last boundary
              pE = c0*Frec%omega(if1)*Frec%omega(if1)
              exu = argz*Frec%omega(if1)
              exu = diexp(exu)

              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz) + &
                                         Atom%phot(iphot)%W1*pE*exu* &
                                         Atom%phot(iphot)%alpha(if1)

              ! ne*Zeta factor
              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz)* &
                                         Saha

            end do ! Heights

          end if ! If absent line

          ! Collect TEI value if MPI
          if (MPID%mpi) then

            ! Sending buffer
            buffer = Atom%phot(iphot)%TEI

            ! Collect
            call MPI_ALLREDUCE(buffer(1),Atom%phot(iphot)%TEI(Rz0), &
                               Rnz,MPI_DOUBLE_PRECISION,MPI_SUM, &
                               MPI_COMM_RT,ierr)
          end if ! If MPI

        end do ! upper levels
      end do ! lower levels

      ! Check if everything is fine
      call control

      return

      end subroutine setphotoTEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocates variables to accelerate epsIphoto.\n
      !!       Atom(Atom_class): Structure with the atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!           T(dfloat(:)): Temperature\n
      !!        MPID(MPI_class): Structure with MPI data
      subroutine ramphoto(Atom,Frec,T,MPID)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class):: Frec
      type(MPI_class), intent(inout):: MPID
      double precision, dimension(:), intent(in):: T

      ! Local

      logical:: alloc

      integer:: ia,iz,ifreq,ilevel,ilevel1,iphot
      integer:: if0,if1

      double precision:: c0,exu

      ! Limits frequency
      if0 = 100000000
      if1 = -1

      alloc = .False.

      ! For each atom
      do ia=1,nA

        ! For every pair of levels
        do ilevel=1,Atom(ia)%nlevel-1
          do ilevel1=ilevel+1,Atom(ia)%nlevel

            ! Check the indexing of photoionizations
            iphot = Atom(ia)%iphot(ilevel,ilevel1)

            ! If there is not, cycle to another pair
            if (iphot.lt.1) cycle

            ! Update limits
            if(Atom(ia)%phot(iphot)%if0.lt.if0) &
              if0 = Atom(ia)%phot(iphot)%if0
            if(Atom(ia)%phot(iphot)%if1.gt.if1) &
              if1 = Atom(ia)%phot(iphot)%if1

          end do ! Upper levels
        end do ! Lower levels
      end do ! Atoms

      if (if1.ge.if0) alloc = .True.

      if (.not.alloc) then
        if (.not.MPID%mpi) then
          allocate(Frec%exu(1,1))
        else
          nullify(Frec%exu)
        end if
        call control
        return
      end if

      ! Allocate
      allocate(Frec%omega3(if0:if1))
      allocate(Frec%exu(if0:if1,Rz0:Rz1))

      ! Master only needs exu
      if (MPID%mpi.and.pid.eq.0) then

        ! For each frequency
        do ifreq=if0,if1

          ! Compute exponential argument constant
          c0 = c2*1d4*Frec%omega(ifreq)

          ! For each height
          do iz=Rz0,Rz1

            exu = c0/T(iz)
            Frec%exu(ifreq,iz) = diexp(exu)

          end do ! Heights
        end do ! Frequencies

        MPID%PRAM = MPID%PRAM + 8d-6* dble(Rnz*(if1 - if0 + 1))

      ! Slaves need both
      else

        ! For each frequency
        do ifreq=if0,if1


          ! Compute frequency cube
          Frec%omega3(ifreq) = Frec%omega(ifreq)*Frec%omega(ifreq)* &
                               Frec%omega(ifreq)

          ! Compute exponential argument constant
          c0 = c2*1d4*Frec%omega(ifreq)

          ! For each height
          do iz=Rz0,Rz1

            exu = c0/T(iz)
            Frec%exu(ifreq,iz) = diexp(exu)

          end do ! Heights
        end do ! Frequencies

        MPID%PRAM = MPID%PRAM + 8d-6* dble(Rnz*(if1 - if0 + 1))

      end if ! Master/slave

      ! Check if everything is fine
      call control

      return

      end subroutine ramphoto

!#####################################################################
!#####################################################################
!#####################################################################

      end module initphotoion_mod

