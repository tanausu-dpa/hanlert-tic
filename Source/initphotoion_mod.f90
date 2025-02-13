      !> Initialization of photoionization quantities
      module initphotoion_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     11/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  setphoto
!    Allocate and initialize photoionization related quantities that
!  are independent of the model atmosphere
!
!  setphotoTEI
!    Calculate the photoionization thermal contribution that goes into
!  the SEE
!
!  ramphoto
!    Allocate and calculate frequency quantities used in photoeps
!  and photoepsI routines
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

      !> Allocate and initialize photoionization related quantities
      !! that are independent of the model atmosphere\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!     freq(double(:)): Frequency array
      subroutine setphoto(Atom,freq)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      double precision, dimension(:), intent(in):: freq

      ! Local

      logical:: llinear

      integer:: nfr,ifreq,iphot,if0,if1,iterm,iterm1
      integer:: ilevel,ilevel1,iJ,iJ1

      double precision:: Z,neff,gbf,gbfe,cfreq
      double precision, dimension(:), allocatable:: sb,sc,sd,sx,sy


      ! For every lower level
      do ilevel=1,Atom%nlevel-1

        ! For every upper level
        do ilevel1=ilevel+1,Atom%nlevel

          ! Check the index of photoionization
          iphot = Atom%iphot(ilevel,ilevel1)

          ! If there is not one, cycle to another pair
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

          ! If Master and doing MPI, store
          if (pid.eq.0.and.nproc.gt.1) then
            Atom%phot(iphot)%MW0 = Atom%phot(iphot)%W0
            Atom%phot(iphot)%MW1 = Atom%phot(iphot)%W1
          end if

          ! Allocate cross-section
          allocate(Atom%phot(iphot)%alpha(if0:if1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(iphot)%alpha)

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

            ! Get x and y input axes
            sx(nfr:1:-1) = 1d2/Atom%phot(iphot)%infreq
            sy(nfr:1:-1) = Atom%phot(iphot)%inalpha

            ! Get the Spline coefficients
            call spline(sx(1:nfr),sy(1:nfr), &
                        sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)

            ! For the relevant frequencies
            do ifreq=if0,if1

              ! Calculate the cross section in cm^2
              Atom%phot(iphot)%alpha(ifreq) = &
                ispline(1d2/freq(ifreq),sx(1:nfr),sy(1:nfr), &
                        sb(1:nfr),sc(1:nfr),sd(1:nfr),nfr)*1d4

              ! If something negative
              if (Atom%phot(iphot)%alpha(ifreq).lt.0d0) then

                ! Master
                if (pid.eq.0) then

                  ! Issue warning
                  write(umsg,'(A,1x,i4,A,1x,i4,1x,A,A,A)') &
                             ' # Spline interpolation of '// &
                             'ionization cross section of '// &
                             'b-f transition',ilevel,' -->', &
                             ilevel1,' of ', &
                             Atom%Element,' atom gave '// &
                             'negative value, doing linear.'
                  call verbose

                end if ! Master

                ! Flag to do linear and leave
                llinear = .True.
                exit

              end if ! Negative cross-section

            end do ! Relevant frequencies

            ! If need to do it linear
            if (llinear) then

              ! For the relevant frequencies
              do ifreq=if0,if1

                ! Calculate the cross section
                call linear(sx(1:nfr),sy(1:nfr), &
                            1d2/freq(ifreq), &
                            Atom%phot(iphot)%alpha(ifreq))

                ! Convert units to cm^2
                Atom%phot(iphot)%alpha(ifreq) = &
                                     Atom%phot(iphot)%alpha(ifreq)*1d4

              end do ! Relevant frequencies

            end if ! Need to interpolate linearly

          ! If the input mode is Hydrogenic
          else

            ! Find the term and J index of the levels involved
            iterm1 = Atom%term(ilevel1)
            iJ1 = Atom%sublevel(ilevel1)
            iterm = Atom%term(ilevel)
            iJ = Atom%sublevel(ilevel)

            ! Charge
            Z = dble(Atom%stage(iterm1) - 1)

            ! Effective principal quantum number
            neff = Z*sqrt(ryd/ &
                   (Atom%FSfreq(iJ1,iterm1) - Atom%FSfreq(iJ,iterm)))

            ! Gaunt factors for b-f at the edge frequency
            gbfe = gHI_bf(Atom%phot(iphot)%edge,neff,Z)

            ! For the relevant frequencies
            do ifreq=if0,if1

              ! Ratio of frequencies (edge and current)
              cfreq = Atom%phot(iphot)%edge/freq(ifreq)

              ! Gaunt factor at this frequency
              gbf  = gHI_bf(freq(ifreq),neff,Z)

              ! Hydrogenic cross-section in cm^2
              Atom%phot(iphot)%alpha(ifreq) = cfreq*cfreq*cfreq* &
                              Atom%phot(iphot)%inalpha(1)*1d4*gbf/gbfe

            end do ! Relevant frequencies

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

      !> Calculate the photoionization thermal contribution that goes
      !! into the SEE\n
      !!       Atom(Atom_class): Structure with atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!           T(double(:)): Temperature\n
      !!          ne(double(:)): Electron number density\n
      !!          free(logical): If allowed to free cross-section data
      subroutine setphotoTEI(Atom,Frec,T,ne,free)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Frequency_class), intent(in):: Frec
      logical, intent(in):: free
      double precision, dimension(:), intent(in):: T,ne

      ! Local

      integer:: ifreq,iphot,if0,if1,iz,ilevel,ilevel1

      double precision:: c0,C1,pE,arg,argz,Saha,exu
      double precision, dimension(Rnz):: buffer


      ! 2h/c^2 in k units, convF cancels with the denominator in the
      ! integral and the 1d3 is to take alpha in cm^2
      ! 2d13 = 2d21 (k units) * 1d-8 (s^-1 -> 10^8 s^-1)
      c0 = 2d13*c*4d0*pi
      ! h/K for the argument of the exponential
      arg = c2*1d4

      ! For every lower level
      do ilevel=1,Atom%nlevel-1

        ! For every upper level
        do ilevel1=ilevel+1,Atom%nlevel

          ! Check the index of the photoionization
          iphot = Atom%iphot(ilevel,ilevel1)

          ! If there is not one, cycle to another pair
          if (iphot.lt.1) cycle

          ! Allocate integral of spontaneous or LTE recombination
          allocate(Atom%phot(iphot)%TEI(Rz0:Rz1))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%phot(iphot)%TEI)

          ! Initialize
          Atom%phot(iphot)%TEI = 0d0

          ! If absent transition or Master in MPI case running
          ! a non-CLE case
          if (Atom%phot(iphot)%absent.or. &
              (nproc.gt.1.and.pid.eq.0.and.run_mode.ne.2)) then

            ! The Master with MPI cannot remove alpha
            if (.not.(pid.eq.0.and.nproc.gt.1)) then

              ! Does not care about the value of alpha and can free it
              if (free) then

                ! Free memory
                MRAMc = MRAMc - 1d-6*sizeof(Atom%phot(iphot)%alpha)
                deallocate(Atom%phot(iphot)%alpha)

              end if ! Can free memory
            end if ! Slave

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

              ! Add contribution
              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz) + &
                                         Atom%phot(iphot)%W0*pE*exu* &
                                         Atom%phot(iphot)%alpha(if0)

              ! Contribution of the last boundary
              pE = c0*Frec%omega(if1)*Frec%omega(if1)
              exu = argz*Frec%omega(if1)
              exu = diexp(exu)

              ! Add contribution
              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz) + &
                                         Atom%phot(iphot)%W1*pE*exu* &
                                         Atom%phot(iphot)%alpha(if1)

              ! ne*Zeta factor
              Atom%phot(iphot)%TEI(iz) = Atom%phot(iphot)%TEI(iz)* &
                                         Saha

            end do ! Heights

          end if ! Absent/present ionization

          ! Collect TEI value if MPI
          if (nproc.gt.1) then

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

      !> Allocate and calculate frequency quantities used in photoeps
      !! and photoepsI routines\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!           T(double(:)): Temperature
      subroutine ramphoto(Atom,Frec,T)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(inout):: Frec
      double precision, dimension(:), intent(in):: T

      ! Local

      logical:: alloc

      integer:: ia,iz,ifreq,ilevel,ilevel1,iphot,if0,if1

      double precision:: c0,exu


      ! Frequency limits
      if0 = 100000000
      if1 = -1

      ! Initialize if allocating
      alloc = .False.

      ! For each atom
      do ia=1,nA

        ! For every lower level
        do ilevel=1,Atom(ia)%nlevel-1

          ! For every ipper level
          do ilevel1=ilevel+1,Atom(ia)%nlevel

            ! Check the index of photoionization
            iphot = Atom(ia)%iphot(ilevel,ilevel1)

            ! If there is not one, cycle to another pair
            if (iphot.lt.1) cycle

            ! Update frequency limits
            if(Atom(ia)%phot(iphot)%if0.lt.if0) &
              if0 = Atom(ia)%phot(iphot)%if0
            if(Atom(ia)%phot(iphot)%if1.gt.if1) &
              if1 = Atom(ia)%phot(iphot)%if1

          end do ! Upper levels
        end do ! Lower levels
      end do ! Atoms

      ! If valid limits, allocate
      if (if1.ge.if0) alloc = .True.

      ! If cannot allocate
      if (.not.alloc) then

        ! Serial
        if (nproc.le.1) then

          ! Dummy allocation
          allocate(Frec%exu(1,1))
          PRAMc = PRAMc + 1d-6*sizeof(Frec%exu)

        ! MPI
        else

          ! Nullify pointer
          nullify(Frec%exu)

        end if ! MPI/serial

        ! Check
        call control
        return

      end if ! Cannot allocate

      ! Allocate
      allocate(Frec%exu(if0:if1,Rz0:Rz1))
      PRAMc = PRAMc + 1d-6*sizeof(Frec%exu)

      ! Master in MPI only needs exu
      if (nproc.gt.1.and.pid.eq.0) then

        ! For each frequency
        do ifreq=if0,if1

          ! Compute exponential argument constant
          c0 = c2*1d4*Frec%omega(ifreq)

          ! For each height
          do iz=Rz0,Rz1

            ! Compute inverse exponential
            exu = c0/T(iz)
            Frec%exu(ifreq,iz) = diexp(exu)

          end do ! Heights
        end do ! Frequencies

      ! Slaves need both
      else

        ! Allocate
        allocate(Frec%omega3(if0:if1))
        PRAMc = PRAMc + 1d-6*sizeof(Frec%omega3)

        ! For each frequency
        do ifreq=if0,if1

          ! Compute frequency cube
          Frec%omega3(ifreq) = Frec%omega(ifreq)*Frec%omega(ifreq)* &
                               Frec%omega(ifreq)

          ! Compute exponential argument constant
          c0 = c2*1d4*Frec%omega(ifreq)

          ! For each height
          do iz=Rz0,Rz1

            ! Compute inverse exponential
            exu = c0/T(iz)
            Frec%exu(ifreq,iz) = diexp(exu)

          end do ! Heights
        end do ! Frequencies

      end if ! Master/slave

      ! Check if everything is fine
      call control

      return

      end subroutine ramphoto

!#####################################################################
!#####################################################################
!#####################################################################

      end module initphotoion_mod

