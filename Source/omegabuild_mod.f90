      !> Routines to generate frequency axes
      module omegabuild_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Contributors:
!     John Dennis (NCAR)
!  Start:
!     18/04/2017
!  Last version:
!     15/06/2026 V4.0.8
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/06/2026:    V4.0.8 - Added the definition of freqR and freqB
!                             variables in Atom_class for hydrogen and
!                             helium (TdPA)
!                           - Bugfix: The wrong variable was being
!                             checked in order to identify the
!                             limits of the loop over transitions in
!                             find_integral_limits (TdPA)
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
!  omegabuild
!    Build the frequency axis of the observer
!
!  omegainitmaster:
!    Initialize frequency related variables needed by the radiative
!  transfer master
!
!  get_transition_out_limit
!    Determine the frequency limits for an output transition in the
!  second order emissivity
!
!  get_input_frequencies
!    Determine the input frequency axis given the transition
!  resonances and the discretization parameters of a pair of output
!  and input transitions
!
!  setmpi_red
!    Split tasks for the calculation of the second order emissivity of
!  a given transition
!
!  find_integral_limits
!    Find the frequency index limits for the calculation of the second
!  order emissivity once the input frequency axes are known
!
!  safe_allocate
!    Allocate a double precision array. If the array is allocated,
!  free the space unless it has the correct dimension already
!
!  omegabuildin
!    Determine the input frequency axis for the calculation of the
!  second order emissivity and split the tasks of the emissivity
!  calculation as evenly as possible. Assume comoving reference
!  frame
!
!  allocate_Warr
!    Allocate space to store the redistribution functions
!
!  omegabuildinI
!    Determine the input frequency axis for the calculation of the
!  second order intensity emissivity and split the tasks of the
!  intensity emissivity calculation as evenly as possible. Assume
!  comoving reference frame
!
!  allocate_WarrI
!    Allocate space to store the intensity redistribution functions
!
!  freqresize:
!    Resize some frequency dependent quantities and adjust indexes
!  for each CPU taking into account the range of frequencies they
!  need to take care of
!
!  refitfrec
!    Resize some frequency dependent quantities and adjust indexes
!  for each CPU taking into account the range of frequencies they
!  need to take care of, and create the output frequency axis, for
!  a CLE synthesis
!
!  index_norm
!    Index and allocate arrays for the normalization data and
!  estimate the minimum RAM neccesary
!
!  index_red
!    Index array for the redistribution quantities and allocate
!  the structures to hold the input frequency axis
!
!  check_nchlt
!    Checks where the non-coherent lower term approximation can be
!  applied
!
!  cleanFrecandRed
!    Safely clean the Frec and Red structures
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use parameters_mod , only : c , TINYA, TINYB, TINYO, TINYWAR , &
                                  resol, resolin , pi , IPI42, jump, &
                                  TINYVEL
      use profile_mod
      use qsort_mod
      use setmpi_mod
      use types_mod

      ! Parameters
      integer, parameter:: ContW = 1
      integer, parameter:: PhotW = 2
      integer, parameter:: LLTEW = 5
      integer, parameter:: LCRDW = 10


      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Build the frequency axis of the observer\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!     Input(Input_class): Structure with configuration data\n
      !!           maxB(double): Maximum magnetic field strength\n
      !!    obs_wave(double(:)): In inversion mode, the data
      !!                         wavelengths, and a dummy array
      !!                         otherwise
      subroutine omegabuild(Frec,Atom,Input,maxB,obs_wave)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Input_class), intent(inout):: Input
      type(Frequency_class), intent(inout):: Frec
      double precision, intent(in):: maxB
      double precision, dimension(:), &
                        allocatable, intent(in):: obs_wave

      ! Local

      logical:: init,warned,MIT
      logical, dimension(:), allocatable:: protect,vMIT

      integer:: ia,itran,jtran,ktran,iterml,itermu
      integer:: ifreq,jfreq,cfreq,it,lnfreq,lnfreqc
      integer:: i,iJl,iJu,nt,ios,nMl,nMu
      integer, dimension(:), allocatable:: flag

      double precision:: norm,norm1,v,dv,O0,O1,maxV,disp,vfacl,vfacr
      double precision, dimension(:), allocatable:: omegaaux
      double precision, dimension(:), allocatable:: omega
      double precision, dimension(:), allocatable:: nut
      double precision, dimension(:), allocatable:: DwT
      double precision, dimension(:), allocatable:: DwTL
      double precision, dimension(:), allocatable:: DwTmin
      double precision, dimension(:), allocatable:: DwTmax
      double precision, dimension(:,:), allocatable:: tmp_lim
      double precision, dimension(:,:), allocatable:: tmp_limL


      ! Routine name
      urou = 'omegabuild'

      !
      ! Read the number of wavelengths from the wavelength files
      ! and add it to the total frequencies
      !

      ! For each wavelength file
      do i=1,Input%NW

        ! Open
        open (200,file=trim(Input%waves(i)%str), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', form='unformatted')

        ! Read size
        read(200,err=1100) jfreq

        ! Add size to the total
        nfreq = nfreq + jfreq

        ! Close
        close(200)

      end do ! For each wavelength file

      ! If forcing inversion data frequencies
      if (Input%force_inv_Freq) then

        ! Add size to the total
        nfreq = nfreq + size(obs_wave)

      end if

      ! Maximum and minimum Doppler widths
      vfacl = 1d0 - Input%maxV
      vfacr = 1d0 + Input%maxV

      ! Compute minimum and maximum Doppler width
      allocate(DwTmax(nA),DwTmin(nA))
      do ia=1,nA
        DwTmax(ia) = Atom(ia)%cDopp*sqrt(Input%maxT)
        DwTmin(ia) = Atom(ia)%cDopp*sqrt(Input%minT)
      end do

      ! If active atoms
      if (nA.gt.0) then
          
        ! Take the Doppler width from the input
        allocate(DwT(nA))

        ! Maximum Doppler width
        if(Input%dws.eq.'MAX')then

          ! Copy
          DwT = DwTmax

        ! Minimum Doppler width
        else if(Input%dws.eq.'MIN')then

          ! Copy
          DwT = DwTmin

        ! Or fixed input Doppler width
        else if(Input%dws.eq.'NUM')then

          ! Compute
          DwT = Input%dw*1d-9/c

        end if
      end if ! Active atoms

      ! If LTE lines
      if (nLTEl.gt.0) then

        ! Take the Doppler width from the input
        allocate(DwTL(nLTEl))

        ! Maximum Doppler width
        if(Input%dws.eq.'MAX')then

          ! For every line
          do ia=1,nLTEl
            DwTL(ia) = Input%LTEline(ia)%cDopp*sqrt(Input%maxT)
          end do

        ! Minimum Doppler width
        else if(Input%dws.eq.'MIN')then

          ! For every line
          do ia=1,nLTEl
            DwTL(ia) = Input%LTEline(ia)%cDopp*sqrt(Input%minT)
          end do

        ! Or fixed input Doppler width
        else if(Input%dws.eq.'NUM')then

          ! Compute
          DwTL = Input%dw*1d-9/c

        end if
      end if ! LTE lines

      ! Correct nfreq if including MIT transitions
      if (Input%MIT_input.ge.0) then

        ! For each atom
        do ia=1,nA

          !
          ! Check if need to include MIT transitions
          !

          ! If multi-level
          if (Atom(ia)%ML) then

            ! No MIT by definition
            MIT = .False.

          ! If multi-term
          else

            ! If forced to consider MIT lines
            if (Input%MIT_input.gt.0) then

              ! Consider MIT
              MIT = .True.

            ! Forced to ignore MIT lines
            else if (Input%MIT_input.lt.0) then

              ! Ignore MIT
              MIT = .False.

            ! If free to choose
            else

              ! Include if magnetic field
              MIT = maxB.gt.TINYB

            end if ! MIT input
          end if ! ML vs MT atom

          ! If not MIT, skip this loop
          if (.not.MIT) cycle

          ! Skipping wavelengths?
          if (Input%skip_wave(ia)) cycle

          ! bound-bound transitions
          do itran=1,Atom(ia)%ntran

            ! If not splitting between components, skip
            if (.not.Atom(ia)%splitf(itran)) cycle

            ! Get terms
            itermu = Atom(ia)%fst(itran)%itermu
            iterml = Atom(ia)%fst(itran)%iterml

            ! Loop fine transitions
            do iJu=1,Atom(ia)%nJ(itermu)
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If transition exist, skip
                if (Atom(ia)%fst(itran)%irad(iJu,iJl).ge.1) cycle

                ! If MIT
                if (MIT) then

                  ! Total number of frequencies
                  lnfreq = nint(dble(Atom(ia)%NfreqT(itran))* &
                                Input%MIT_node)

                  ! Ensure odd
                  if (mod(lnfreq,2).eq.0) lnfreq = lnfreq + 1

                  ! Core number of frequencies
                  lnfreqc = nint(dble(Atom(ia)%NfreqTc(itran))* &
                                 Input%MIT_node)

                  ! Ensure odd
                  if (mod(lnfreqc,2).eq.0) lnfreqc = lnfreqc + 1

                  ! Ensure enough frequencies to make a line
                  if (lnfreqc.lt.0) lnfreqc = 3
                  if (lnfreq.lt.lnfreqc) lnfreq = lnfreqc + 2

                  ! Add frequencies total
                  nfreq = nfreq + lnfreq
                  Atom(ia)%nfreq = Atom(ia)%nfreq + lnfreq

                end if ! If including MIT

              end do ! Lower J level
            end do ! Upper J level
          end do ! L-L transition
        end do ! Atom

      end if ! If possible MIT transitions

      ! Allocate temporal limits array
      if (nxtran.gt.0) allocate(tmp_lim(2,nxtran))
      if (nLTEl.gt.0) allocate(tmp_limL(2,nLTEl))

      ! Set counter of frequencies to 0
      cfreq = 0

      !
      ! Get maximum number of FS transitions
      !

      ! Initialize counter
      nMu = 0

      ! For every atom
      do ia=1,nA

        ! Get maximum number of J components
        nMl = maxval(Atom(ia)%nJ)
        if (nMl.gt.nMu) nMu = nMl

      end do

      ! Absolute maximum is the square
      nMu = nMu*nMu


      !
      ! Allocations
      !

      ! Initial vector
      allocate(omega(nfreq))
      ! Fine structure frequencies
      if (nMu.gt.0) allocate(nut(nMu))
      ! Vector for MIT flag
      if (nMu.gt.0) allocate(vMIT(nMU))

      ! Initialize quantities to check velocity shifts
      if (dyn) maxV = Input%maxV

      ! Initialize flag of existing warning
      warned = .False.


      !
      ! Collect the frequencies
      !

      ! For each atom
      do ia=1,nA

        ! Make sure that there is space for atomic frequencies
        if (.not.allocated(omegaaux))then
          allocate(omegaaux(Atom(ia)%nfreq))
        else
          if(size(omegaaux).lt.Atom(ia)%nfreq)then
            deallocate(omegaaux)
            allocate(omegaaux(Atom(ia)%nfreq))
          end if
        end if

        !
        ! Check if need to include MIT transitions
        !

        ! If multi-level
        if (Atom(ia)%ML) then

          ! Do not include
          MIT = .False.

        ! Multi-term
        else

          ! Forced to include them
          if (Input%MIT_input.gt.0) then

            ! Include
            MIT = .True.

          ! Forced to neglect them
          else if (Input%MIT_input.lt.0) then

            ! Neglect them
            MIT = .False.

          ! Free to determine
          else

            ! Include if magnetic field
            MIT = maxB.gt.TINYB

          end if ! Input about MIT
        end if ! ML vs MT

        ! Bound-bound transitions
        do itran=1,Atom(ia)%ntran

          ! Apply atom shift
          jtran = itran + Atom(ia)%tshift

          ! Initialize the limits of the line
          tmp_lim(1,jtran) =  1D99
          tmp_lim(2,jtran) = -1D99

          ! Get terms
          itermu = Atom(ia)%fst(itran)%itermu
          iterml = Atom(ia)%fst(itran)%iterml

          ! If splitting between components
          if (Atom(ia)%splitf(itran)) then

            !
            ! Count fine transitions
            !

            ! Initialize count
            nt = 0

            ! For each level in upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If transition does not exist
                if (Atom(ia)%fst(itran)%irad(iJu,iJl).lt.1) then

                  ! If including MIT
                  if (MIT) then

                    ! Add it anyways
                    nt = nt + 1
                    vMIT(nt) = .True.

                  ! If not MIT
                  else

                    ! Just cycle
                    cycle

                  end if

                ! Permitted transition
                else

                  ! Add to count
                  nt = nt + 1
                  vMIT(nt) = .False.

                end if ! Allowed transition

                ! Get resonance
                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)

              end do ! Lower level in term
            end do ! Upper level in term

          ! Not splitting between components
          else

            ! Single big transition
            nt = 1
            nut(1) = Atom(ia)%Dfreq(itran)
            vMIT(1) = .False.

          end if

          !
          ! Add frequencies for each FS transition
          !

          ! For each FS transition
          do it=1,nt

            ! If MIT transition
            if (vMIT(it)) then

              ! Total number of frequencies
              lnfreq = nint(dble(Atom(ia)%NfreqT(itran))* &
                            Input%MIT_node)

              ! Ensure odd
              if (mod(lnfreq,2).eq.0) lnfreq = lnfreq + 1

              ! Core number of frequencies
              lnfreqc = nint(dble(Atom(ia)%NfreqTc(itran))* &
                             Input%MIT_node)

              ! Ensure odd
              if (mod(lnfreqc,2).eq.0) lnfreqc = lnfreqc + 1

              ! Ensure minimum amount of frequencies
              if (lnfreqc.lt.0) lnfreqc = 3
              if (lnfreq.lt.lnfreqc) lnfreq = lnfreqc + 2

            ! If not MIT
            else

              ! Total number of frequencies
              lnfreq = Atom(ia)%NfreqT(itran)

              ! Core number of frequencies
              lnfreqc = Atom(ia)%NfreqTc(itran)

            end if ! MIT or not

            !
            ! Linear part of the axis
            !

            ! Initialize center
            omegaaux(1) = 0d0

            ! To the right in the core
            do ifreq=2,lnfreqc/2+1

              ! Next frequency
              omegaaux(ifreq) = omegaaux(ifreq-1) + &
                                Atom(ia)%Dwvlc(itran)/ &
                                dble(lnfreqc/2)
            end do

            !
            ! Logarithmic part
            !

            ! Initial and step
            v = log10(Atom(ia)%Dwvlc(itran))
            dv = log10(Atom(ia)%Dwvl(itran)) - v

            ! For wing frequencies to the right
            do ifreq=lnfreqc/2+2,lnfreq/2+1

              ! Next frequency
              v = v + dV/dble(lnfreq/2 - lnfreqc/2)
              omegaaux(ifreq) = 1d1**v

            end do

            ! Build symmetric axis
            omegaaux(lnfreq/2+2:lnfreq) = omegaaux(2:lnfreq/2+1)
            omegaaux(1:lnfreq/2) = -omegaaux(lnfreq/2+1:2:-1)
            omegaaux(lnfreq/2+1) = 0d0

            ! Check limits in actual frequency
            O0 = minval(omegaaux(1:lnfreq))* &
                 Atom(ia)%Dfreq(itran)*DwT(ia) + nut(it)
            O1 = maxval(omegaaux(1:lnfreq))* &
                 Atom(ia)%Dfreq(itran)*DwT(ia) + nut(it)

            ! Update limits for term-term transition
            if (O0.lt.tmp_lim(1,jtran)) tmp_lim(1,jtran) = O0
            if (O1.gt.tmp_lim(2,jtran)) tmp_lim(2,jtran) = O1

            ! Skip wavelengths for this atom?
            if (Input%skip_wave(ia)) cycle

            ! Build real contribution and add to total axis
            do ifreq=1,lnfreq

              ! Add to frequency axis
              omega(ifreq + cfreq) = nut(it) + &
                                     omegaaux(ifreq)* &
                                     Atom(ia)%Dfreq(itran)*DwT(ia)

            end do ! Atomic frequencies

            ! Accumulate the defined frequencies
            cfreq = cfreq + lnfreq

          end do ! FS transition

          !
          ! Check that the range is adecuated to the velocity imposed
          ! in this frame
          !

          ! If dynamic and master
          if (dyn.and.gpid.eq.0) then

            ! For each transition
            do it=1,nt

            ! Compute maximum Doppler shift
              disp = maxV*nut(it)/DwT(ia)/Atom(ia)%Dfreq(itran)

              ! If goes to half of the range
              if (2d0*disp.gt.Atom(ia)%Dwvl(itran)) then

                ! No warning yet
                if (.not.warned) then

                  ! Warning header
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose

                  ! Do only once
                  warned = .True.

                end if

                ! Issue warning
                write(umsg,'(A,i4,",",i4,3A,2(f6.1,A))') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than half of the total width specified (', &
                  disp,'>',Atom(ia)%Dwvl(itran),'/2)'
                call verbose

              ! Fully out of core
              else if (2d0*disp.gt.Atom(ia)%Dwvlc(itran)) then

                ! No warning yet
                if (.not.warned) then

                  !i Warning header
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose

                  ! Do only once
                  warned = .True.

                end if

                ! Issue warning
                write(umsg,'(A,i4,",",i4,3A,2(f6.1,A))') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than half of the core width specified (', &
                  disp,'>',Atom(ia)%Dwvlc(itran),'/2)'
                call verbose

              ! If goes out of the core region
              else if (disp.gt.Atom(ia)%Dwvlc(itran)) then

                ! No warning yet
                if (.not.warned) then

                  ! Warning header
                  umsg = '###'
                  call verbose
                  umsg = '### IMPORTANT WARNING'
                  call verbose
                  umsg = '###'
                  call verbose

                  ! Do only once
                  warned = .True.

                end if

                ! Issue warning
                write(umsg,'(A,i4,",",i4,3A,2(f6.1,A))') &
                  ' - Warning: transition ',itran,it,' in ', &
                  Atom(ia)%Element,' atom can be shifted more '// &
                  'than the core width specified (', &
                  disp,'>',Atom(ia)%Dwvlc(itran),')'
                call verbose

              end if ! Big displacements

            end do ! Line components

          end if ! Dynamic and master

        end do ! b-b Transition

        ! Bound-free transitions
        do itran=1,Atom(ia)%nphot

          ! Skip this atom?
          if (Input%skip_wave(ia)) cycle

          ! If it is explicit
          if (Atom(ia)%phot(itran)%mode.eq.0) then

            ! Just add frequencies to the axis
            do ifreq=1,Atom(ia)%phot(itran)%nfreq
              omega(ifreq+cfreq) = Atom(ia)%phot(itran)%infreq(ifreq)
            end do

          ! If it is hydrogenic
          else

            ! Start from the edge
            omega(1 + cfreq) = Atom(ia)%phot(itran)%edge

            ! Get index for maximum frequency
            nt = (1 - Atom(ia)%phot(itran)%mode)* &
                 Atom(ia)%phot(itran)%nfreq + &
                 Atom(ia)%phot(itran)%mode

            ! Get step
            dv = (Atom(ia)%phot(itran)%infreq(nt) - omega(1+cfreq))/ &
                 dble(Atom(ia)%phot(itran)%nfreq - 1)

            ! For rest of photoionization frequencies
            do ifreq=2,Atom(ia)%phot(itran)%nfreq

              ! Linear steps
              omega(ifreq + cfreq) = omega(ifreq - 1 + cfreq) + dv

            end do ! Frequencies

          end if ! Explicit or hydrogenic

          ! Accumulate the defined frequencies
          cfreq = cfreq + Atom(ia)%phot(itran)%nfreq

        end do ! b-f Transiton
      end do ! Atom

      !
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If no frequencies, do not bother
        if (Input%LTEline(ia)%nfreq.le.0) cycle

        ! Make sure that there is memory to work with
        if(.not.allocated(omegaaux))then
          allocate(omegaaux(Input%LTEline(ia)%nfreq))
        else
          if(size(omegaaux).lt.Input%LTEline(ia)%nfreq)then
            deallocate(omegaaux)
            allocate(omegaaux(Input%LTEline(ia)%nfreq))
          end if
        end if

        ! Initialize the limits of the line
        tmp_limL(1,ia) =  1D99
        tmp_limL(2,ia) = -1D99

        !
        ! Linear part of the axis
        !

        ! Start with the center
        omegaaux(1) = 0d0

        ! For the rest of core frequencies to the right
        do ifreq=2,Input%LTEline(ia)%nfreqc/2+1

          ! Linear step
          omegaaux(ifreq) = omegaaux(ifreq-1) + &
                            Input%LTEline(ia)%Dwvlc/ &
                            dble(Input%LTEline(ia)%nfreqc/2)
        end do

        !
        ! Logarithmic part
        !

        ! Get initial and step
        v = log10(Input%LTEline(ia)%Dwvlc)
        dv = log10(Input%LTEline(ia)%Dwvl) - v

        ! For the wing frequencies
        do ifreq=Input%LTEline(ia)%nfreqc/2+2, &
                 Input%LTEline(ia)%nfreq/2+1

          ! Logarithmic steps
          v = v + dV/dble(Input%LTEline(ia)%nfreq/2 - &
                          Input%LTEline(ia)%nfreqc/2)
          omegaaux(ifreq) = 1d1**v

        end do

        ! Build symmetric axis
        omegaaux(Input%LTEline(ia)%nfreq/2+2: &
                 Input%LTEline(ia)%nfreq) = &
                               omegaaux(2:Input%LTEline(ia)%nfreq/2+1)
        omegaaux(1:Input%LTEline(ia)%nfreq/2) = &
                           -omegaaux(Input%LTEline(ia)%nfreq/2+1:2:-1)
        omegaaux(Input%LTEline(ia)%nfreq/2+1) = 0d0

        ! Check limits
        tmp_limL(1,ia) = minval(omegaaux(1:Input%LTEline(ia)%nfreq))*&
                         Input%LTEline(ia)%Dfreq*DwTL(ia) + &
                         Input%LTEline(ia)%Dfreq
        tmp_limL(2,ia) = maxval(omegaaux(1:Input%LTEline(ia)%nfreq))*&
                           Input%LTEline(ia)%Dfreq*DwTL(ia) + &
                           Input%LTEline(ia)%Dfreq

        ! Skip if this line does not add wavelengths
        if (Input%LTEline(ia)%nowave) cycle

        !
        ! Build real contribution
        !

        ! For all line frequencies
        do ifreq=1,Input%LTEline(ia)%nfreq

          ! Add to frequency axis
          omega(ifreq + cfreq) = Input%LTEline(ia)%Dfreq + &
                                 omegaaux(ifreq)* &
                                 Input%LTEline(ia)%Dfreq*DwTL(ia)

        end do ! Atomic frequencies

        ! Accumulate the defined frequencies
        cfreq = cfreq + Input%LTEline(ia)%nfreq

      end do ! LTE lines


      !
      ! Wavelength files
      !

      ! For each wavelength file
      do i=1,Input%NW

        ! Open file
        open (200,file=trim(Input%waves(i)%str), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', form='unformatted')

        ! Read size of array
        read(200,err=1100) jfreq

        ! Read wavelengths
        read(200,err=1100) omega(cfreq+1:cfreq+jfreq)

        ! Convert
        omega(cfreq+1:cfreq+jfreq) = 1d2/omega(cfreq+1:cfreq+jfreq)

        ! Update last index
        cfreq = cfreq + jfreq

        ! Close file
        close(200)

      end do ! Wavelength files

      !
      ! Inversion data wavelengths
      !

      ! If forcing observed frequencies
      if (Input%force_inv_Freq) then

        ! Add frequencies
        omega(cfreq+1:cfreq+size(obs_wave)) = 1d2/obs_wave

        ! For each data wavelengths
        do i=1,size(obs_wave)

          ! Check if already in
          do ifreq=1,cfreq

            ! If a close frequency exist, negate it
            if(abs(1d2/omega(ifreq)-obs_wave(i)).lt.resol) &
              omega(ifreq) = -1d0

          end do ! Synthesis freqs.
        end do ! Data freqs.

        ! Update size
        cfreq = cfreq + size(obs_wave)

      end if ! Force inversion wavelengths


      !
      ! Check for duplicates
      !

      ! Allocate a flag of valid frequencies
      allocate(flag(cfreq))
      flag = 1

      ! For each frequency
      do ifreq=1,cfreq

        ! If it has been flagged, we already checked
        if (flag(ifreq).lt.1) cycle

        ! If nosense
        if (omega(ifreq).le.0d0) then

          ! Flag and continue
          flag(ifreq) = 0
          cycle

        end if

        ! Check the following ones
        do jfreq=ifreq+1,cfreq

          ! If it has been flagged, we already checked
          if (flag(jfreq).lt.1) cycle

          ! If some of them are repeated, flag them to be removed
          if (abs(1d2/omega(ifreq)-1d2/omega(jfreq)).lt.resol) &
            flag(jfreq) = 0

        end do ! jfreq
      end do ! ifreq

      ! Reset the running real index
      jfreq = 0

      ! For each frequency in the vector
      do ifreq=1,cfreq

        ! If it is flagged correct
        if(flag(ifreq).gt..5)then

          ! Add to real vector
          jfreq = jfreq + 1
          omega(jfreq) = omega(ifreq)

        end if ! Correct

      end do ! All frequencies

      ! The number of frequencies is the number of admitted
      ! frequencies in the flag vector
      cfreq = sum(flag)
      nfreq = cfreq

      ! Deallocate the flag
      deallocate(flag)

      ! At this points we need sanity check
      if (nfreq.lt.1) then

        ! No frequencies!
        urou = 'omegabuild'
        umsg = 'Error when creating wavelength axis, there '// &
               'are no wavelengths to allocate'
        call abortedS(umsg,urou,.True.,.True.)
        call control
        return

      end if


      !
      ! Allocate the true axes
      !

      ! Frequencies
      allocate(Frec%omega(nfreq))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%omega)
      ! Integration weights
      allocate(Frec%W_freq(nfreq))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%W_freq)
      ! Node weights
      allocate(Frec%IW_freq(nfreq))
      Frec%IW_freq = ContW
      ! If CLE
      if (run_mode.eq.2) then
        ! Input weights
        allocate(Frec%IW_freq_in(nfreq))
        Frec%IW_freq_in = 0
      end if

      ! Take only the valid frequencies and deallocate the auxiliar
      Frec%omega = omega(1:nfreq)
      deallocate(omega)

      ! Free
      if (allocated(nut)) deallocate(nut,vMIT)

      ! Order the frequencies in the axis
      call QsortC(Frec%omega)


      !
      ! Check the presence of transitions at each frequency
      !

      ! Reset the ranges that the transitions holds
      Frec%lif0 = 100000000
      Frec%lif1 = -1
      Frec%pif0 = 100000000
      Frec%pif1 = -1

      ! Initialize size of profile variable
      Frec%ntfreq = 0
      Frec%ntfreqi = 0
      Frec%npfreq = 0

      ! For each atom
      do ia=1,nA

        ! If H or He, allocate freqR and freqB
        if (Atom(ia)%Element.eq.' H'.or.Atom(ia)%Element.eq.'HE') then
          allocate(Atom(ia)%freqR(Atom(ia)%ntran))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%freqR)
          allocate(Atom(ia)%freqB(Atom(ia)%ntran))
          MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%freqB)
        end if

        ! For each b-b transitions
        do itran=1,Atom(ia)%ntran

          ! Apply atomic shift
          ktran = itran + Atom(ia)%tshift

          ! Get terms
          itermu = Atom(ia)%fst(itran)%itermu
          iterml = Atom(ia)%fst(itran)%iterml

          ! Initialize limits
          Atom(ia)%if0(itran) = 100000000
          Atom(ia)%if1(itran) = -1

          ! Check for each frequency if it is within limits
          do ifreq=1,nfreq

            ! Line present
            if (Frec%omega(ifreq).ge.tmp_lim(1,ktran).and. &
                Frec%omega(ifreq).le.tmp_lim(2,ktran)) then

              ! Add weight CRD transition
              Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + LCRDW

              ! If CLE
              if (run_mode.eq.2) &
                Frec%IW_freq_in(ifreq) = Frec%IW_freq_in(ifreq) + &
                                         LCRDW

              ! Update lower limit
              if (ifreq.lt.Atom(ia)%if0(itran)) &
                Atom(ia)%if0(itran) = ifreq

              ! Update upper limit
              if (ifreq.gt.Atom(ia)%if1(itran)) &
                Atom(ia)%if1(itran) = ifreq

            end if ! Line present

          end do ! All frequencies

          ! If only one frequency
          if (Atom(ia)%if1(itran).le.Atom(ia)%if0(itran)) then

            ! Remove transition
            Atom(ia)%fflag(itran)%absent = .True.
            Atom(ia)%if0(itran) = 1
            Atom(ia)%if1(itran) = 0
            Atom(ia)%W0(itran) = 0d0
            Atom(ia)%W1(itran) = 0d0

          ! Multiple frequencies
          else

            ! Line present
            Atom(ia)%fflag(itran)%absent = .False.

          end if

          ! Store real Master frequency (for slaves)
          ! The r ones  will be overwritten later, maybe
          ! The t ones are persistent for sure
          Atom(ia)%rif0(itran) = Atom(ia)%if0(itran)
          Atom(ia)%rif1(itran) = Atom(ia)%if1(itran)
          Atom(ia)%tif0(itran) = Atom(ia)%if0(itran)
          Atom(ia)%tif1(itran) = Atom(ia)%if1(itran)

          ! If absent, skip
          if (Atom(ia)%fflag(itran)%absent) cycle

          ! If H or He, save freqR and freqB
          if (Atom(ia)%Element.eq.' H'.or.Atom(ia)%Element.eq.'HE') then
            Atom(ia)%freqR = Frec%omega(Atom(ia)%if0(itran))
            Atom(ia)%freqB = Frec%omega(Atom(ia)%if1(itran))
          end if

          ! Compute weights for the limits
          Atom(ia)%W0(itran) = .5d5* &
                               (Frec%omega(Atom(ia)%if0(itran)+1) - &
                                Frec%omega(Atom(ia)%if0(itran)))
          Atom(ia)%W1(itran) = .5d5* &
                               (Frec%omega(Atom(ia)%if1(itran)) - &
                                Frec%omega(Atom(ia)%if1(itran)-1))

          ! Update limits for the range with any line
          if (Atom(ia)%if0(itran).lt.Frec%lif0) Frec%lif0 = &
                                                   Atom(ia)%if0(itran)
          if (Atom(ia)%if1(itran).gt.Frec%lif1) Frec%lif1 = &
                                                   Atom(ia)%if1(itran)

          ! Add frequencies to count for profile dimensions
          Frec%ntfreq = Frec%ntfreq + Atom(ia)%if1(itran) - &
                        Atom(ia)%if0(itran) + 1
          Frec%ntfreqi = Frec%ntfreqi + (Atom(ia)%if1(itran) - &
                         Atom(ia)%if0(itran) + 1)* &
                         Atom(ia)%fst(itran)%nt

        end do ! Bound-bound transitions

        ! For each b-f transitions
        do itran=1,Atom(ia)%nphot

          ! Get the maximum frequency
          nt = (1 - Atom(ia)%phot(itran)%mode)* &
               Atom(ia)%phot(itran)%nfreq + &
               Atom(ia)%phot(itran)%mode
          v = Atom(ia)%phot(itran)%infreq(nt)

          ! Initialize as present
          Atom(ia)%phot(itran)%absent = .False.

          ! Reset the ranges that the transition holds
          Atom(ia)%phot(itran)%if0 = -1
          Atom(ia)%phot(itran)%if1 = -1

          ! If explicit
          if (Atom(ia)%phot(itran)%mode.eq.0) then

            ! For each frequency
            do ifreq=1,nfreq

              ! If we are below range, keep searching
              if (Frec%omega(ifreq).lt. &
                  Atom(ia)%phot(itran)%infreq(1)) cycle

              ! If we are above maximum, identify the index and go out
              if (Frec%omega(ifreq).gt.v) then
                Atom(ia)%phot(itran)%if1 = ifreq - 1
                exit
              end if

              ! If we are above the minimum and not initialized
              if (Atom(ia)%phot(itran)%if0.lt.0) then

                ! Update index
                if (Frec%omega(ifreq).ge. &
                    Atom(ia)%phot(itran)%infreq(1)) &
                  Atom(ia)%phot(itran)%if0 = ifreq

              end if

            end do ! Frequencies

          ! If hydrogenic
          else

            ! For each frequency
            do ifreq=1,nfreq

              ! If we are below the edge, keep searching
              if (Frec%omega(ifreq).lt.Atom(ia)%phot(itran)%edge) &
                cycle

              ! If we are above maximum, identify the index and go out
              if (Frec%omega(ifreq).gt.v) then
                Atom(ia)%phot(itran)%if1 = ifreq - 1
                exit
              end if

              ! If we are above the minimum and initial index not set
              if (Atom(ia)%phot(itran)%if0.lt.0) then

                ! Get initial index
                if (Frec%omega(ifreq).ge.Atom(ia)%phot(itran)%edge) &
                  Atom(ia)%phot(itran)%if0 = ifreq

              end if

            end do ! Frequencies

          end if ! Type of b-f transition

          ! If you didn't find the maximum, that means it was the
          ! last frequency itself
          if (Atom(ia)%phot(itran)%if1.lt.0) &
            Atom(ia)%phot(itran)%if1 = nfreq

          ! Limits for the range with photoionizations
          if (Atom(ia)%phot(itran)%if0.lt.Frec%pif0) Frec%pif0 = &
                                              Atom(ia)%phot(itran)%if0
          if (Atom(ia)%phot(itran)%if1.gt.Frec%pif1) Frec%pif1 = &
                                              Atom(ia)%phot(itran)%if1

          ! Add frequencies to count for profile dimensions
          Frec%npfreq = Frec%npfreq + &
                        (Atom(ia)%phot(itran)%if1 - &
                         Atom(ia)%phot(itran)%if0 + 1)

          ! For each frequency
          do ifreq=Atom(ia)%phot(itran)%if0,Atom(ia)%phot(itran)%if1

            ! Add photoionization weight to the node
            Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + PhotW

          end do

          ! If CLE
          if (run_mode.eq.2) then

            ! For each frequency
            do ifreq=Atom(ia)%phot(itran)%if0,Atom(ia)%phot(itran)%if1

              ! Add weight to the node of the input axis
              Frec%IW_freq_in(ifreq) = Frec%IW_freq_in(ifreq) + PhotW

            end do ! Frequencies

          end if ! CLE

        end do ! Bound-free transitions
      end do ! Atoms

      ! For each LTE line
      do ia=1,nLTEl

        ! Initialize limits
        Input%LTEline(ia)%if0 = nfreq+1
        Input%LTEline(ia)%if1 = -1

        ! Check for each frequency if it is within limits
        do ifreq=1,nfreq

          ! Line present
          if (Frec%omega(ifreq).ge.tmp_limL(1,ia).and. &
              Frec%omega(ifreq).le.tmp_limL(2,ia)) then

            ! Add weight of LTE line
            Frec%IW_freq(ifreq) = Frec%IW_freq(ifreq) + LLTEW

            ! Update lower limit
            if (ifreq.lt.Input%LTEline(ia)%if0) &
              Input%LTEline(ia)%if0 = ifreq

            ! Update upper limit
            if (ifreq.gt.Input%LTEline(ia)%if1) &
              Input%LTEline(ia)%if1 = ifreq

          end if ! Line present

        end do ! Frequencies

        ! If only one frequency
        if (Input%LTEline(ia)%if1.le.Input%LTEline(ia)%if0) then

          ! Remove the transition
          Input%LTEline(ia)%absent = .True.
          Input%LTEline(ia)%if0 = 1
          Input%LTEline(ia)%if1 = 0

        ! More frequencies
        else

          ! Line is present
          Input%LTEline(ia)%absent = .False.

        end if ! Number of frequencies

      end do ! LTE lines


      !
      ! Protect frequencies against jumps due to big range
      !

      ! Allocate and initialize protection
      allocate(protect(Nfreq))
      protect = .False.

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! If between terms, we cannot protect
          if (Atom(ia)%fst(itran)%nt.gt.1) cycle

          ! For every frequency in this transition
          do ifreq=Atom(ia)%if0(itran)+1,Atom(ia)%if1(itran)

            ! Protect it
            protect(ifreq) = .True.

          end do ! Each frequency
        end do ! b-b transitions

        ! For each b-f transitions
        do itran=1,Atom(ia)%nphot

          ! For each frequency
          do ifreq=Atom(ia)%phot(itran)%if0+1, &
                   Atom(ia)%phot(itran)%if1

            ! Protect ionizations
            protect(ifreq) = .True.

          end do ! Each frequency
        end do ! b-b transitions
      end do ! Atom


      !
      ! Define the integration weights
      !

      ! The first point is special in compound trapezoidal rule
      Frec%W_freq(1) = .5d0*(Frec%omega(2) - Frec%omega(1))

      ! Initialize the integral to normalize the weights
      norm1 = Frec%W_freq(1)

      ! The initial lower limit is the first point
      O0 = Frec%omega(1)

      ! This is the pointer to the first element of the current
      ! interval, we are pointing to the first element
      cfreq = 1

      ! Flag that says that the point 2 is not the initial point
      ! of the interval (because 1 is the initial point)
      init = .False.

      ! For the rest of frequencies except the last
      do ifreq=2,Nfreq-1

        ! If ifreq is the initial point of an interval
        if (init) then

          ! The first point is special in compound trapezoidal rule
          Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq+1) - &
                                       Frec%omega(ifreq))

          ! The next point cannot be a first point
          init = .False.

          ! Initialize the integral to normalize the weights
          norm1 = Frec%W_freq(ifreq)

          ! Pointer is now in this frequency
          cfreq = ifreq

          ! And it is the beginning of the current interval
          O0 = Frec%omega(ifreq)

        ! If ifreq is not the initial point of an interval
        else

          ! If there is a big jump (parameter) and the frequency is
          ! not protected
          if (abs(1d2/Frec%omega(ifreq+1) - &
                  1d2/Frec%omega(ifreq)).gt.jump.and. &
              .not.protect(ifreq+1)) then

            ! Master
            if (gpid.eq.0) then

              ! Notify of this jump
              write(umsg,'(A,1x,f10.3,1x,A,1x,f10.3,1x,A)') &
                  ' # Wavelengths ',1d2/Frec%omega(ifreq+1),'and', &
                  1d2/Frec%omega(ifreq),'are considered to be'// &
                  ' in different ranges.'
              call verbose

            end if

            ! The last point is special in compound trapezoidal rule
            Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq) - &
                                         Frec%omega(ifreq-1))

            ! It is the end of the current interval
            O1 = Frec%omega(ifreq)

            ! Add to the integral
            norm1 = norm1 + Frec%W_freq(ifreq)

            ! We know that the integral must be
            ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
            norm = 1d5*(O1 - O0)

            ! The normalizing factor is thus
            norm = norm/norm1

            ! For each frequency in this interval
            do jfreq=cfreq,ifreq

              ! Normalize the weights of this interval
              Frec%W_freq(jfreq) = Frec%W_freq(jfreq)*norm

            end do ! Frequencies in interval

            ! The next point is the first point of its interval
            init = .True.

          ! If ifreq is not the last point of an interval
          else

            ! Check if ifreq could be the last point of an interval
            ! but it is protected
            if (abs(1d2/Frec%omega(ifreq+1) - &
                    1d2/Frec%omega(ifreq)).gt.jump.and. &
                protect(ifreq+1).and.gpid.eq.0) then

              ! Issue message
              write(umsg,'(A,1x,f10.3,1x,A,1x,f10.3,1x,A)') &
                  ' # Wavelengths ',1d2/Frec%omega(ifreq+1),'and', &
                  1d2/Frec%omega(ifreq),'are considered to be'// &
                  ' in different ranges, but they are '// &
                  'protected from the splitting.'
              call verbose

            end if

            ! Compound trapezoidal rule weight
            Frec%W_freq(ifreq) = .5d0*(Frec%omega(ifreq+1) - &
                                       Frec%omega(ifreq-1))

            ! Add to the integral
            norm1 = norm1 + Frec%W_freq(ifreq)

          end if ! Last point
        end if ! Initial point

      end do ! Frequencies

      ! The last point is special in compound trapezoidal rule
      Frec%W_freq(nfreq) = .5d0*(Frec%omega(nfreq) - &
                                 Frec%omega(nfreq-1))

      ! It is the end of the interval
      O1 = Frec%omega(Nfreq)

      ! Add to the integral
      norm1 = norm1 + Frec%W_freq(Nfreq)

      ! We know that the integral must be
      ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
      norm = 1d5*(O1 - O0)

      ! The normalizing factor is thus
      norm = norm/norm1

      ! For all frequencies in this interval
      do ifreq=cfreq,nfreq

        ! Normalize the weights of this interval
        Frec%W_freq(ifreq) = Frec%W_freq(ifreq)*norm

      end do ! Frequencies in the interval

      ! Deallocate arrays
      if (allocated(tmp_lim)) deallocate(tmp_lim)
      deallocate(protect)

      ! Check if everything is fine
      call control

      return

1000  write(umsg,'(A,1x,i2)') 'Error opening wave file',i
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  write(umsg,'(A,1x,i2)') 'Error reading wave file',i
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine omegabuild

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize frequency related variables needed by the
      !! radiative transfer master\n
      !!  Atom(Atom_class(:)): Structures with atomic data
      subroutine omegainitmaster(Atom)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom

      ! Local

      integer:: ia,itran


      ! Only the Máster with MPI
      if (pid.eq.0.and.nproc.gt.1) then

        ! For each atom
        do ia=1,nA

          ! For each b-b transition
          do itran=1,Atom(ia)%ntran

            ! Allocate and initialize the presence flag for master
            allocate(Atom(ia)%fflag(itran)%Mabsent(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%fflag(itran)%Mabsent)
            Atom(ia)%fflag(itran)%Mabsent = .False.

            ! Initialize
            Atom(ia)%Mif0(itran,:) = Atom(ia)%if0(itran)
            Atom(ia)%Mif1(itran,:) = Atom(ia)%if1(itran)
            Atom(ia)%MW0(itran,:) = Atom(ia)%W0(itran)
            Atom(ia)%MW1(itran,:) = Atom(ia)%W1(itran)

          end do ! b-b transitions

          ! For each b-f tarnsition
          do itran=1,Atom(ia)%nphot

            ! Allocate and initialize the presence flag for master
            allocate(Atom(ia)%phot(itran)%Mabsent(0:nproc-1))
            MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%phot(itran)%Mabsent)
            Atom(ia)%phot(itran)%Mabsent = .False.

            ! Initialize master
            Atom(ia)%phot(itran)%Mif0 = Atom(ia)%phot(itran)%if0
            Atom(ia)%phot(itran)%Mif1 = Atom(ia)%phot(itran)%if1

          end do ! b-f transitions
        end do ! Atoms

      end if ! Máster with MPI

      end subroutine omegainitmaster

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine the frequency limits for an output transition in
      !! the second order emissivity\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Red(Red_class): Structure with redistribution input
      !!                    frequency data, redistribution function
      !!                    data, and profile or normalization data\n
      !!       ia(integer): Index of current atom\n
      !!    jtran(integer): Index of current transition\n
      !!       nt(integer): Number of components to check\n
      !!      if0(integer): Initial frequency index for this
      !!                    transition\n
      !!      if1(integer): Final frequency index for this
      !!                    transition\n
      !!  red_negl(double): Doppler width distance to neglect
      !!                    frequencies in the second order
      !!                    emissivity\n
      !!               DwN: Nominal Doppler width\n
      !!    nut(double(:)): Array of output frequencies\n
      !!  omega(double(:)): Frequency axis
      subroutine get_transition_out_limit(Atmo,Red,ia,jtran,nt, &
                                          if0,if1,red_negl, &
                                          DwN,nut,omega)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Red_class), intent(inout):: Red
      integer, intent(in):: ia,jtran,nt,if0,if1
      double precision, intent(in):: red_negl,DwN
      double precision, dimension(:), intent(in):: nut,omega

      ! Local

      logical:: lskip,skip

      integer:: iz,indx,it,ifreq,nran,bf0,bf1,np

      double precision:: vfacm,vfacp,Dw,red_neglW


      ! Get index
      indx = Red%izao(jtran,ia,Rz0)

      ! Doppler shift factor
      vfacm = 1d0
      vfacp = 1d0

      ! If dynamic
      if (dyn) then

        ! Maximum velocity
        Dw = maxval(sqrt(Atmo%vx*Atmo%vx + &
                         Atmo%vy*Atmo%vy + &
                         Atmo%vz*Atmo%vz))

        ! If big enough
        if (Dw.gt.TINYVEL) then

          ! Calculate Doppler shift factors
          vfacp = 1d0 + Dw
          vfacm = 1d0 - Dw

        end if
      end if ! Dynamic

      ! Transform the searching parameters from normalized
      ! to proper frequency units
      red_neglW = red_negl*DwN

      ! Initialize search parameters
      nran = 0
      bf0 = -1
      bf1 = -2

      ! Only one output
      if (if0.eq.if1) then

        ! For each line component
        do it=1,nt

          ! Check if we are close to a transition frequency
          if (abs(omega(if0)*vfacm-nut(it)).lt.red_neglW.or. &
              abs(omega(if0)*vfacp-nut(it)).lt.red_neglW) then

            ! Only trivial range
            bf0 = if0
            bf1 = if1
            nran = 1
            exit

          end if ! Close to a transition
        end do ! Transitions

      ! More than one output
      else

        ! Reset logical variable
        lskip = .True.

        ! Look for the limits, for each output frequency
        do ifreq=if0,if1

          ! Initialize the flag
          skip = .True.

          ! For each transition component
          do it=1,nt

            ! Check if we are close to a transition frequency
            if (abs(omega(ifreq)*vfacm-nut(it)).lt.red_neglW.or. &
                abs(omega(ifreq)*vfacp-nut(it)).lt.red_neglW) then

              ! We cannot skip this one
              skip = .False.
              exit

            end if ! Close to a transition component

          end do ! Transition components

          ! If we skip this frequency or it is the last
          if (skip.or.ifreq.eq.if1) then

            ! If we have found anything before
            if (.not.lskip) then

              ! If it is the last
              if (ifreq.eq.if1) then

                ! End of range is this
                bf1 = ifreq

              ! Not the last
              else

                ! Previous one was end of range
                bf1 = ifreq - 1

              end if ! If last index

              ! Add a new range
              nran = nran + 1

            end if ! We have found anything before

          ! If we cannot skip this frequency
          else

            ! If this is the first index for this
            ! transition
            if (bf0.lt.0) bf0 = ifreq

          end if ! Can skip

          ! Update status of last frequency
          lskip = skip

        end do ! output frequencies

      end if ! Number of outputs

      ! Store number of ranges in structure
      Red%ao(indx)%nran = nran

      ! If no ranges
      if (Red%ao(indx)%nran.lt.1) then

        ! No frequencies
        Red%ao(indx)%nran = 0
        Red%ao(indx)%nfreq = 0
        Red%ao(indx)%Mi0 = 0
        Red%ao(indx)%Mi1 = -1
        Red%ao(indx)%nn = 0
        return

      end if

      ! Allocate the range indexes
      allocate(Red%ao(indx)%if0(nran))
      Red%ao(indx)%if0 = bf0
      allocate(Red%ao(indx)%if1(nran))
      Red%ao(indx)%if1 = bf1

      ! Only one frequency
      if (if0.eq.if1) then

        ! Single range
        nran = 1
        np = 1
        Red%ao(indx)%if0(nran) = if0
        Red%ao(indx)%if1(nran) = if0

      ! More than one frequency
      else

        ! Only one output
        if (bf1.eq.bf0) then

          ! Initialize the flag
          skip = .True.

          ! For each transition component
          do it=1,nt

            ! Check if we are close to a transition frequency
            if (abs(omega(bf0)*vfacm-nut(it)).lt.red_neglW.or. &
                abs(omega(bf0)*vfacp-nut(it)).lt.red_neglW) then

              ! Cannot skip this one
              skip = .False.
              exit

            end if ! Close to a transition component

          end do ! Transition components

          ! No transition
          if (skip) then

            ! Zero out
            Red%ao(indx)%nran = 0
            deallocate(Red%ao(indx)%if0)
            deallocate(Red%ao(indx)%if1)
            Red%ao(indx)%Mi0 = 0
            Red%ao(indx)%Mi1 = -1
            Red%ao(indx)%nn = 0
            return

          end if ! No transition

          ! Save limits
          Red%ao(indx)%if0(1) = bf0
          Red%ao(indx)%if1(1) = bf0
          np = 1
          nran = 1

        ! Several outputs
        else

          ! Reset logical variable
          lskip = .True.

          ! Initialize counters
          nran = 0
          np = 0

          ! For each output frequency in the pre-checked limits
          do ifreq=bf0,bf1

            ! Initialize the flag
            skip = .True.

            ! For each transition component
            do it=1,nt

              ! Check if we are close to a transition frequency
              if (abs(omega(ifreq)*vfacm-nut(it)).lt.red_neglW.or. &
                  abs(omega(ifreq)*vfacp-nut(it)).lt.red_neglW) then

                ! Cannot skip this frequency
                skip = .False.
                exit

              end if ! Close to a transition component

            end do ! Transition components

            ! If we skip this frequency or it is the last
            if (skip.or.ifreq.eq.bf1) then

              ! Did not skip the last
              if (.not.lskip) then

                ! It is the last
                if (ifreq.eq.bf1) then

                  ! Last is this one then
                  Red%ao(indx)%if1(nran) = ifreq
                  np = np + 1

                ! It is not the last
                else

                  ! Previous one closes
                  Red%ao(indx)%if1(nran) = ifreq - 1

                end if ! Last index
              end if ! Did skip the last

            ! If we cannot skip this frequency
            else

              ! Advance index
              np = np + 1

              ! If skipped last
              if (lskip) then

                ! Initialize next range
                nran = nran + 1
                Red%ao(indx)%if0(nran) = ifreq

              end if ! Skipped last
            end if ! Skip this frequency

            ! Update skipped status
            lskip = skip

          end do ! output frequencies

        end if ! Number of outputs
      end if ! Number of outputs

      ! Set global limits
      Red%ao(indx)%gf0 = minval(Red%ao(indx)%if0)
      Red%ao(indx)%gf1 = maxval(Red%ao(indx)%if1)
      Red%ao(indx)%tgf0 = Red%ao(indx)%gf0
      Red%ao(indx)%tgf1 = Red%ao(indx)%gf1

      ! Count frequencies
      Red%ao(indx)%nfreq = np

      ! Initialize for heights
      do iz=Rz0,Rz1_PRD

        ! Get index
        indx = Red%izao(jtran,ia,iz)

        ! Initialize limits
        Red%zao(indx)%ggf0 = 10000000
        Red%zao(indx)%ggf1 = -1

      end do ! Heights

      end subroutine get_transition_out_limit

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine the input frequency axis given the transition
      !! resonances and the discretization parameters of a pair of
      !! output and input transitions\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Redd_aux(Redc_class): Auxiliar to store the input frequency
      !!                        axis for the calculation of the
      !!                        second order emissivity\n
      !!        icohw(logical): If coherent scattering in the wings\n
      !!           ia(integer): Index for current atom\n
      !!        jtran(integer): Index for current output transition\n
      !!          iti(integer): Rolling index for current input
      !!                       transition\n
      !!           nr(integer): Number of Raman resonances\n
      !!          ntj(integer): Number of output transition
      !!                        components\n
      !!           nt(integer): Number of output+input transition
      !!                        components\n
      !!        Dfreqi(double): Input transition frequency\n
      !!        Dfreqo(double): Output transition frequency\n
      !!         dcohw(double): Doppler width distance to start
      !!                        considering coherent scattering\n
      !!   red_pars(double(:)): Discretization parameters to
      !!                        determine the input frequency axis\n
      !!      omega(double(:)): Frequency axis\n
      !!        dnl(double(:)): Raman resonances\n
      !!        nut(double(:)): Transition components frequencies
      subroutine get_input_frequencies(Atmo,Red,Redd_aux,icohw,ia, &
                                       jtran,iti,nr,ntj,nt,np0, &
                                       Dfreqi,Dfreqo,dcohw,red_pars, &
                                       omega,dnl,nut,cDopp,vphv, &
                                       vplv,vphve,vplve,vpr,vpp,flag)
      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Red_class), intent(in):: Red
      type(Redc_class), dimension(:,:), &
                        allocatable, intent(inout):: Redd_aux
      logical, intent(in):: icohw
      integer, dimension(:), allocatable, intent(inout):: flag
      integer, intent(inout):: np0
      integer, intent(in):: ia,jtran,iti,nr,ntj,nt
      double precision, intent(in):: Dfreqi,Dfreqo,dcohw,cDopp
      double precision, dimension(:), intent(in):: omega,dnl,nut
      double precision, dimension(:), intent(out):: vphv,vplv
      double precision, dimension(:), intent(out):: vphve,vplve,vpr
      double precision, dimension(:), intent(in):: red_pars
      double precision, dimension(:), allocatable, intent(inout):: vpp

      ! Local

      logical:: cohw,core,reset,nfound,init

      integer:: iz,indx,jndx,iifreq,iran,ifreq,it,ir,jfreq,kfreq,cfreq
      integer:: np,ni,nie,ip,ipp,npp,nn

      double precision:: vfacm,vfacp,DwT,Dw1,Dw
      double precision:: dnlmin,dnlmax,vph,vpl,norm,O1,norm1,O0
      double precision:: red_cohwW,red_resoW,red_neglW,red_coreW
      double precision:: red_rangwW1,red_vlarwW1,red_fstpwW1
      double precision:: red_mstpwW1,red_rangcW1,red_vlarcW1
      double precision:: red_fstpcW1,red_mstpcW1
      double precision:: red_rangW1,red_vlarW1
      double precision:: red_fstpW1,red_mstpW1
      double precision, dimension(:), allocatable:: Wvpp

      ! Pointers

      type(dbabox_class), pointer:: bomega, bw_freq, bdaux


      ! Get limits for Raman resonance
      dnlmin = minval(dnl(1:nr))
      dnlmax = maxval(dnl(1:nr))

      ! Get index
      jndx = Red%izao(jtran,ia,Rz0)

      ! If the CPU does not have ranges, leave
      if (Red%ao(jndx)%nran.lt.1) return

      ! For each height node
      do iz=Rz0,Rz1_PRD

        ! Get index
        indx = Red%izao(jtran,ia,iz)

        ! Init
        vfacm = 1d0
        vfacp = 1d0

        ! If dynamic
        if (dyn) then

          ! Compute velocity
          Dw = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                    Atmo%vy(iz)*Atmo%vy(iz) + &
                    Atmo%vz(iz)*Atmo%vz(iz))

          ! Big enough
          if (Dw.gt.TINYVEL) then

            ! Get Doppler shift
            vfacp = 1d0 + Dw
            vfacm = 1d0 - Dw

          end if ! Enough velocity
        end if ! Dynamic

        !
        ! Calculate the Doppler width of both input and output
        ! transitions
        !

        ! Thermal common part
        DwT = cDopp*sqrt(Atmo%T(iz))

        ! Input
        Dw1 = Dfreqi*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        ! Output
        Dw  = Dfreqo*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        ! Transform the searching parameters from normalized to proper
        ! frequency units
        red_cohwW = dcohw*Dw
        red_resoW = red_pars(2)*Dw
        red_neglW = red_pars(3)*Dw
        red_coreW = red_pars(7)*Dw
        red_rangwW1 = red_pars(1)*Dw1
        red_vlarwW1 = red_pars(4)*Dw1
        red_fstpwW1 = red_pars(5)*Dw1
        red_mstpwW1 = red_fstpwW1*red_pars(6)
        red_rangcW1 = red_pars(8)*Dw1
        red_vlarcW1 = red_pars(9)*Dw1
        red_fstpcW1 = red_pars(10)*Dw1
        red_mstpcW1 = red_fstpcW1*red_pars(11)

        ! Allocate input frequency size
        np = Red%ao(jndx)%nfreq
        allocate(Redd_aux(iti,iz)%mfreq(np))

        ! Allocate and initialize pointers and back dimension trace
        allocate(bomega, bw_freq)
        nullify(bomega%next)
        nullify(bomega%prev)
        nullify(bw_freq%next)
        nullify(bw_freq%prev)
        bomega%nback = 0
        bw_freq%nback = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Red%ao(jndx)%nran
          do ifreq=Red%ao(jndx)%if0(iran), &
                   Red%ao(jndx)%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

      !
      ! Reset indentation
      !

      ! Reset number of frequencies
      Redd_aux(iti,iz)%mfreq(iifreq) = 0

      ! Coherent wings?
      if (icohw) then

        ! Initialize flag
        cohw = .True.

        ! For all output transition components
        do it=1,ntj

          ! Check if close to the transition component
          if (abs(nut(it) - omega(ifreq)*vfacm).lt.red_cohwW.or. &
              abs(nut(it) - omega(ifreq)*vfacp).lt.red_cohwW) then

            ! No possible to apply coherent scattering
            cohw=.False.
            exit

          end if ! Close to a transition component

        end do ! Output transition components

      ! Non-coherent wings
      else

        ! Flag non-coherent wings
        cohw = .False.

      end if ! Coherent wings?

      ! If Coherent wing
      if (cohw) then

        ! Definitely not a core frequency (it has to be wing!)
        core = .False.

        !
        ! Store the determined vector
        !

        ! Current boxes are empty
        if (.not.allocated(bomega%A)) then

          ! Allocate single frequency
          allocate(bomega%A(1))
          allocate(bw_freq%A(1))

        ! Need new boxes
        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(1))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(1))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if ! Empty boxes

        ! Store the fequency axis
        bomega%A = omega(ifreq) + dnlmin

        ! Store the dimension of the axis
        bomega%mfreq = 0
        Redd_aux(iti,iz)%mfreq(iifreq) = 0

        ! Store indexes
        bomega%ifreq = iifreq

      ! Non-coherent wing
      else

        ! Reset core flag
        core = .False.

        ! For every output transition component
        do it=1,ntj

          ! Flag as core if close enough to any resonance
          if (abs(nut(it) - omega(ifreq)*vfacm).le.red_coreW.or. &
              abs(nut(it) - omega(ifreq)*vfacp).le.red_coreW) &
            core=.True.

          ! Initialize limits
          vphv(it) = -1D0
          vplv(it) = -1D0

        end do ! Transition components

        ! If it is core
        if (core) then

          ! Get rest of parameters from core inputs
          red_rangW1 = red_rangcW1
          red_vlarW1 = red_vlarcW1
          red_fstpW1 = red_fstpcW1
          red_mstpW1 = red_mstpcW1

        ! If it is wing
        else

          ! Get rest of parameters from wing inputs
          red_rangW1 = red_rangwW1
          red_vlarW1 = red_vlarwW1
          red_fstpW1 = red_fstpwW1
          red_mstpW1 = red_mstpwW1

        end if

        ! For every input transition component
        do it=ntj+1,nt

          ! Initialize limits
          vphv(it) = nut(it) + red_rangW1
          vplv(it) = nut(it) - red_rangW1

        end do

        ! For every resonance
        do ir=1,nr

          ! Find real resonance
          vpr(ir) = omega(ifreq) + dnl(ir)

          ! Initialize range limits
          vphv(nt + ir) = vpr(ir) + red_rangW1
          vplv(nt + ir) = vpr(ir) - red_rangW1

        end do ! Resonances


        !
        ! Define the true limits:
        ! The resonances specify the true limits, but if
        ! they are close to a transition, we expand the limit
        !

        ! Take the furthest resonances
        vph = omega(ifreq) + dnlmax
        vpl = omega(ifreq) + dnlmin

        ! For each input transition component
        do it=ntj+1,nt

          ! If we are close to the transition
          if (abs(vph - nut(it)).lt.red_resoW.or. &
              abs(vpl - nut(it)).lt.red_resoW) then

            ! Expand suitable limit to include the resonance
            if (nut(it).lt.vpl) vpl = nut(it)
            if (nut(it).gt.vph) vph = nut(it)

          end if ! Close to the transition

        end do ! Output components

        ! Now expand the range from the parameters
        vph = vph + red_rangW1
        vpl = vpl - red_rangW1


        !
        ! Flag lines and resonances out of limits
        !

        ! The total number of resonances that we had
        ni = nr + nt

        ! We are going to change ni, so store it because
        ! we need the original value later
        np = ni

        ! For each resonance and transition
        do ir=1,np

          ! Check if it is out of range
          if (vphv(ir).lt.vpl.or.vplv(ir).gt.vph) then

            ! Remove the resonance
            vplv(ir) = vpl - 1
            ni = ni - 1

          ! It is not out of range, but the lower limit is
          ! out
          else if(vplv(ir).lt.vpl)then

            ! Adjust lower limit
            vplv(ir) = vpl

          ! It is not out of range, but the upper limit is
          ! out
          else if(vphv(ir).gt.vph)then

            ! Adjust upper limit
            vphv(ir) = vph

          end if ! Compare ranges

        end do ! Resonances and transitions


        !
        ! Check lines and resonances holding the same range
        !

        ! For each transition and resonance
        do ir=np,2,-1

          ! Skip if already flagged out
          if (vplv(ir).lt.vpl) cycle

          ! Any other transition and resonance not already checked as
          ! a pair
          do it=ir-1,1,-1

            ! Skip if already flagged out
            if (vplv(it).lt.vpl) cycle

            ! If the ranges overlap
            if ((vphv(ir).ge.vplv(it).and.vphv(ir).le.vphv(it)).or. &
                (vplv(ir).ge.vplv(it).and.vplv(ir).le.vphv(it))) then

              ! Combine ranges into just one
              vplv(it) = min(vplv(ir),vplv(it))
              vphv(it) = max(vphv(ir),vphv(it))

              ! Flag the second out
              vplv(ir) = vpl - 1d0
              ni = ni - 1
              exit

            end if ! Ranges overlap

          end do ! Pair transition/resonance
        end do ! Transition components and resonances


        !
        ! Shift the individual limits moving the valid ones
        ! to the first part of the vector
        !

        ! If we have changed the number of ranges from the
        ! beginning
        if (ni.ne.np) then

          ! For each pair of resonances
          do ir=1,np-1

            ! If it is flagged out
            if (vplv(ir).lt.vpl) then

              ! Reset the added index
              ip = 1

              ! For all the resonances in front of this one
              do it=ir+1,np

                ! If it is flagged, skip it
                if (vplv(it).lt.vpl) then

                  ! Advance and do nothing
                  ip = ip + 1

                ! If it is not flagged
                else

                  ! Move it to the position of the flagged one
                  vplv(it-ip) = vplv(it)
                  vphv(it-ip) = vphv(it)

                  ! And flag this position
                  vplv(it) = vpl - 1
                  exit

                end if ! it flagged

              end do ! All resonances and transitions in front

            end if ! ir flagged

          end do ! All resonances and transitions

        end if ! if something flagged


        !
        ! Order the individual limits
        !

        ! Lower limits
        call QsortC(vplv(1:ni))
        ! Upper limits
        call QsortC(vphv(1:ni))


        !
        ! Define extended limits
        !

        ! Reset the index of extended limits
        nie = 0

        ! For each normal limit
        do ir=1,ni

          ! Get an extended range to the left
          nie = nie + 1
          vplve(nie) = vplv(ir) - red_vlarW1
          vphve(nie) = vplv(ir) - red_fstpW1

          ! Get an extended range to the right
          nie = nie + 1
          vplve(nie) = vphv(ir) + red_fstpW1
          vphve(nie) = vphv(ir) + red_vlarW1

        end do


        !
        ! Check lines and resonances holding the same range
        !

        ! Store the number of original limits
        np = nie

        ! For each transition and resonance
        do ir=np,2,-1

          ! Skip if already flagged out
          if (vplve(ir).lt.0d0) cycle

          ! Each other transition and resonance not already paired
          do it=ir-1,1,-1

            ! Skip if already flagged out
            if (vplve(it).lt.0d0) cycle

            ! If the ranges overlap
            if ((vphve(ir).ge.vplve(it).and. &
                 vphve(ir).le.vphve(it)).or. &
                (vplve(ir).ge.vplve(it).and. &
                 vplve(ir).le.vphve(it))) then

              ! Combine them into just one
              vplve(it) = min(vplve(ir),vplve(it))
              vphve(it) = max(vphve(ir),vphve(it))

              ! Flag the second limit out
              vplve(ir) = -1d0
              nie = nie - 1
              exit

            end if ! Overlap

          end do ! Any other transition and resonance
        end do ! Transition components and resonances


        !
        ! Shift the individual limits moving the valid ones
        ! to the first part of the vector
        !

        ! If we have changed the number of ranges from the
        ! beginning
        if (nie.ne.np) then

          ! For each transition or resonance
          do ir=1,np-1

            ! If it is flagged out
            if (vplve(ir).lt.0d0) then

              ! Reset the added index
              ip = 1

              ! For all the resonances in front of this one
              do it=ir+1,np

                ! If it is flagged
                if (vplve(it).lt.0d0) then

                  ! Advance position
                  ip = ip + 1

                ! If it is not flagged
                else

                  ! Move it to the position of the flagged one
                  vplve(it-ip) = vplve(it)
                  vphve(it-ip) = vphve(it)

                  ! And flag this one
                  vplve(it) = -1d0
                  exit

                end if ! Flagged

              end do ! Resonances in front

            end if ! Flagged

          end do ! Resonances

        end if ! if something flagged


        !
        ! Order the individual limits
        !

        ! Lower limits
        call QsortC(vplve(1:nie))
        ! Upper limits
        call QsortC(vphve(1:nie))


        !
        ! Build the vector of input frequencies from the
        ! limits
        !

        ! Flag to reset in case we run out of space while
        ! building a vector
        reset = .False.

        ! Do until we are finished
        do while (.True.)

          ! If flagged, we need more space for the vector
          if (reset) then

            ! Duplicate needed space
            np0 = np0*2

            ! Reallocate
            deallocate(vpp)
            allocate(vpp(np0))
            deallocate(flag)
            allocate(flag(np0))

            ! Do not reset again unless running out of space again
            reset = .False.

          end if ! Ran out of space and had to reset

          ! Reset index counter
          np = 0

          !
          ! Build for the short limits
          !

          ! For each short range
          do it=1,ni

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if (np.gt.np0) then
              reset = .True.
              exit
            end if

            ! Start with the lower limit
            vpp(np) = vplv(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_fstpW1

              ! If we are over the range, we are done
              if (vpp(np).ge.vphv(it)) then
                vpp(np) = vphv(it)
                exit
              end if

            end do ! Until finished

            ! If we have no space, we need to reset
            if (reset) exit

          end do ! Short ranges

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Build for the extended limits
          !

          ! For each long range
          do it=1,nie

            ! Advance the index
            np = np + 1

            ! If we ran out of space, we have to reset
            if(np.gt.np0)then
              reset = .True.
              exit
            end if

            ! We start with the lower limit
            vpp(np) = vplve(it)

            ! Do until finished
            do while(.True.)

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Next frequency is the previous plus the step
              vpp(np) = vpp(np-1) + red_mstpW1

              ! If we are over the range, we are done
              if (vpp(np).ge.vphve(it)) then
                vpp(np) = vphve(it)
                exit
              end if

            end do ! Until done

            ! If we have no space, we need to allocate it above
            if (reset) exit

          end do ! Extended ranges

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          !
          ! Add the resonance frequencies to the vector
          !

          ! The number of frequencies we already have
          npp = np

          ! For each resonance
          do ir=1,nr

            ! Reset the flag
            nfound = .True.

            ! Run over the existing frequencies
            do ip=1,npp

              ! If the frequency is there
              if (abs(1d2/vpp(ip) - 1d2/vpr(ir)).lt.resolin) then

                ! Flag found because we are not adding it
                nfound = .False.
                exit

              end if ! Resonance already in list

            end do ! Existing frequencies

            ! If we did not find it
            if (nfound) then

              ! Advance the index
              np = np + 1

              ! If we ran out of space, we have to reset
              if(np.gt.np0)then
                reset = .True.
                exit
              end if

              ! Add the frequency
              vpp(np) = vpr(ir)

            end if ! Resonance not in axis

            ! If we have no space, we need to allocate it
            ! above
            if (reset) exit

          end do ! Resonances

          ! If we have no space, we need to allocate it above
          if (reset) cycle

          ! If we get to this point, we have everything
          exit

        end do ! Dummy loop to allow for resets


        !
        ! Check for duplicates
        !

        ! Reset the flag
        flag(1:np) = 1

        ! For each frequency
        do ip=1,np

          ! If it has been flagged, we already checked
          if (flag(ip).lt.1) cycle

          ! Check the other frequencies
          do ipp=ip+1,np

            ! If it has been flagged, we already checked
            if (flag(ipp).lt.1) cycle

            ! If same frequency, flag second to remove
            if (abs(1d2/vpp(ip)-1d2/vpp(ipp)).lt.resolin) &
              flag(ipp) = 0

          end do ! Other frequencies
        end do ! Frequencies

        ! Reset the running real index
        ipp = 0

        ! For each frequency in the vector
        do ip=1,np

          ! If it is flagged correct
          if (flag(ip).gt..5) then

            ! Advance index and add to real vector
            ipp = ipp + 1
            vpp(ipp) = vpp(ip)

          end if ! Flagged correct

        end do ! All frequencies

        ! The number of frequencies is the last value of ipp
        np = ipp

        ! Order the frequency axis
        call QsortC(vpp(1:np))

        !
        ! Store the determined vector
        !

        ! Empty box
        if (.not.allocated(bomega%A)) then

          ! Allocate space in current box
          allocate(bomega%A(np))
          allocate(bw_freq%A(np))

        ! Filled box, advance the boxes
        else

          ! bomega
          bdaux => bomega
          allocate(bomega%next)
          bomega => bomega%next
          bomega%prev => bdaux
          allocate(bomega%A(np))
          nullify(bomega%next)
          nullify(bdaux)
          bomega%nback = bomega%prev%nback + bomega%prev%mfreq

          ! bw_freq
          bdaux => bw_freq
          allocate(bw_freq%next)
          bw_freq => bw_freq%next
          bw_Freq%prev => bdaux
          allocate(bw_freq%A(np))
          nullify(bw_freq%next)
          nullify(bdaux)

        end if ! Empty box

        ! Allocate weight
        call safe_allocate(Wvpp,np)

        ! Store the fequency axis
        bomega%A = vpp(1:np)

        ! Store the dimension of the axis
        bomega%mfreq = np
        Redd_aux(iti,iz)%mfreq(iifreq) = np

        ! Store indexes
        bomega%ifreq = iifreq


        !
        ! Define the integration weights (same as omegabuild)
        !

        ! The first point is special in compound trapezoidal rule
        Wvpp(1) = .5d0*(vpp(2) - vpp(1))

        ! Initialize the integral to normalize the weights
        norm1 = Wvpp(1)

        ! The initial lower limit is the first point
        O0 = vpp(1)

        ! This is the pointer to the first element of the current
        ! interval, we are pointing to the first element
        cfreq = 1

        ! Flag that says that the point 2 is not the initial point of
        ! the interval (because 1 is the initial point)
        init = .False.

        ! For the rest of frequencies except the last
        do jfreq=2,np-1

          ! If ifreq is the initial point of an interval
          if (init) then

            ! The first point is special in compound trapezoidal rule
            Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq))

            ! The next point cannot be a first point
            init = .False.

            ! Initialize the integral to normalize the weights
            norm1 = Wvpp(jfreq)

            ! Pointer is now in this frequency
            cfreq = jfreq

            ! And it is the beginning of the current interval
            O0 = vpp(jfreq)

          ! If jfreq is not the initial point of an interval
          else

            ! Check if ifreq is the last point of an interval
            if (abs(vpp(jfreq+1) - vpp(jfreq)).gt.red_neglW) then

              ! The last point is special in compound trapezoidal rule
              Wvpp(jfreq) = .5d0*(vpp(jfreq) - vpp(jfreq-1))

              ! It is the end of the current interval
              O1 = vpp(jfreq)

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

              ! We know that the integral must be
              ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
              norm = 1d5*(O1 - O0)/norm1

              ! For each frequency in this interval
              do kfreq=cfreq,jfreq

                ! Normalize the weights of this interval
                Wvpp(kfreq) = Wvpp(kfreq)*norm

              end do

              ! The next point is the first point of its interval
              init = .True.

            ! If jfreq is not the last point of an interval
            else

              ! Compound trapezoidal rule weight
              Wvpp(jfreq) = .5d0*(vpp(jfreq+1) - vpp(jfreq-1))

              ! Add to the integral
              norm1 = norm1 + Wvpp(jfreq)

            endif ! Last point
          end if ! Initial point

        end do ! Frequencies

        ! The last point is special in compound trapezoidal rule
        Wvpp(np) = .5d0*(vpp(np) - vpp(np-1))

        ! It is the end of the interval
        O1 = vpp(np)

        ! Add to the integral
        norm1 = norm1 + Wvpp(np)

        ! We know that the integral must be
        ! NOTICE THE 1D5, IT IS IN PROPER cm^-1
        norm = 1d5*(O1 - O0)/norm1

        ! For each frequency in this interval
        do jfreq=cfreq,np

          ! Normalize the weights of this interval
          Wvpp(jfreq) = Wvpp(jfreq)*norm

        end do

        ! Store the weights
        bw_freq%A = Wvpp(1:np)

      end if ! Coherent scattering wing

            !
            ! Restore indentation
            !

          end do ! Output frequencies
        end do ! Output frequency ranges


        !
        ! Properly store and index the data
        !

        ! Total dimension of omega and w_freq
        nn = bomega%nback + bomega%mfreq

        ! Allocate omega and W_freq
        allocate(Redd_aux(iti,iz)%omega(nn))
        allocate(Redd_aux(iti,iz)%w_freq(nn))

        ! Determine size
        Redd_aux(iti,iz)%osize = nn

        ! Go backwards in the linked lists until done
        do while (.True.)

          ! Current index
          iifreq = bomega%ifreq

          ! Current limits
          ip = bomega%nback + 1
          ipp = ip + bomega%mfreq - 1

          ! If valid limits
          if (ipp.ge.ip) then

            ! Save vectors in boxes into auxiliar structure
            Redd_aux(iti,iz)%omega(ip:ipp) = bomega%A
            Redd_aux(iti,iz)%W_freq(ip:ipp)= bw_freq%A

          end if ! Valid limits

          ! Deallocate arrays in box
          deallocate(bomega%A,bw_freq%A)

          ! If last one
          if (.not.associated(bomega%prev)) then

            ! Clean boxes and leave
            deallocate(bomega,bw_freq)
            nullify(bomega,bw_freq)
            exit

          ! Not done with the list
          else

            ! Clean current box and go to the previous one
            bomega => bomega%prev
            bw_freq => bw_freq%prev
            nullify(bomega%next%prev,bw_freq%next%prev)
            deallocate(bomega%next,bw_freq%next)
            nullify(bomega%next,bw_freq%next)

          end if ! Last box

        end do ! Run backwards over the boxes with the frequencies

        ! End of array constructions

      end do ! Height

      end subroutine get_input_frequencies

!#####################################################################
!#####################################################################
!#####################################################################

      !> Split tasks for the calculation of the second order
      !! emissivity of a given transition\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Redd_aux(Redc_class): Auxiliar to store the input frequency
      !!                        axis for the calculation of the
      !!                        second order emissivity\n
      !!       MPID(MPI_class): Structure with MPI data\n
      !!           ia(integer): Index of current atom\n
      !!        jtran(integer): Index of current output transition\n
      !!          nti(integer): Number of input transitions for the
      !!                        current output transition\n
      !!    Lfield(logical(:)): Array flagging flags with fields if
      !!                        there is mix of yes and no fields
      subroutine setmpi_red(Atom,Red,Redd_aux,MPID,ia,jtran,nti, &
                            Lfield)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Red_class), intent(inout):: Red
      type(Redc_class), dimension(:,:), &
                        allocatable, intent(inout):: Redd_aux
      type(MPI_class), intent(in):: MPID
      logical, dimension(:), allocatable, intent(inout):: Lfield
      integer, intent(in):: ia,jtran,nti

      ! Local

      integer:: iaux,extra,i0,i1,np,nn,ntask,lnz,iz0,iz1,ff
      integer:: iz,indx,jndx,iproc,iti,jjfreq,kkfreq,iifreq,iran,ifreq
      integer, dimension(:), allocatable:: IW_freq,fB,fNB


      ! Atom -- transition index
      jndx = Red%izao(jtran,ia,Rz0)

      ! If no frequencies to calculate
      if (Red%ao(jndx)%nran.lt.1) return

      ! Number of heights
      lnz = Rz1_PRD - Rz0 + 1

      ! Number of tasks
      ntask = Red%ao(jndx)%nfreq*lnz

      ! If MPI
      if (nproc.gt.1) then

        ! If more processes than frequencies
        if (nproc.ge.ntask) then

          ! One frequency each (until done)
          iaux = 0
          do iproc=0,ntask-1
            Red%ao(jndx)%Mi0(iproc) = iproc+1
            Red%ao(jndx)%Mi1(iproc) = iproc+1
            Red%ao(jndx)%nn(iproc) = 1
          end do

          ! Rest do nothing
          do iproc=ntask,nproc-1
            Red%ao(jndx)%Mi0(iproc) = 0
            Red%ao(jndx)%Mi1(iproc) = -1
            Red%ao(jndx)%nn(iproc) = 0
          end do

        ! All participate (more tasks than CPU)
        else

          ! Only master
          if (pid.eq.0) then

            ! Basic amount
            iaux = ntask/nproc
            extra = ntask - iaux*nproc

            ! Initialize
            i0 = 0
            i1 = 0

            ! First CPUs have extra
            do iproc=0,extra-1
              Red%ao(jndx)%Mi0(iproc) = i1+1
              i1 = i1 + iaux + 1
              Red%ao(jndx)%Mi1(iproc) = i1
              Red%ao(jndx)%nn(iproc) = iaux + 1
            end do

            ! Rest have nominal
            do iproc=extra,nproc-1
              Red%ao(jndx)%Mi0(iproc) = i1+1
              i1 = i1 + iaux
              Red%ao(jndx)%Mi1(iproc) = i1
              Red%ao(jndx)%nn(iproc) = iaux
            end do

            ! Allocate and initialize work weights
            allocate(IW_freq(ntask))
            IW_freq = 1

            ! If there is information about the magnetic field
            if (allocated(Lfield)) then

              ! Allocate factors
              allocate(fB(nti),fNB(nti))

              ! For each input transition
              do iti=1,nti

                ! Calculate sizes in components
                fB(iti) = maxval(Atom%trano(jtran)%trani(iti)%indB)
                fNB(iti) = maxval(Atom%trano(jtran)%trani(iti)%indNB)

              end do ! Input transitions

              ! For each height
              do iz=Rz0,Rz1_PRD

                ! Limits
                i0 = (iz - Rz0)*Red%ao(jndx)%nfreq + 1
                i1 = i0 + Red%ao(jndx)%nfreq - 1

                ! For each input transition
                do iti=1,nti

                  ! Choose factor
                  if (Lfield(iz)) then
                    ff = fB(iti)
                  else
                    ff = fNB(iti)
                  end if

                  ! Add input frequencies to weight
                  IW_freq(i0:i1) = IW_freq(i0:i1) + &
                                   Redd_aux(iti,iz)%mfreq*ff

                end do ! Input transitions
              end do ! Heights

            ! No information about magnetic field
            else

              ! For each height
              do iz=Rz0,Rz1_PRD

                ! Limits
                i0 = (iz - Rz0)*Red%ao(jndx)%nfreq + 1
                i1 = i0 + Red%ao(jndx)%nfreq - 1

                ! For each input transition
                do iti=1,nti

                  ! Add input frequencies to weight
                  IW_freq(i0:i1) = IW_freq(i0:i1) + &
                                   Redd_aux(iti,iz)%mfreq

                end do ! Input transitions
              end do ! Heights

            end if ! Magnetic field information

            ! Control overflow
            do while (sum(IW_freq).lt.maxval(IW_freq).and. &
                      minval(IW_freq).gt.1)

              ! Reduce by factor 2
              IW_freq = IW_freq/2

            end do ! If overflow

            ! Optimize split
            call weighted_split(nproc,1000, &
                                IW_freq, &
                                Red%ao(jndx)%Mi0, &
                                Red%ao(jndx)%Mi1, &
                                Red%ao(jndx)%nn)

            ! Deallocate weights
            deallocate(IW_freq)

          end if ! Master

          ! Share split with everyone
          call MPI_BCAST(Red%ao(jndx)%nn(0), nproc, &
                         MPI_INTEGER, 0, MPI_COMM_RT, ierr)
          call MPI_BCAST(Red%ao(jndx)%Mi0(0), nproc, &
                         MPI_INTEGER, 0, MPI_COMM_RT, ierr)
          call MPI_BCAST(Red%ao(jndx)%Mi1(0), nproc, &
                         MPI_INTEGER, 0, MPI_COMM_RT, ierr)

        end if ! Type of distribution

      ! Serial
      else

        ! Everything for single CPU
        Red%ao(jndx)%Mi0(0) = 1
        Red%ao(jndx)%Mi1(0) = ntask
        Red%ao(jndx)%nn(0) = ntask

      end if ! MPI/serial

      ! If CPU has tasks
      if (Red%ao(jndx)%nn(pid).gt.0) then

        ! Get z limits
        iz0 = Rz0 + (Red%ao(jndx)%Mi0(pid)-1)/Red%ao(jndx)%nfreq
        iz1 = Rz0 + (Red%ao(jndx)%Mi1(pid)-1)/Red%ao(jndx)%nfreq

        ! For each height
        do iz=iz0,iz1

          ! Get index
          indx = Red%izao(jtran,ia,iz)

          ! Get left frequency indexes
          if (iz.eq.iz0) then
            i0 = Red%ao(jndx)%Mi0(pid) - (iz0-Rz0)*Red%ao(jndx)%nfreq
          else
            i0 = 1
          end if

          ! Get right frequency index
          if (iz.eq.iz1) then
            i1 = Red%ao(jndx)%Mi1(pid) - (iz1-Rz0)*Red%ao(jndx)%nfreq
          else
            i1 = Red%ao(jndx)%nfreq
          end if

          ! Reinitialize indexes
          Red%zao(indx)%Igf0 = nfreq+1
          Red%zao(indx)%Igf1 = 0

          !
          ! Now take only what is needed
          !

          ! Initialize maximum number of frequencies for this CPU
          Red%zao(indx)%mxfreq = 0

          ! For each input transition
          do iti=1,nti

            ! Initialize running index
            jjfreq = 0
            kkfreq = 0

            ! Allocate local sizes
            allocate(Red%zao(indx)%trani(iti)%mfreq(i0:i1))

            ! Copy number of frequencies
            Red%zao(indx)%trani(iti)%mfreq = &
                                         Redd_aux(iti,iz)%mfreq(i0:i1)

            ! Total number of frequencies for this CPU
            np = sum(Red%zao(indx)%trani(iti)%mfreq)

            ! Maximum
            Red%zao(indx)%mxfreq = max(Red%zao(indx)%mxfreq, &
                                       maxval(Red%zao(indx)% &
                                                  trani(iti)%mfreq))
            ! Size
            Red%zao(indx)%trani(iti)%osize = np

            ! Allocate local arrays
            allocate(Red%zao(indx)%trani(iti)%omega(np))
            allocate(Red%zao(indx)%trani(iti)%W_freq(np))

            ! For each output frequency
            iifreq = 0
            do iran=1,Red%ao(jndx)%nran
              do ifreq=Red%ao(jndx)%if0(iran), &
                       Red%ao(jndx)%if1(iran)

                ! Advance index
                iifreq = iifreq + 1

                ! If out of range above, done
                if (iifreq.gt.i1) exit

                ! If out of range below
                if (iifreq.lt.i0) then

                  ! Advance and skip
                  jjfreq = jjfreq + Redd_aux(iti,iz)%mfreq(iifreq)
                  cycle

                end if ! Out of range below

                ! Skip if 0 size
                if (Redd_aux(iti,iz)%mfreq(iifreq).lt.1) cycle

                ! Copy omega
                Red%zao(indx)%trani(iti)% &
                    omega(kkfreq+1: &
                          kkfreq+Redd_aux(iti,iz)%mfreq(iifreq)) = &
                  Redd_aux(iti,iz)%&
                    omega(jjfreq+1: &
                          jjfreq+Redd_aux(iti,iz)%mfreq(iifreq))

                ! Copy weight
                Red%zao(indx)%trani(iti)% &
                    W_freq(kkfreq+1: &
                           kkfreq+Redd_aux(iti,iz)%mfreq(iifreq)) = &
                  Redd_aux(iti,iz)%&
                    W_freq(jjfreq+1: &
                           jjfreq+Redd_aux(iti,iz)%mfreq(iifreq))

                ! Advance indexes
                jjfreq = jjfreq + Redd_aux(iti,iz)%mfreq(iifreq)
                kkfreq = kkfreq + Redd_aux(iti,iz)%mfreq(iifreq)

              end do ! Output frequencies
            end do ! Output ranges

            ! Get size of omega to add
            nn = np

          end do ! Input transitions

          ! For each output frequency
          iifreq = 0
          do iran=1,Red%ao(jndx)%nran
            do ifreq=Red%ao(jndx)%if0(iran), &
                     Red%ao(jndx)%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! If out of range
              if (iifreq.gt.i1) exit
              if (iifreq.lt.i0) cycle

              ! Update limits
              if (ifreq.lt.Red%zao(indx)%Igf0) &
                Red%zao(indx)%Igf0 = ifreq
              if (ifreq.gt.Red%zao(indx)%Igf1) &
                Red%zao(indx)%Igf1 = ifreq

            end do
          end do
        end do ! Height

        ! Determine frequency range to save
        Red%ao(jndx)%gf0 = max(Red%ao(jndx)%tgf0,MPID%if0(pid))
        Red%ao(jndx)%gf1 = min(Red%ao(jndx)%tgf1,MPID%if1(pid))

      end if ! Tasks to deal with

      ! For every considered height
      do iz=Rz0,Rz1_PRD

        ! For every input transition
        do iti=1,nti

          ! Free space in auxiliar structure
          deallocate(Redd_aux(iti,iz)%omega, &
                     Redd_aux(iti,iz)%W_freq, &
                     Redd_aux(iti,iz)%mfreq)

        end do ! Input transitions
      end do ! Heights

      ! Free auxiliar
      deallocate(Redd_aux)

      end subroutine setmpi_red

!#####################################################################
!#####################################################################
!#####################################################################

      !> Find the frequency index limits for the calculation of the
      !! second order emissivity once the input frequency axes are
      !! known\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Red(Red_class): Structure with redistribution input
      !!                    frequency data, redistribution function
      !!                    data, and profile or normalization data\n
      !!  omega(double(:)): Frequency axis\n
      !!       ia(integer): Index for current atom\n
      !!      pol(logical): If polarized case
      subroutine find_integral_limits(Atom,Atmo,Red,omega,ia,pol)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Red_class), intent(inout):: Red
      logical, intent(in):: pol
      integer, intent(in):: ia
      double precision, dimension(:), intent(in):: omega


      ! Local

      integer:: iz0,iz1,jndx,i0,i1
      integer:: iz,ktran,jtran,fjtran,indx,jtrano,itran,fitran,ffitran
      integer:: itermu,itermf,iterml,iJu,iJf,iJl,ii,iti
      integer:: jjfreq,iifreq,iran,ifreq,lifreq,jfreq,jufreq,ibfreq
      integer:: tif0,tif1,np

      double precision:: vfacm,vfacp,dnl,dnlmin,dnlmax,O0,O1

      ! For each output transition level
      do ktran=lbound(Red%izao,1),ubound(Red%izao,1)

        ! If polarization
        if (pol) then

          ! Skip if too big index
          if (ktran.gt.Atom%ntran) cycle

          ! Get transition term index
          jtran = ktran

          ! If no PRD line, skip
          if (.not.Atom%lemiss2(jtran)) cycle

          ! Other index
          jtrano = jtran

        ! If intensity
        else

          ! Skip if too big index
          if (ktran.gt.Atom%nftran) cycle

          ! Get transition term index
          jtran = Atom%ifst(ktran)

          ! If no PRD line, skip
          if (.not.Atom%lemiss2(jtran)) cycle

          ! Other index
          jtrano = Atom%itrano(ktran)

        end if

        ! Atom -- transition index
        jndx = Red%izao(ktran,ia,Rz0)

        ! Skip if no assigned frequencies
        if (Red%ao(jndx)%nran.lt.1) cycle
        if (Red%ao(jndx)%nn(pid).lt.1) cycle

        ! Get z limits
        iz0 = Rz0 + (Red%ao(jndx)%Mi0(pid)-1)/Red%ao(jndx)%nfreq
        iz1 = Rz0 + (Red%ao(jndx)%Mi1(pid)-1)/Red%ao(jndx)%nfreq

        ! For each height
        do iz=iz0,iz1

          ! Initialize Doppler displacement factor
          vfacm = 1d0
          vfacp = 1d0

          ! If dynamic
          if (dyn) then

            ! Compute velocity
            dnl = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))

            ! Big enough velocity
            if (dnl.gt.TINYVEL) then

              ! Compute displacements
              vfacp = 1d0 - dnl
              vfacp = 1d0/vfacp
              vfacm = 1d0 + dnl
              vfacm = 1d0/vfacm

            end if ! Significant velocity
          end if ! Dynamic

          ! Get index
          indx = Red%izao(ktran,ia,iz)

          ! Get left frequency indexes
          if (iz.eq.iz0) then
            i0 = Red%ao(jndx)%Mi0(pid) - (iz0-Rz0)*Red%ao(jndx)%nfreq
          else
            i0 = 1
          end if

          ! Get right frequency index
          if (iz.eq.iz1) then
            i1 = Red%ao(jndx)%Mi1(pid) - (iz1-Rz0)*Red%ao(jndx)%nfreq
          else
            i1 = Red%ao(jndx)%nfreq
          end if

          ! If intensity case
          if (.not.pol) then

            ! For every transition component
            do ii=1,Atom%fst(jtran)%nt

              ! If rolling index corresponds to this one
              if (Atom%ifst_ij(ii,jtran).eq.ktran) then

                ! Set internal index
                fjtran = ii
                exit

              end if ! Same rolling index

            end do ! Transition components

            ! Get level indexes
            iJu = Atom%fst(jtran)%ilevelu(fjtran)
            iJf = Atom%fst(jtran)%ilevell(fjtran)

          end if ! Polarization/intensity

          ! Get term indexes
          itermu = Atom%fst(jtran)%itermu
          itermf = Atom%fst(jtran)%iterml

          ! For every input transition
          do iti=1,size(Red%zao(indx)%trani)

            ! If polarized case
            if (pol) then

              ! Get transition and term indexes
              itran = Atom%trano(jtran)%indT(iti)
              iterml = Atom%fst(itran)%iterml

              ! Resonance
              dnlmax = -1d99
              dnlmin =  1d99

              ! For each level in the lower term
              do iJl=1,Atom%nJ(iterml)

                ! For each level in the final term
                do iJf=1,Atom%nJ(itermf)

                  ! Raman resonance
                  dnl = Atom%FSfreq(iJf,itermf) - &
                        Atom%FSfreq(iJl,iterml)

                  ! Check extrema
                  if (dnl.gt.dnlmax) dnlmax = dnl
                  if (dnl.lt.dnlmin) dnlmin = dnl

                end do ! Final levels
              end do ! Lower levels

            ! Intensity case
            else

              ! Get indexes for term and transition
              ffitran = Atom%tranoI(jtrano)%indT(iti)
              itran = Atom%ifst(ffitran)
              fitran = Atom%ifstj(ffitran)
              iterml = Atom%fst(itran)%iterml

              ! Get J index
              iJl = Atom%fst(itran)%ilevell(fitran)

              ! Resonance
              dnlmin = Atom%FSfreq(iJf,itermf) - &
                       Atom%FSfreq(iJl,iterml)
              dnlmax = dnlmin

            end if ! Polarization/intensity

            ! Get atomic limits
            tif0 = Atom%tif0(itran)
            tif1 = Atom%tif1(itran)


            !
            ! Reset identation
            !

      !
      ! Initialize input frequency index
      jjfreq = 0

      ! For each output frequency
      iifreq = 0
      do iran=1,Red%ao(jndx)%nran
        do ifreq=Red%ao(jndx)%if0(iran), &
                 Red%ao(jndx)%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! In MPI range
          if (iifreq.gt.i1) exit
          if (iifreq.lt.i0) cycle

          ! Input frequency number
          np = Red%zao(indx)%trani(iti)%mfreq(iifreq)

          ! If empty
          if (np.lt.1) then

            ! Left and right limits from resonance
            lifreq = Atom%tif0(itran)
            jfreq = Atom%tif1(itran)
            O0 = omega(ifreq) + dnlmin
            O1 = omega(ifreq) + dnlmax

            ! Correct
            if (O0.lt.omega(lifreq)) O0 = omega(lifreq)
            if (O1.gt.omega(jfreq)) O1 = omega(jfreq)

            !
            ! Look for the indexes

            !
            ! Left

            ! Only if not beyond already
            if (omega(lifreq)*vfacm.lt.O0) then

              ! Search
              do while (.True.)

                ! Check if next inside
                if (omega(lifreq+1)*vfacm.gt.O0) exit

                ! Advance
                lifreq = lifreq + 1
                if ((lifreq+1).gt.nfreq) exit

              end do ! Search

            end if ! Need to search

            !
            ! Right

            ! Only if not beyond already
            if (omega(jfreq)*vfacp.gt.O1) then

              ! Search
              do while (.True.)

                ! Check if already inside
                if (omega(jfreq-1)*vfacp.lt.O1) exit

                ! Advance
                jfreq = jfreq - 1
                if ((jfreq-1).lt.1) exit

              end do

            end if ! Need to search

            ! Update global limits
            if (lifreq.lt.Red%zao(indx)%ggf0) &
              Red%zao(indx)%ggf0 = lifreq
            if (jfreq.gt.Red%zao(indx)%ggf1) &
              Red%zao(indx)%ggf1 = jfreq

            ! Skip rest
            cycle

          end if ! Empty axis, only resonance

          ! Add to size
          Red%zao(indx)%trani(iti)%isize = np + &
                                        Red%zao(indx)%trani(iti)%isize

          !
          ! Reset identation
          !

      ! Compute jump to check only extrema
      lifreq = Atom%tif0(itran)
      jufreq = np
      if (np.gt.1) jufreq = jufreq - 1

      ! First and last frequencies
      do jfreq=jjfreq+1,jjfreq+np,jufreq

        ! If out of range
        if (Red%zao(indx)%trani(iti)%omega(jfreq)*vfacm.le. &
            (omega(tif0)+TINYO)) then

          ! We are in the first frequency, take boundary value
          if (Red%zao(indx)%ggf0.gt.tif0) &
            Red%zao(indx)%ggf0 = tif0
          if (Red%zao(indx)%ggf1.lt.tif0) &
            Red%zao(indx)%ggf1 = tif0

        ! If out of range
        else if (Red%zao(indx)%trani(iti)%omega(jfreq)*vfacp.ge. &
                 (omega(tif1) - TINYO)) then

          ! We are in the last frequency, take boundary value
          if (Red%zao(indx)%ggf0.gt.tif1) &
            Red%zao(indx)%ggf0 = tif1
          if (Red%zao(indx)%ggf1.lt.tif1) &
            Red%zao(indx)%ggf1 = tif1

        ! If within the boundaries
        else

          ! Search between the last found frequency and
          ! all but the boundary for left displacement
          do ibfreq=lifreq,tif1-1

            ! If the input is between this output and the next
            if (Red%zao(indx)%trani(iti)%omega(jfreq)*vfacm.ge. &
                                         omega(ibfreq).and. &
                Red%zao(indx)%trani(iti)%omega(jfreq)*vfacm.le. &
                                         omega(ibfreq+1)) then

              ! Update lifreq
              lifreq = ibfreq

              ! Found frequency, update limits
              if (Red%zao(indx)%ggf0.gt.ibfreq) &
                Red%zao(indx)%ggf0 = ibfreq
              if (Red%zao(indx)%ggf1.lt.ibfreq+1) &
                Red%zao(indx)%ggf1 = ibfreq+1

              exit

            end if ! Check output frequency

          end do ! Run output frequencies

          ! Search between the last found frequency and
          ! all but the boundary for right displacement
          do ibfreq=lifreq,tif1-1

            ! If the input is between this output and the next
            if (Red%zao(indx)%trani(iti)%omega(jfreq)*vfacp.ge. &
                                         omega(ibfreq).and. &
                Red%zao(indx)%trani(iti)%omega(jfreq)*vfacp.le. &
                                         omega(ibfreq+1)) then

              ! Found frequency
              if (Red%zao(indx)%ggf0.gt.ibfreq) &
                Red%zao(indx)%ggf0 = ibfreq
              if (Red%zao(indx)%ggf1.lt.ibfreq+1) &
                Red%zao(indx)%ggf1 = ibfreq+1

              exit

            end if ! Check output frequency

          end do ! Run output frequencies

        end if ! Check if out of limits

      end do ! Run input frequencies

      ! Advance index
      jjfreq = jjfreq + np

          !
          ! Restore indentation
          !


          end do ! Output frequencies
        end do ! Output frequency ranges

            !
            ! Restore identation
            !

          end do ! Input transition
        end do ! Output transition
      end do ! height nodes

      end subroutine find_integral_limits

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate a double precision array. If the array is allocated,
      !! free the space unless it has the correct dimension already\n
      !!  array(double(:)): Array to allocate\n
      !!      siz(integer): Size to allocate
      subroutine safe_allocate(array,siz)

      ! I/O

      integer, intent(in):: siz
      double precision, dimension(:), allocatable:: array


      ! If not allocated
      if (.not.allocated(array)) then

        ! Allocate
        allocate(array(siz))

      ! Already allocated
      else

        ! If smaller than requested
        if (size(array).lt.siz) then

          ! Make bigger
          deallocate(array)
          allocate(array(siz))

        end if ! Smaller than requested
      end if ! Allocated or not

      end subroutine safe_allocate

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine the input frequency axis for the calculation of the
      !! second order emissivity and split the tasks of the emissivity
      !! calculation as evenly as possible. Assume comoving reference
      !! frame\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Bstrength(double(:)): Magnetic field strength\n
      !!     Input(Input_class): Structure with configuration data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!         ofram(logical): If reached the RAM limit
      subroutine omegabuildin(Frec,Red,Atom,Atmo,Bstrength,Input, &
                              MPID,ofram)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(MPI_class), intent(in):: MPID
      double precision, dimension(:), intent(in):: Bstrength
      logical, intent(out):: ofram

      ! Local

      type(Redc_class), dimension(:,:), allocatable:: Redd_aux

      logical:: skip,Yfield,Nfield,RAMOF,ldipole
      logical, dimension(:), allocatable:: Lfield

      integer:: iz,ia,itran,jtran,if0,if1
      integer:: iJl,iJu,iJf,itermf,itermu,iterml
      integer:: nr,nt,ntj,ntk,ni,np0,nti,iti
      integer, dimension(:), allocatable:: flag

      double precision:: dnlmin,dnlmax,SRAM,DwN
      double precision, dimension(:), allocatable:: dnl, nut
      double precision, dimension(:), allocatable:: vphv, vplv, vpr
      double precision, dimension(:), allocatable:: vphve, vplve
      double precision, dimension(:), allocatable:: vpp


      ! Routine name
      urou = 'omegabuildin'

      ! Initialize out of RAM
      ofram = .False.

      ! Flags
      Yfield = .False.
      Nfield = .False.
      RAMOF = .False.

      ! Small field
      Nfield = minval(Bstrength(Rz0:Rz1_PRD)).le.TINYB
      Yfield = maxval(Bstrength(Rz0:Rz1_PRD)).gt.TINYB

      ! If found both
      if (Yfield.and.Nfield) then

        ! Allocate vector
        allocate(Lfield(Rz0:Rz1_PRD))

        ! For each considered height
        do iz=Rz0,Rz1_PRD

          ! Flag
          Lfield(iz) = Bstrength(iz).gt.TINYB

        end do ! Heights

      end if ! Found both magnetic and non-magnetic

      ! Flag if we need to check selection rules
      ldipole = .not.Yfield.and.Input%MIT_input.lt.0

      ! Preliminar allocation of vpp, auxiliar vector to store
      ! frequencies
      np0 = 10000
      allocate(vpp(np0))
      allocate(flag(np0))

      ! For each atom
      do ia=1,nA

        ! Initialize as skipping
        skip = .True.

        ! For every transition
        do jtran=1,Atom(ia)%ntran

          ! If it is PRD
          if (Atom(ia)%lemiss2(jtran)) then

            ! Flag as not skipping and stop
            skip = .False.
            exit

          end if ! PRD transition

        end do ! Transitions

        ! There are no PRD lines for this atom
        if (skip) cycle

        !
        ! Calculate nominal Doppler width
        !

        ! Maximum temperature
        if (Input%dws.eq.'MAX') then

          ! Calculate
          DwN = Atom(ia)%cDopp*sqrt(Input%maxT)

        ! Minimum temperature
        else if (Input%dws.eq.'MIN') then

          ! Calculate
          DwN = Atom(ia)%cDopp*sqrt(Input%minT)

        ! Fixed input Doppler width
        else if (Input%dws.eq.'NUM') then

          ! Calculate
          DwN = Input%dw*1d-9/c

        end if ! Type of nominal Doppler width

        ! For each output transition
        do jtran=1,Atom(ia)%ntran

          ! If no PRD line, skip
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get terms
          itermf = Atom(ia)%fst(jtran)%iterml
          itermu = Atom(ia)%fst(jtran)%itermu

          ! Master
          if (pid.eq.0) then

            ! Get limits
            if0 = Atom(ia)%if0(jtran)
            if1 = Atom(ia)%if1(jtran)

          end if

          ! If MPI
          if (nproc.gt.1) then

            ! Share limits
            call MPI_BCAST(if0,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)
            call MPI_BCAST(if1,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)

          end if

          ! Number of input transitions
          nti = Atom(ia)%trano(jtran)%nt

          ! Count input transitions
          ntk = 0

          ! For each input transition
          do iti=1,Atom(ia)%trano(jtran)%nt

            ! Get transition index and lower term
            itran = Atom(ia)%trano(jtran)%indT(iti)
            iterml = Atom(ia)%fst(itran)%iterml

            !
            ! Count input components
            !

            ! For each level in upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If checking dipole rules, do it
                if (ldipole.and. &
                    (abs(Atom(ia)%rJval(iJu,itermu) - &
                         Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                     Atom(ia)%rJval(iJu,itermu) + &
                     Atom(ia)%rJval(iJl,iterml).lt..4d0)) cycle

                ! Add to count
                ntk = ntk + 1

              end do ! Lower levels
            end do ! Upper levels
          end do ! Input transitions


          !
          ! Count transitions
          !

          ! Reset index of input transition
          iti = 0
          nt = 0

          !
          ! Count the FS components in the output transition
          !

          ! For each level in the upper term
          do iJu=1,Atom(ia)%nJ(itermu)

            ! For each level in the final term
            do iJf=1,Atom(ia)%nJ(itermf)

              ! If checking dipole rules, do it
              if(ldipole.and. &
                 (abs(Atom(ia)%rJval(iJu,itermu) - &
                      Atom(ia)%rJval(iJf,itermf)).gt.1.or. &
                  Atom(ia)%rJval(iJu,itermu) + &
                  Atom(ia)%rJval(iJf,itermf).lt..4d0)) cycle

              ! Add to count
              nt = nt + 1

            end do ! Final levels
          end do ! Upper levels

          ! Allocate component frequency array
          call safe_allocate(nut,nt+ntk)


          !
          ! Save transitions
          !

          ! Reset index
          nt = 0

          ! For each level in the upper term
          do iJu=1,Atom(ia)%nJ(itermu)

            ! For each level in the final term
            do iJf=1,Atom(ia)%nJ(itermf)

              ! If checking dipole rules, do it
              if(ldipole.and. &
                 (abs(Atom(ia)%rJval(iJu,itermu) - &
                      Atom(ia)%rJval(iJf,itermf)).gt.1.or. &
                  Atom(ia)%rJval(iJu,itermu) + &
                  Atom(ia)%rJval(iJf,itermf).lt..4d0)) cycle

              ! Advance count and save frequency
              nt = nt + 1
              nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                        Atom(ia)%FSfreq(iJf,itermf)

            end do ! Final levels
          end do ! Upper levels

          ! Save number of output transition components
          ntj = nt

          ! For each other transition that shares upper term
          do iti=1,Atom(ia)%trano(jtran)%nt

            ! Get transition and term indexes
            itran = Atom(ia)%trano(jtran)%indT(iti)
            iterml = Atom(ia)%fst(itran)%iterml

            ! For each level in the upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in the lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If checking dipole rules, do it
                if(ldipole.and. &
                   (abs(Atom(ia)%rJval(iJu,itermu) - &
                        Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                    Atom(ia)%rJval(iJu,itermu) + &
                    Atom(ia)%rJval(iJl,iterml).lt..4d0)) cycle

                  ! Advance count and save frequency
                  nt = nt + 1
                  nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                            Atom(ia)%FSfreq(iJl,iterml)

              end do ! Lower levels
            end do ! Upper levels
          end do ! Lower transitions

          !
          ! Determine frequency limits for the output transition
          call get_transition_out_limit(Atmo,Red,ia, &
                                        jtran,nt,if0,if1, &
                                        Input%red_pars(3), &
                                        DwN,nut,Frec%omega)

          ! Allocate auxiliar structure to store input frequencies
          allocate(Redd_aux(nti,Rz0:Rz1_PRD))

          ! For each input transition
          do iti=1,Atom(ia)%trano(jtran)%nt

            ! Get transition and term indexes
            itran = Atom(ia)%trano(jtran)%indT(iti)
            iterml = Atom(ia)%fst(itran)%iterml

            !
            ! Count the number of Raman resonances. We have as many
            ! resonances as combinations between lower levels (the
            ! repeated ones will be removes later)
            !

            ! Reset the index
            nr = 0

            ! For each level in lower term
            do iJl=1,Atom(ia)%nJ(iterml)

              ! For each level in final term
              do iJf=1,Atom(ia)%nJ(itermf)

                ! Add resonance
                nr = nr + 1

              end do ! Final levels
            end do ! Lower levels

            ! Allocate space for resonances
            call safe_allocate(dnl,nr)

            ! Reset number of resonances
            nr = 0

            ! Reset ranges
            dnlmax = -1D99
            dnlmin = 1D99

            ! For each level in lower term
            do iJl=1,Atom(ia)%nJ(iterml)

              ! For each level in final term
              do iJf=1,Atom(ia)%nJ(itermf)

                ! Add the resonance and store frequency
                nr = nr + 1
                dnl(nr) = Atom(ia)%FSfreq(iJf,itermf) - &
                          Atom(ia)%FSfreq(iJl,iterml)

                ! And check that we know what are the futhest ones
                if (dnl(nr).gt.dnlmax) dnlmax = dnl(nr)
                if (dnl(nr).lt.dnlmin) dnlmin = dnl(nr)

              end do ! Lower levels
            end do ! Final levels

            ! Reset to number of output transitions
            nt = ntj

            ! For each level in upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If checking dipole rules, do it
                if (ldipole.and. &
                    (abs(Atom(ia)%rJval(iJu,itermu) - &
                         Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                     Atom(ia)%rJval(iJu,itermu) + &
                     Atom(ia)%rJval(iJl,iterml).lt..4d0)) cycle

                ! Advance count
                nt = nt + 1

              end do ! Lower levels
            end do ! Upper levels

            ! The total number of limits to consider is the sum of
            ! lines and resonances
            ni = nr + nt

            !
            ! Allocate
            !
            ! Upper ranges
            call safe_allocate(vphv,ni)
            ! Lower ranges
            call safe_allocate(vplv,ni)
            ! Extended upper range
            call safe_allocate(vphve,ni*2)
            ! Extended lower range
            call safe_allocate(vplve,ni*2)
            ! Resonances
            call safe_allocate(vpr,nr)


            !
            ! Store the frequencies of the FS transitions
            !

            ! Reset index
            nt = 0

            ! For each level in the upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in the final term
              do iJf=1,Atom(ia)%nJ(itermf)

                ! Skip PRD MIT outputs if not electric dipole
                if (.True..and. &
                    (abs(Atom(ia)%rJval(iJu,itermu) - &
                         Atom(ia)%rJval(iJf,itermf)).gt.1.or. &
                     Atom(ia)%rJval(iJu,itermu) + &
                     Atom(ia)%rJval(iJf,itermf).lt..4d0)) cycle

                ! Advance index and store frequency
                nt = nt + 1
                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJf,itermf)

              end do ! Final levels
            end do ! Upper levels

            ! For each level in the upper term
            do iJu=1,Atom(ia)%nJ(itermu)

              ! For each level in the lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! If checking dipole rules, do it
                if(ldipole.and. &
                   (abs(Atom(ia)%rJval(iJu,itermu) - &
                        Atom(ia)%rJval(iJl,iterml)).gt.1.or. &
                    Atom(ia)%rJval(iJu,itermu) + &
                    Atom(ia)%rJval(iJl,iterml).lt..4d0)) cycle

                ! Advance index and store frequency
                nt = nt + 1
                nut(nt) = Atom(ia)%FSfreq(iJu,itermu) - &
                          Atom(ia)%FSfreq(iJl,iterml)

              end do ! Lower levels
            end do ! Upper levels

            ! Get input frequency axis and its weights
            call get_input_frequencies(Atmo,Red,Redd_aux, &
                                       Input%cohw, &
                                       ia,jtran,iti,nr,ntj,nt,np0, &
                                       Atom(ia)%Dfreq(itran), &
                                       Atom(ia)%Dfreq(jtran), &
                                       Input%dcohw,Input%red_pars, &
                                       Frec%omega,dnl,nut, &
                                       Atom(ia)%cDopp, &
                                       vphv,vplv,vphve,vplve,vpr, &
                                       vpp,flag)

          end do ! Input transition

          ! Split tasks by amount of frequencies
          call setmpi_red(Atom(ia),Red,Redd_aux,MPID,ia,jtran, &
                          Atom(ia)%trano(jtran)%nt,Lfield)

        end do ! Output transition

        ! Find frequency limits for the frequency range needed for
        ! integrals
        call find_integral_limits(Atom(ia),Atmo,Red, &
                                  Frec%omega,ia,.True.)

      end do ! Atoms

      ! Count memory
      call cram_red_frec(Red,SRAM)
      FRAMc = SRAM

      ! Check if everything is fine
      call control

      return

      end subroutine omegabuildin

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate space to store the redistribution functions\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!  Bstrength(double(:)): Magnetic field strength\n
      !!        ofram(logical): If reached the RAM limit
      subroutine allocate_Warr(Atom,Red,Geom,Bstrength,ofram)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Red_class), intent(inout):: Red
      type(Geometry_class), intent(inout):: Geom
      double precision, dimension(:), intent(in):: Bstrength
      logical, intent(out):: ofram

      ! Local

      logical:: Yfield,Nfield,lNCHLT

      integer:: indx,jndx,iz0,iz1
      integer:: iz,ia,jtran,iti,itran
      integer:: nn,iYNF,iYYF,iNF,iDF,iDFR

      double precision:: RAM,SRAM,dnn


      ! Initialize RAM counters
      WRAMc = 0d0

      ! Initialize
      ofram = .False.

      ! Generate scattering angles if angle-dependent
      if (.not.AV) &
        call get_scattering(Geom)

      ! If not storing, go back
      if (.not.PRAM) return

      ! Get current RAM status
      RAM = cram_add(1)

      ! Allocate rzao
      allocate(Red%rzao(Red%nzao))

      ! For each index
      do indx=1,Red%nzao

        ! Add size of logical
        RAM = RAM + 4d-6*dble(size(Red%zao(indx)%trani))
        nullify(Red%rzao(indx)%trani)

      end do ! Indexes

      ! Flags
      Yfield = .False.
      Nfield = .False.

      ! For each considered height
      do iz=Rz0,Rz1_PRD

        ! If too small field
        if (Bstrength(iz).lt.TINYB) then

          ! Flag there is a no field height
          Nfield = .True.

        ! Significant enough field
        else

          ! Flag there is a yes field height
          Yfield = .True.

        end if ! Field strength

        ! If found both magnetic and unmagnetized heights, stop
        if (Yfield.and.Nfield) exit

      end do ! Heights

      ! For each atom
      do ia=1,nA

        ! For each output transition
        do jtran=1,Atom(ia)%ntran

          ! Skip CRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get indexes
          jndx = Red%izao(jtran,ia,Rz0)

          ! Skip no work
          if (Red%ao(jndx)%nran.lt.1) cycle
          if (Red%ao(jndx)%nn(pid).lt.1) cycle

          ! Get z limits
          iz0 = Rz0 + (Red%ao(jndx)%Mi0(pid)-1)/Red%ao(jndx)%nfreq
          iz1 = Rz0 + (Red%ao(jndx)%Mi1(pid)-1)/Red%ao(jndx)%nfreq

          ! For every height
          do iz=iz0,iz1

            ! Get indexes
            indx = Red%izao(jtran,ia,iz)

            ! Allocate
            allocate(Red%rzao(indx)%trani(Atom(ia)%trano(jtran)%nt))

            ! For each input transition
            do iti=1,Atom(ia)%trano(jtran)%nt

              ! Add RAM
              RAM = RAM + 1d-6*sizeof(Red%rzao(indx)%trani(iti))

            end do ! Input transitions
          end do ! Heights

          ! For each input transition
          do iti=1,Atom(ia)%trano(jtran)%nt

            ! Get sizes in terms of components
            if (Yfield) then
              iYNF = maxval(Atom(ia)%trano(jtran)%trani(iti)%indB)
              iYYF = iYNF - Atom(ia)%trano(jtran)%trani(iti)%nchlt
            end if
            if (Nfield) &
              iNF = maxval(Atom(ia)%trano(jtran)%trani(iti)%indNB)

            ! Get transition index
            itran = Atom(ia)%trano(jtran)%indT(iti)

            ! Assume coherent
            lNCHLT = .False.

            ! If non-coherent input and magnetic field
            if (Yfield.and.NCHLT) then

              ! For each height
              do iz=iz0,iz1

                ! If non-coherent for this transition
                if (Atom(ia)%NCHLT(iz,itran)) then

                  ! Flag and leave
                  lNCHLT = .True.
                  exit

                end if ! Non-coherent for this transition

              end do ! Heights

            end if ! Magnetic and non-coherent input

            ! For every considered height
            do iz=iz0,iz1

              ! Get indexes
              indx = Red%izao(jtran,ia,iz)

              ! Initialize as not stored
              Red%rzao(indx)%trani(iti)%RAM = .False.

              ! Non magnetic
              if (Bstrength(iz).le.TINYB) then

                ! Get component sizes
                iDF = iNF
                iDFR = iNF

              ! Magnetic
              else

                ! Get component sizes
                iDF = iYNF

                ! Check if NCHLT
                if (lNCHLT) then

                  ! If NCHLT applies
                  if (Atom(ia)%NCHLT(iz,itran)) then

                    ! Get reduces size
                    iDFR = iYYF

                  ! NCHLT does not apply
                  else

                    ! Get full size
                    iDFR = iYNF

                  end if ! NCHLT

                ! No NCHLT
                else

                  ! Just full size
                  iDFR = iYNF

                end if ! NHCLT
              end if ! Magnetic field or not

              ! Predict size of next block
              nn = sum(Red%zao(indx)%trani(iti)%mfreq)*iDFR
              dnn = dble(sum(Red%zao(indx)%trani(iti)%mfreq))* &
                    dble(iDFR)

              ! Skip if no size
              if (nn.le.0) cycle

              ! Skip if overflown
              if (dnn - dble(nn).gt.0.9d0) cycle

              ! Correct if AD
              if (.not.AV) then

                ! Skip coherent one
                if (jtran.eq.itran) then

                  ! Scattering angles to consider
                  nn = nn*(Geom%nScatt-1)
                  dnn = dnn*dble(Geom%nScatt-1)

                ! Not Rayleigh
                else

                  ! Full scattering angles
                  nn = nn*Geom%nScatt
                  dnn = dnn*dble(Geom%nScatt)

                end if ! Rayleigh scattering
              end if ! Angle-dependent

              ! Skip if overflown
              if (dnn - dble(nn).gt.0.9d0) cycle

              ! Get size in MB (single precision) + logical
              SRAM = 8d-6*dble(nn) + 4d-6*iDF

              ! If no more space
              if (floor(RAM+SRAM).gt.RLIM.or.SRAM.le.0d0) then

                ! Full memory, do not allocate
                ofram = .True.
                Red%rzao(indx)%trani(iti)%RAM = .False.
                cycle

              end if

              ! Add to RAM count
              RAM = RAM + SRAM

              ! Saving
              Red%rzao(indx)%trani(iti)%RAM = .True.

              ! Allocate and initialize initialization flag
              allocate(Red%rzao(indx)%trani(iti)%iPPRD(iDF))
              Red%rzao(indx)%trani(iti)%iPPRD = .True.

              ! Allocate space
              allocate(Red%rzao(indx)%trani(iti)%PWarr2(nn))

            end do ! heights
          end do ! Input transitions
        end do ! Output transitions
      end do ! Atoms

      ! Count RAM in redistribution
      call cram_red_warr(Red,RAM)
      WRAMc = RAM

      end subroutine allocate_Warr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Determine the input frequency axis for the calculation of the
      !! second order intensity emissivity and split the tasks of the
      !! intensity emissivity calculation as evenly as possible.
      !! Assume comoving reference frame\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!     Input(Input_class): Structure with configuration data\n
      !!        MPID(MPI_class): Structure with MPI data\n
      !!         ofram(logical): If reached the RAM limit
      subroutine omegabuildinI(Frec,Red,Atom,Atmo,Input,MPID,ofram)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(Input_class), intent(in):: Input
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(MPI_class), intent(in):: MPID
      logical, intent(out):: ofram

      ! Local

      type(Redc_class), dimension(:,:), allocatable:: Redd_aux

      logical:: skip
      logical, dimension(:), allocatable:: dummy

      integer:: ia,itran,jtran,ni,np0,nti,if0,if1
      integer:: fitran,fjtran,ffitran,ffjtran,ffktran
      integer:: iJl,iJu,iJf,itermf,itermu,iterml,iti
      integer:: nr,nt,ntj,ntk
      integer, dimension(:), allocatable:: flag

      double precision:: RAM,DwN
      double precision, dimension(1):: nutout, dnl, vpr
      double precision, dimension(2):: nut
      double precision, dimension(3):: vphv, vplv
      double precision, dimension(6):: vphve, vplve
      double precision, dimension(:), allocatable:: vpp


      ! Routine name
      urou = 'omegabuildinI'

      ! Initialize out of RAM
      ofram = .False.

      ! Preliminar allocation of vpp, auxiliar vector to store
      ! frequencies
      np0 = 10000
      allocate(vpp(np0))
      allocate(flag(np0))

      ! For each atom
      do ia=1,nA

        ! Initialize as skipping this atom
        skip = .True.

        ! For every transition
        do jtran=1,Atom(ia)%ntran

          ! If it is PRD
          if (Atom(ia)%lemiss2(jtran)) then

            ! Cannot skip
            skip = .False.
            exit

          end if ! PRD line

        end do ! Transitions

        ! There are no PRD lines for this atom
        if (skip) cycle

        !
        ! Calculate nominal Doppler width
        !

        ! Maximum temperature
        if (Input%dws.eq.'MAX') then

          ! Calculate
          DwN = Atom(ia)%cDopp*sqrt(Input%maxT)

        ! Minimum temperature
        else if (Input%dws.eq.'MIN') then

          ! Calculate
          DwN = Atom(ia)%cDopp*sqrt(Input%minT)

        ! Fixed input Doppler width
        else if (Input%dws.eq.'NUM') then

          ! Calculate
          DwN = Input%dw*1d-9/c

        end if ! Type of nominal Doppler width

        ! For each output transition
        do jtran=1,Atom(ia)%ntran

          ! Skip CRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get terms
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml

          ! Master
          if (pid.eq.0) then

            ! Get frequency limits
            if0 = Atom(ia)%if0(jtran)
            if1 = Atom(ia)%if1(jtran)

          end if

          ! If MPI
          if (nproc.gt.1) then

            ! Share limits
            call MPI_BCAST(if0,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)
            call MPI_BCAST(if1,1,MPI_INTEGER,0,MPI_COMM_RT,ierr)

          end if ! MPI

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Get indexes
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJf = Atom(ia)%fst(jtran)%ilevell(fjtran)
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
            ffktran = Atom(ia)%itrano(ffjtran)

            ! Number of input transitions
            nti = Atom(ia)%tranoI(ffktran)%nt

            ! Count input transitions
            ntk = Atom(ia)%tranoI(ffktran)%nt

            ! Output transition frequency
            nutout(1) = Atom(ia)%FSfreq(iJu,itermu) - &
                        Atom(ia)%FSfreq(iJf,itermf)

      !
      ! Reset indent
      !

      ! Get output transition ranges
      call get_transition_out_limit(Atmo,Red,ia,ffjtran, &
                                    1,if0,if1, &
                                    Input%redi_pars(3), &
                                    DwN,nutout,Frec%omega)

      ! Allocate auxiliar structure to hold input frequencies
      allocate(Redd_aux(nti,Rz0:Rz1_PRD))

      ! Reset index of input transition
      iti = 0

      ! For each input transition
      do iti=1,Atom(ia)%tranoI(ffktran)%nt

        ! Get indexes
        ffitran = Atom(ia)%tranoI(ffktran)%indT(iti)
        itran = Atom(ia)%ifst(ffitran)
        fitran = Atom(ia)%ifstj(ffitran)

        ! Get indexes
        iterml = Atom(ia)%fst(itran)%iterml
        iJl = Atom(ia)%fst(itran)%ilevell(fitran)

        !
        ! Store resonance
        !

        ! Single resonance
        nr = 1
        dnl(1) = Atom(ia)%FSfreq(iJf,itermf) - &
                 Atom(ia)%FSfreq(iJl,iterml)

        !
        ! Store the frequencies of transitions
        !

        ! All transitions
        ntj = 1
        nt = 2
        nut(1) = Atom(ia)%FSfreq(iJu,itermu) - &
                 Atom(ia)%FSfreq(iJf,itermf)
        nut(2) = Atom(ia)%FSfreq(iJu,itermu) - &
                 Atom(ia)%FSfreq(iJl,iterml)

        ! Total number of limits
        ni = nr + nt

        ! Get input frequency and weights
        call get_input_frequencies(Atmo,Red,Redd_aux,Input%cohwi, &
                                   ia,ffjtran,iti,nr,ntj,nt,np0, &
                                   nut(2),nut(1), &
                                   Input%dcohwi,Input%redi_pars, &
                                   Frec%omega,dnl,nut, &
                                   Atom(ia)%cDopp, &
                                   vphv,vplv,vphve,vplve,vpr, &
                                   vpp,flag)

      end do ! Input transitions

      ! Split tasks weighted by number of frequencies
      call setmpi_red(Atom(ia),Red,Redd_aux,MPID,ia,ffjtran, &
                      Atom(ia)%tranoI(ffktran)%nt,dummy)


            !
            ! Restore identation
            !

          end do ! FS transition (out)
        end do ! Output transition

        ! Find frequency limits for the frequency range needed for
        ! integrals
        call find_integral_limits(Atom(ia),Atmo,Red, &
                                  Frec%omega,ia,.False.)

      end do ! Atoms

      ! Count memory
      call cram_red_frec(Red,RAM)
      FRAMc = RAM

      ! Check if everything is fine
      call control

      return

      end subroutine omegabuildinI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Allocate space to store the intensity redistribution
      !! functions\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!        ofram(logical): If reached the RAM limit
      subroutine allocate_WarrI(Atom,Red,Geom,ofram)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Red_class), intent(inout):: Red
      type(Geometry_class),intent(inout):: Geom
      logical, intent(out):: ofram

      ! Local

      integer:: nn,indx,jndx,iz0,iz1
      integer:: iz,ia,jtran,fjtran,ffjtran,ffktran,ffitran,iti

      double precision:: RAM,SRAM


      ! Initialize RAM counters
      WRAMc = 0d0

      ! Initialize
      ofram = .False.

      ! Generate scattering angles if angle-dependent
      if (.not.AVI) &
        call get_scattering(Geom)

      ! If not storing, go back
      if (.not.IRAM) return

      ! Get current RAM status
      RAM = cram_add(1)

      ! Allocate rzao
      allocate(Red%rzao(Red%nzao))

      ! For each index
      do indx=1,Red%nzao

        ! Count memory due to the bool
        RAM = RAM + 4d-6*dble(size(Red%zao(indx)%trani))
        nullify(Red%rzao(indx)%trani)

      end do ! Indexes

      ! For each atom
      do ia=1,nA

        ! For output transition
        do jtran=1,Atom(ia)%ntran

          ! Skip CRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Get indexes
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
            ffktran = Atom(ia)%itrano(ffjtran)
            jndx = Red%izao(ffjtran,ia,Rz0)

            ! Skip no work
            if (Red%ao(jndx)%nran.lt.1) cycle
            if (Red%ao(jndx)%nn(pid).lt.1) cycle

            ! Get z limits
            iz0 = Rz0 + (Red%ao(jndx)%Mi0(pid)-1)/Red%ao(jndx)%nfreq
            iz1 = Rz0 + (Red%ao(jndx)%Mi1(pid)-1)/Red%ao(jndx)%nfreq

            ! For each height
            do iz=iz0,iz1

              ! Get index
              indx = Red%izao(ffjtran,ia,iz)

              ! Allocate input transitions
              allocate(Red%rzao(indx)%trani(Atom(ia)% &
                                            tranoI(ffktran)%nt))

              ! For each input transition
              do iti=1,Atom(ia)%tranoI(ffktran)%nt

                ! Add RAM
                RAM = RAM + 1d-6*sizeof(Red%rzao(indx)%trani(iti))

              end do ! Input transitions
            end do ! Heights

            ! For each height
            do iz=iz0,iz1

              ! Get index
              indx = Red%izao(ffjtran,ia,iz)

              ! For each input transition
              do iti=1,Atom(ia)%tranoI(ffktran)%nt

                ! Initialize
                Red%rzao(indx)%trani(iti)%RAM = .False.

                ! Get indexes
                ffitran = Atom(ia)%tranoI(ffktran)%indT(iti)

                ! Predict size of next block
                nn = sum(Red%zao(indx)%trani(iti)%mfreq)

                ! Skip empty
                if (nn.le.0) cycle

                ! If angle-dependent
                if (.not.AVI) then

                  ! If Rayleigh scattering
                  if (ffjtran.eq.ffitran) then

                    ! All scattering angles but forward
                    nn = nn*(Geom%nScatt - 1)

                  ! Raman scattering
                  else

                    ! All scattering angles
                    nn = nn*Geom%nScatt

                  end if ! Rayleigh scattering
                end if ! Angle-dependent

                ! Get size in MB (single precision)
                SRAM = 4d-6*dble(nn)

                ! If no more space
                if (floor(RAM+SRAM).gt.RLIM.or.SRAM.le.0d0) then

                  ! Out of RAM and do not allocate
                  ofram = .True.
                  Red%rzao(indx)%trani(iti)%RAM = .False.
                  cycle

                end if ! Out of RAM

                ! Add to RAM count
                RAM = RAM + SRAM

                ! Saving
                Red%rzao(indx)%trani(iti)%RAM = .True.

                ! Allocate
                allocate(Red%rzao(indx)%trani(iti)%IWarr2(nn))

              end do ! Input transition
            end do ! Height
          end do ! Output transition (FS)
        end do ! Output transition
      end do ! Atoms

      ! Count RAM
      call cram_red_warr(Red,RAM)
      WRAMc = RAM

      endsubroutine allocate_WarrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resize some frequency dependent quantities and adjust indexes
      !! for each CPU taking into account the range of frequencies
      !! they need to take care of\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!     Input(Input_class): Structure with configuration data\n
      !!        MPID(MPI_class): Structure with MPI data
      subroutine frecresize(Frec,Atom,Input,MPID)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Frequency_class), intent(inout):: Frec
      type(Input_class), intent(inout):: Input
      type(MPI_class), intent(in):: MPID

      ! Local

      integer:: ia,itran,iproc,if0,if1


      ! If no MPI
      if (.not.MPID%mpi) then

        ! Ensure a minimum dimension for profile arrays later
        if (Frec%ntfreqi.lt.1) Frec%ntfreqi = 1
        if (Frec%npfreq.lt.1) Frec%npfreq = 1
        if (Frec%ntfreq.lt.1) Frec%ntfreq = 1

        ! And return
        return

      end if ! Serial


      ! Allocate size for profiles
      allocate(Frec%Mntfreq(0:nproc-1))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mntfreq)
      allocate(Frec%Mntfreqi(0:nproc-1))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mntfreqi)
      allocate(Frec%Mnpfreq(0:nproc-1))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mnpfreq)


      !
      ! If master
      !
      if (pid.eq.0) then

        ! Allocate and initialize limits for photoionizations
        allocate(Frec%Mlif0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(Frec%Mlif0)
        allocate(Frec%Mlif1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(Frec%Mlif1)
        allocate(Frec%Mpif0(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(Frec%Mpif0)
        allocate(Frec%Mpif1(0:nproc-1))
        MRAMc = MRAMc + 1d-6*sizeof(Frec%Mpif1)
        Frec%Mlif0(0) = 0
        Frec%Mlif1(0) = 0
        Frec%Mpif0(0) = 0
        Frec%Mpif1(0) = 0

        ! Reset Mntfreq, Mnpfreq, and Mnpfrei arays
        Frec%Mntfreq(0) = Frec%ntfreq
        Frec%Mntfreqi(0) = Frec%ntfreqi
        Frec%Mnpfreq(0) = Frec%npfreq
        Frec%Mntfreq(1:nproc-1) = 0
        Frec%Mntfreqi(1:nproc-1) = 0
        Frec%Mnpfreq(1:nproc-1) = 0

        ! For each CPU
        do iproc=1,nproc-1

          ! Get limits of that CPU
          if0 = MPID%if0(iproc)
          if1 = MPID%if1(iproc)

          ! Initialize master limits
          Frec%Mlif0(iproc) = 100000000
          Frec%Mlif1(iproc) = -1
          Frec%Mpif0(iproc) = 100000000
          Frec%Mpif1(iproc) = -1

          ! For each atom
          do ia=1,nA

            ! For each bound-bound transition
            do itran=1,Atom(ia)%ntran

              !
              ! Rearrange the limits taking into account the CPU
              ! limits
              !

              ! The line is totally absent in this CPU
              if (Atom(ia)%Mif0(itran,iproc).gt.if1.or. &
                  Atom(ia)%Mif1(itran,iproc).lt.if0) then

                ! Set no line data for this CPU
                Atom(ia)%Mif0(itran,iproc) = if1
                Atom(ia)%Mif1(itran,iproc) = if0-1
                Atom(ia)%MW0(itran,iproc) = 0d0
                Atom(ia)%MW1(itran,iproc) = 0d0
                Atom(ia)%fflag(itran)%Mabsent(iproc) = .True.

              ! The line is present
              else

                ! If the lower limit is out of range
                if (Atom(ia)%Mif0(itran,iproc).lt.if0) then

                  ! Adjust left limit
                  Atom(ia)%Mif0(itran,iproc) = if0
                  Atom(ia)%MW0(itran,iproc) = Frec%W_freq(if0)

                end if ! Lower limit out of CPU range

                ! If the upper limit is out of range
                if (Atom(ia)%Mif1(itran,iproc).gt.if1) then

                  ! Adjust right limit
                  Atom(ia)%Mif1(itran,iproc) = if1
                  Atom(ia)%MW1(itran,iproc) = Frec%W_freq(if1)

                end if ! Upper limit out of CPU range

                ! If there is only one frequency in this CPU,
                ! nullify second weight
                if (Atom(ia)%Mif0(itran,iproc).eq. &
                    Atom(ia)%Mif1(itran,iproc)) &
                  Atom(ia)%MW1(itran,iproc) = 0d0

                ! Add frequencies in this process to counters
                Frec%Mntfreq(iproc) = Frec%Mntfreq(iproc) + 1 + &
                                      Atom(ia)%Mif1(itran,iproc) - &
                                      Atom(ia)%Mif0(itran,iproc)
                Frec%Mntfreqi(iproc) = Frec%Mntfreqi(iproc) + (1 + &
                                       Atom(ia)%Mif1(itran,iproc) - &
                                       Atom(ia)%Mif0(itran,iproc))* &
                                       Atom(ia)%fst(itran)%nt

              end if ! Line presence

              ! Update line ranges
              if (Atom(ia)%Mif0(itran,iproc).lt.Frec%Mlif0(iproc)) &
                Frec%Mlif0(iproc) = Atom(ia)%Mif0(itran,iproc)
              if (Atom(ia)%Mif1(itran,iproc).gt.Frec%Mlif1(iproc)) &
                Frec%Mlif1(iproc) = Atom(ia)%Mif1(itran,iproc)

            end do ! bound-bound Transition

            ! For each bound-free transition
            do itran=1,Atom(ia)%nphot

              ! The line is totally absent in this CPU
              if (Atom(ia)%phot(itran)%Mif0(iproc).gt.if1.or. &
                  Atom(ia)%phot(itran)%Mif1(iproc).lt.if0) then

                ! Set no line data for this CPU
                Atom(ia)%phot(itran)%Mif0(iproc) = if1
                Atom(ia)%phot(itran)%Mif1(iproc) = if0-1
                Atom(ia)%phot(itran)%MW0(iproc) = 0d0
                Atom(ia)%phot(itran)%MW1(iproc) = 0d0
                Atom(ia)%phot(itran)%Mabsent(iproc) = .True.

              ! The line is present
              else

                ! If the lower limit is out of range
                if (Atom(ia)%phot(itran)%Mif0(iproc).lt.if0) then

                  ! Adjust left limit
                  Atom(ia)%phot(itran)%Mif0(iproc) = if0
                  Atom(ia)%phot(itran)%MW0(iproc) = Frec%W_freq(if0)

                end if ! Lower limit out of CPU range

                ! If the upper limit is out of range
                if (Atom(ia)%phot(itran)%Mif1(iproc).gt.if1) then

                  ! Adjust right limit
                  Atom(ia)%phot(itran)%Mif1(iproc) = if1
                  Atom(ia)%phot(itran)%MW1(iproc) = Frec%W_freq(if1)

                end if ! Upper limit out of CPU range

                ! If there is only one frequency in this CPU,
                ! nullify second weight
                if (Atom(ia)%phot(itran)%Mif0(iproc).eq. &
                    Atom(ia)%phot(itran)%Mif1(iproc)) &
                  Atom(ia)%phot(itran)%MW1(iproc) = 0d0

                ! Add frequencies in this process to counter
                Frec%Mnpfreq(iproc) = Frec%Mnpfreq(iproc) + 1 + &
                                  Atom(ia)%phot(itran)%Mif1(iproc) - &
                                  Atom(ia)%phot(itran)%Mif0(iproc)

              end if ! Line presence

              ! Update photoionization range
              if (Atom(ia)%phot(itran)%Mif0(iproc).lt. &
                  Frec%Mpif0(iproc)) &
                Frec%Mpif0(iproc) = Atom(ia)%phot(itran)%Mif0(iproc)
              if (Atom(ia)%phot(itran)%Mif1(iproc).gt. &
                  Frec%Mpif1(iproc)) &
                Frec%Mpif1(iproc) = Atom(ia)%phot(itran)%Mif1(iproc)

            end do ! bound-free Transition
          end do ! Atom

          ! Ensure minimum dimensionality in profiles
          if (Frec%Mntfreqi(iproc).lt.1) Frec%Mntfreqi(iproc) = 1
          if (Frec%Mnpfreq(iproc).lt.1) Frec%Mnpfreq(iproc) = 1
          if (Frec%Mntfreq(iproc).lt.1) Frec%Mntfreq(iproc) = 1

        end do ! CPU

      !
      ! If slave with jobs
      !
      else if (MPID%nf(pid).ge.1) then

        ! Store MPI scope limits in short variables
        if0 = MPID%if0(pid)
        if1 = MPID%if1(pid)

        ! If weights have more data than needed
        if (size(Frec%W_freq).ne.MPID%nf(pid)) then

          ! For each atom
          do ia=1,nA

            ! For each bound-bound transition
            do itran=1,Atom(ia)%ntran

              ! Skip if already absent
              if (Atom(ia)%fflag(itran)%absent) cycle

              ! If upper line limit above lower CPU limit
              if (Atom(ia)%if1(itran).lt.MPID%if0(pid)) then

                ! Flag absent
                Atom(ia)%fflag(itran)%absent = .True.
                cycle

              end if ! Line out of scope

              ! If lower line limit above upper CPU limit
              if (Atom(ia)%if0(itran).gt.MPID%if1(pid)) then

                ! Flag absent
                Atom(ia)%fflag(itran)%absent = .True.
                cycle

              end if ! Line out of scope

            end do ! bound-bound Transition
          end do ! Atom

        end if ! Need to resize

        ! Reset number of frequencies for profiles
        Frec%ntfreq = 0
        Frec%ntfreqi = 0
        Frec%npfreq = 0

        ! For each atom
        do ia=1,nA

          ! For each bound-bound transition
          do itran=1,Atom(ia)%ntran

            !
            ! Rearrange the limits taking into account the CPU limits
            !

            ! The line is totally absent
            if (Atom(ia)%fflag(itran)%absent) then

              ! Set empty data
              Atom(ia)%if0(itran) = if1
              Atom(ia)%if1(itran) = if0-1
              Atom(ia)%W0(itran) = 0d0
              Atom(ia)%W1(itran) = 0d0

            ! The line is present
            else

              ! If the lower limit is out of range
              if (Atom(ia)%if0(itran).lt.if0) then

                ! Adjust left limit
                Atom(ia)%if0(itran) = if0
                Atom(ia)%W0(itran) = Frec%W_freq(if0)

              end if ! Lower limit out of CPU range

              ! If the upper limit is out of range
              if (Atom(ia)%if1(itran).gt.if1) then

                ! Adjust right limit
                Atom(ia)%if1(itran) = if1
                Atom(ia)%W1(itran) = Frec%W_freq(if1)

              end if ! Upper limit out of CPU range

              ! If there is only one frequency in this CPU,
              ! nullify second weight
              if (Atom(ia)%if0(itran).eq.Atom(ia)%if1(itran)) &
                Atom(ia)%W1(itran) = 0d0

              ! Add frequencies to counters
              Frec%ntfreq = Frec%ntfreq + 1 + &
                            Atom(ia)%if1(itran) - Atom(ia)%if0(itran)
              Frec%ntfreqi = Frec%ntfreqi + (1 + &
                              Atom(ia)%if1(itran) - &
                              Atom(ia)%if0(itran))* &
                             Atom(ia)%fst(itran)%nt

            end if ! Line presence

          end do ! bound-bound Transition

          ! For each bound-free transition
          do itran=1,Atom(ia)%nphot

            ! The line is totally absent
            if (Atom(ia)%phot(itran)%if0.gt.if1.or. &
                Atom(ia)%phot(itran)%if1.lt.if0) then

              ! Set empty data
              Atom(ia)%phot(itran)%absent = .True.
              Atom(ia)%phot(itran)%if0 = if1
              Atom(ia)%phot(itran)%if1 = if0-1
              Atom(ia)%phot(itran)%W0 = 0d0
              Atom(ia)%phot(itran)%W1 = 0d0

            ! The line is present
            else

              ! If the lower limit is out of range
              if (Atom(ia)%phot(itran)%if0.lt.if0) then

                ! Adjust left limit
                Atom(ia)%phot(itran)%if0 = if0
                Atom(ia)%phot(itran)%W0 = Frec%W_freq(if0)

              end if ! Lower limit out of CPU range

              ! If the upper limit is out of range
              if (Atom(ia)%phot(itran)%if1.gt.if1) then

                ! Adjust right limit
                Atom(ia)%phot(itran)%if1 = if1
                Atom(ia)%phot(itran)%W1 = Frec%W_freq(if1)

              end if ! Upper limit out of CPU range

              ! If there is only one frequency in this CPU,
              ! nullify second weight
              if (Atom(ia)%phot(itran)%if0.eq. &
                  Atom(ia)%phot(itran)%if1) &
                Atom(ia)%phot(itran)%W1 = 0d0

              ! Add frequencies to counter
              Frec%npfreq = Frec%npfreq + 1 + &
                            Atom(ia)%phot(itran)%if1 - &
                            Atom(ia)%phot(itran)%if0

            end if ! Line presence

          end do ! bound-free Transition
        end do ! Atom

        ! Ensure minimum dimensionality
        if (Frec%ntfreqi.lt.1) Frec%ntfreqi = 1
        if (Frec%npfreq.lt.1) Frec%npfreq = 1
        if (Frec%ntfreq.lt.1) Frec%ntfreq = 1

        ! LTE lines
        do ia=1,nLTEl

          ! The line is totally absent
          if (Input%LTEline(ia)%absent) then

            ! Fake out of bound indexes
            Input%LTEline(ia)%if0 = if1
            Input%LTEline(ia)%if1 = if0-1

          ! The line is present
          else

            ! If the lower limit is out of range, adjust
            if (Input%LTEline(ia)%if0.lt.if0) &
              Input%LTEline(ia)%if0 = if0

            ! If the upper limit is out of range, adjust
            if (Input%LTEline(ia)%if1.gt.if1) &
              Input%LTEline(ia)%if1 = if1

          end if ! Line presence

        end do ! LTE lines

      !
      ! If slave without job
      !
      else

        ! Free weights
        deallocate(Frec%W_freq)
        allocate(Frec%W_freq(1))

        ! Reset number of frequencies for profiles
        Frec%ntfreq = 1
        Frec%ntfreqi = 1
        Frec%npfreq = 1

        ! For each atom
        do ia=1,nA

          ! For each bound-bound transition
          do itran=1,Atom(ia)%ntran

            ! Flag absent
            Atom(ia)%fflag(itran)%absent = .True.
            Atom(ia)%if0(itran) = 0
            Atom(ia)%if1(itran) = -1
            Atom(ia)%W0(itran) = 0d0
            Atom(ia)%W1(itran) = 0d0

          end do ! bound-bound Transition

          ! For each bound-free transition
          do itran=1,Atom(ia)%nphot

            ! Set empty data
            Atom(ia)%phot(itran)%absent = .True.
            Atom(ia)%phot(itran)%if0 = 0
            Atom(ia)%phot(itran)%if1 = -1
            Atom(ia)%phot(itran)%W0 = 0d0
            Atom(ia)%phot(itran)%W1 = 0d0

          end do ! bound-free Transition
        end do ! Atom

      end if ! Needs resize

      ! Check if everything is fine
      call control

      return

      end subroutine frecresize

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resize some frequency dependent quantities and adjust indexes
      !! for each CPU taking into account the range of frequencies
      !! they need to take care of, and create the output frequency
      !! axis, for a CLE synthesis\n
      !!     Input(Input_class): Structure with configuration data\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!        MPID(MPI_class): Structure with MPI data
      subroutine refitfrec(Input,Frec,Atom,MPID)

      ! I/O

      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Frequency_class), intent(inout):: Frec

      ! Local

      integer:: ia,iran,ifreq,jfreq,itran
      integer:: if0,if1,iif0,iif1,jf0,jf1,kf0,kf1
      integer, dimension(:), allocatable:: invmapping

      double precision, dimension(nfreq):: daux


      ! Store MPI limits in short variables
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)
      iif0 = MPID%iif0(pid)
      iif1 = MPID%iif1(pid)

      ! If size different to what the CPU holds
      if(size(Frec%W_freq).ne.MPID%inf(pid))then

        ! Store weights into a temporal array
        daux = Frec%W_freq

        ! Resize weights
        MRAMc = MRAMc - 1d-6*sizeof(Frec%W_freq)
        deallocate(Frec%W_freq)
        allocate(Frec%W_freq(iif0:iif1))
        MRAMc = MRAMc + 1d-6*sizeof(Frec%W_freq)

        ! Recover the data
        Frec%W_freq = daux(iif0:iif1)

      end if ! Need to resize

      !
      ! Create indexing axis
      !

      ! Set omega3
      allocate(Frec%omega3(nfreq))
      PRAMc = PRAMc + 1d-6*sizeof(Frec%omega3)
      Frec%omega3 = Frec%omega*Frec%omega*Frec%omega

      ! Allocate
      allocate(Frec%mapping(Input%lim_stk%nn))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%mapping)
      allocate(Frec%omega_ou(Input%lim_stk%nn))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%omega_ou)
      allocate(Frec%omega3_ou(Input%lim_stk%nn))
      PRAMc = PRAMc + 1d-6*sizeof(Frec%omega3_ou)

      ! Allocate "Master" line and photoionization limits
      allocate(Frec%Mlif0(pid:pid))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mlif0)
      allocate(Frec%Mlif1(pid:pid))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mlif1)
      allocate(Frec%Mpif0(pid:pid))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mpif0)
      allocate(Frec%Mpif1(pid:pid))
      MRAMc = MRAMc + 1d-6*sizeof(Frec%Mpif1)

      ! And initialize these limits
      Frec%Mlif0 = 10000000
      Frec%Mlif1 = -1
      Frec%Mpif0 = 10000000
      Frec%Mpif1 = -1
      Frec%lif0 = 10000000
      Frec%lif1 = -1
      Frec%pif0 = 10000000
      Frec%pif1 = -1

      ! For each atom
      do ia=1,nA

        ! Allocate input limits
        allocate(Atom(ia)%ilf0(Atom(ia)%ntran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%ilf0)
        allocate(Atom(ia)%ilf1(Atom(ia)%ntran))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%ilf1)
        allocate(Atom(ia)%ipf0(Atom(ia)%nphot))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%ipf0)
        allocate(Atom(ia)%ipf1(Atom(ia)%nphot))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%ipf1)

        ! For each bound-bound transition
        do itran=1,Atom(ia)%ntran

          ! Copy the limits
          Atom(ia)%ilf0(itran) = Atom(ia)%if0(itran)
          Atom(ia)%ilf1(itran) = Atom(ia)%if1(itran)

          !
          ! Rearrange the limits taking into account the CPU limits
          !

          ! The line is totally absent
          if (Atom(ia)%ilf1(itran).lt.iif0.or. &
              Atom(ia)%ilf0(itran).gt.iif1) then

            ! Set empty ranges
            Atom(ia)%ilf0(itran) = iif1
            Atom(ia)%ilf1(itran) = iif0-1
            Atom(ia)%W0(itran) = 0d0
            Atom(ia)%W1(itran) = 0d0

          ! The line is present
          else

            ! If the lower limit is out of range
            if (Atom(ia)%ilf0(itran).lt.iif0) then

              ! Adjust left limit
              Atom(ia)%ilf0(itran) = iif0
              Atom(ia)%W0(itran) = Frec%W_freq(iif0)

            end if ! Lower limit out of CPU range

            ! If the upper limit is out of range
            if (Atom(ia)%ilf1(itran).gt.iif1) then

              ! Adjust right limit
              Atom(ia)%ilf1(itran) = iif1
              Atom(ia)%W1(itran) = Frec%W_freq(iif1)

            end if ! Upper limit out of CPU range

            ! If there is only one frequency in this CPU, nullify
            ! second weight
            if (Atom(ia)%ilf0(itran).eq.Atom(ia)%ilf1(itran)) &
              Atom(ia)%W1(itran) = 0d0

          end if ! Line presence

          ! Update the line range
          if (Atom(ia)%ilf0(itran).lt.Frec%Mlif0(pid)) &
            Frec%Mlif0(pid) = Atom(ia)%ilf0(itran)
          if (Atom(ia)%ilf1(itran).gt.Frec%Mlif1(pid)) &
            Frec%Mlif1(pid) = Atom(ia)%ilf1(itran)

        end do ! bound-bound Transition

        ! For each bound-free transition
        do itran=1,Atom(ia)%nphot

          ! Copy
          Atom(ia)%ipf0(itran) = Atom(ia)%phot(itran)%if0
          Atom(ia)%ipf1(itran) = Atom(ia)%phot(itran)%if1

          ! The line is totally absent
          if (Atom(ia)%ipf0(itran).gt.iif1.or. &
              Atom(ia)%ipf1(itran).lt.iif0) then

            ! Set empty limits
            Atom(ia)%ipf0(itran) = iif1
            Atom(ia)%ipf1(itran) = iif0-1
            Atom(ia)%phot(itran)%W0 = 0d0
            Atom(ia)%phot(itran)%W1 = 0d0

          ! The line is present
          else

            ! If the lower limit is out of range
            if (Atom(ia)%ipf0(itran).lt.iif0) then

              ! Adjust left limit
              Atom(ia)%ipf0(itran) = iif0
              Atom(ia)%phot(itran)%W0 = Frec%W_freq(if0)

            end if ! Lower limit out of CPU range

            ! If the upper limit is out of range
            if (Atom(ia)%ipf1(itran).gt.iif1) then

              ! Adjust right limit
              Atom(ia)%ipf1(itran) = iif1
              Atom(ia)%phot(itran)%W1 = Frec%W_freq(if1)

            end if ! Upper limit out of CPU range

            ! If there is only one frequency in this CPU,
            ! nullify second weight
            if (Atom(ia)%ipf0(itran).eq.Atom(ia)%ipf1(itran)) &
              Atom(ia)%phot(itran)%W1 = 0d0

            ! Add number of frequencies to counter
            Frec%npfreq = Frec%npfreq + 1 + &
                          Atom(ia)%ipf1(itran) - &
                          Atom(ia)%ipf0(itran)

          end if ! Line presence

          ! Update the photoionization ranges
          if (Atom(ia)%ipf0(itran).lt.Frec%Mpif0(pid)) &
            Frec%Mpif0(pid) = Atom(ia)%ipf0(itran)
          if (Atom(ia)%ipf1(itran).gt.Frec%Mpif1(pid)) &
            Frec%Mpif1(pid) = Atom(ia)%ipf1(itran)

        end do ! bound-free Transition
      end do ! Atom

      ! If wrong line limits
      if (Frec%Mlif0(pid).gt.Frec%Mlif1(pid)) then

        ! Get fake limits
        Frec%Mlif0(pid) = iif1
        Frec%Mlif1(pid) = iif0-1

      end if ! Wrong line limits

      ! If wrong photoionization limits
      if (Frec%Mpif0(pid).gt.Frec%Mpif1(pid)) then

        ! Get fake limits
        Frec%Mpif0(pid) = iif1
        Frec%Mpif1(pid) = iif0-1

      end if ! Wrong photoionization limits

      ! If no limitations in output frequencies
      if (Input%lim_stk%nran.le.0) then

        ! For all frequencies
        do ifreq=1,nfreq

          ! Just copy the relevant data
          Frec%mapping(ifreq) = ifreq

        end do

        ! Copy vectors
        Frec%omega_ou = Frec%omega
        Frec%omega3_ou = Frec%omega3

        ! For each atom
        do ia=1,nA

          ! For each bound-bound transition
          do itran=1,Atom(ia)%ntran

            ! Copy input data and set presence
            Atom(ia)%if0(itran) = Atom(ia)%ilf0(itran)
            Atom(ia)%if1(itran) = Atom(ia)%ilf1(itran)
            Atom(ia)%fflag(itran)%absent = &
                                      Atom(ia)%if1(itran).lt.if0.or. &
                                      Atom(ia)%if0(itran).gt.if1

          end do ! bound-bound transitions

          ! For each bound-free transition
          do itran=1,Atom(ia)%nphot

            ! Copy input data and set presence
            Atom(ia)%phot(itran)%if0 = Atom(ia)%ipf0(itran)
            Atom(ia)%phot(itran)%if1 = Atom(ia)%ipf1(itran)
            Atom(ia)%phot(itran)%absent = &
                                 Atom(ia)%phot(itran)%if1.lt.if0.or. &
                                 Atom(ia)%phot(itran)%if0.gt.if1

          end do ! bound-free transitions
        end do ! Atoms

        ! And copy limits
        Frec%lif0 = Frec%Mlif0(pid)
        Frec%lif1 = Frec%Mlif1(pid)
        Frec%pif0 = Frec%Mpif0(pid)
        Frec%pif1 = Frec%Mpif1(pid)

      ! If limiting the output axis
      else

        ! Allocate inverse mapping
        allocate(invmapping(nfreq))
        invmapping = -1

        !
        ! Initialize everything to absent
        !

        ! For each atom
        do ia=1,nA

          ! For each bound-bound transition
          do itran=1,Atom(ia)%ntran

            ! Set absent
            Atom(ia)%fflag(itran)%absent = .True.

          end do ! Bound-bound transitions

          ! For each bound-free transition
          do itran=1,Atom(ia)%nphot

            ! Set absent
            Atom(ia)%phot(itran)%absent = .True.

          end do ! Bound-free transitions
        end do ! Atoms

        ! Initialize rolling index
        jfreq = 0

        ! For each output range
        do iran=1,Input%lim_stk%nran

          ! For each frequency in range
          do ifreq=Input%lim_stk%indx(1,iran), &
                   Input%lim_stk%indx(2,iran)

            ! Advance rolling index
            jfreq = jfreq + 1

            ! Copy correct frequency
            Frec%mapping(jfreq) = ifreq
            Frec%omega_ou(jfreq) = Frec%omega(ifreq)
            Frec%omega3_ou(jfreq) = Frec%omega3(ifreq)

            ! If below CPU limit, skip
            if (jfreq.lt.MPID%if0(pid)) cycle
            ! If above CPU limit, skip
            if (jfreq.gt.MPID%if1(pid)) cycle

            ! Inverse mapping if in this CPU
            invmapping(ifreq) = jfreq

          end do ! Frequency in range
        end do ! Ranges

        ! For each range (again)
        do iran=1,Input%lim_stk%nran

          ! Translate range limits
          iif0 = invmapping(Input%lim_stk%indx(1,iran))
          iif1 = invmapping(Input%lim_stk%indx(2,iran))

          ! If this range is out of mine, skip rest
          if (iif0.gt.MPID%if1(pid)) cycle
          if (iif1.lt.MPID%if0(pid)) cycle

          ! Get limits
          if0 = max(iif0,MPID%if0(pid))
          if1 = min(iif1,MPID%if1(pid))

          !
          ! Check lines present
          !

          ! For each atom
          do ia=1,nA

            ! For each bound-bound transition transition
            do itran=1,Atom(ia)%ntran

              ! If present, skip
              if (.not.Atom(ia)%fflag(itran)%absent) cycle

              ! Get limits
              jf0 = Atom(ia)%if0(itran)
              jf1 = Atom(ia)%if1(itran)

              ! Translate
              kf0 = invmapping(jf0)
              kf1 = invmapping(jf1)

              ! Find lower
              do while (kf0.lt.1)
                jf0 = jf0 + 1
                if (jf0.gt.nfreq) then
                  jf0 = -1
                  kf0 = -1
                  exit
                end if
                kf0 = invmapping(jf0)
              end do

              ! If bad, skip
              if (jf0.lt.1) cycle

              ! Find upper
              do while (kf1.lt.1)
                jf1 = jf1 - 1
                if (jf1.lt.1) then
                  jf1 = -1
                  kf1 = -1
                  exit
                end if
                kf1 = invmapping(jf1)
              end do

              ! If bad, skip
              if (jf1.lt.if0) cycle
              if (kf1.lt.if0) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.lt.jf0) cycle
              if (kf1.lt.kf0) cycle

              ! Save limit
              Atom(ia)%if0(itran) = kf0
              Atom(ia)%if1(itran) = kf1

              ! Is present here!
              Atom(ia)%fflag(itran)%absent = .False.

              ! Update totals
              if (Atom(ia)%if0(itran).lt.Frec%lif0) &
                Frec%lif0 = Atom(ia)%if0(itran)
              if (Atom(ia)%if1(itran).gt.Frec%lif1) &
                Frec%lif1 = Atom(ia)%if1(itran)

            end do ! Bound-bound transitions

            ! b-f transition
            do itran=1,Atom(ia)%nphot

              ! If present, skip
              if (.not.Atom(ia)%phot(itran)%absent) cycle

              ! Get limits
              jf0 = Atom(ia)%phot(itran)%if0
              jf1 = Atom(ia)%phot(itran)%if1

              ! Translate
              kf0 = invmapping(jf0)
              kf1 = invmapping(jf1)

              ! Find lower
              do while (kf0.lt.1)
                jf0 = jf0 + 1
                kf0 = invmapping(jf0)
                if (jf0.gt.nfreq) then
                  jf0 = -1
                  kf0 = -1
                  exit
                end if
              end do

              ! Bad limit, skip
              if (jf0.lt.1) cycle

              ! Find upper
              do while (kf1.lt.1)
                jf1 = jf1 - 1
                kf1 = invmapping(jf1)
                if (jf1.lt.1) then
                  jf1 = -1
                  kf1 = -1
                  exit
                end if
              end do

              ! If bad ranges, skip
              if (jf1.lt.if0) cycle
              if (kf1.lt.if0) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.gt.nfreq) cycle
              if (jf1.lt.jf0) cycle
              if (kf1.lt.kf0) cycle

              ! Translate limits
              Atom(ia)%phot(itran)%if0 = kf0
              Atom(ia)%phot(itran)%if1 = kf1

              ! Is present here!
              Atom(ia)%phot(itran)%absent = .False.

              ! Update totals
              if (Atom(ia)%phot(itran)%if0.lt.Frec%pif0) &
                Frec%pif0 = Atom(ia)%phot(itran)%if0
              if (Atom(ia)%phot(itran)%if1.gt.Frec%pif1) &
                Frec%pif1 = Atom(ia)%phot(itran)%if1

            end do ! Bound-free transitions
          end do ! Atoms
        end do ! Ranges

        ! Free
        deallocate(invmapping)

      end if ! Limited output

      ! If valid photoionization ranges
      if (Frec%pif1.ge.Frec%pif0) then

        ! Allocate photoionization exponential
        allocate(Frec%exu(Frec%pif0:Frec%pif1,1))
        PRAMc = PRAMc + 1d-6*sizeof(Frec%exu)

      end if ! Valid photoionization range

      ! Check if everything is fine
      call control

      return

      end subroutine refitfrec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Index and allocate arrays for the normalization data and
      !! estimate the minimum RAM neccesary\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!                pol(logical): If this is for polarization
      subroutine index_norm(Atom,LTElines,Atmo,Bfield,Geom,Red,pol)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Geometry_class), intent(in):: Geom
      type(Red_class), intent(inout):: Red
      logical, intent(in):: pol

      ! Local

      logical:: lvel,field

      integer:: njdir,ntran,idir,iz,ia,nl,ncom
      integer:: jtran,ktran,fjtran,ffjtran,ffktran

      double precision:: vel


      !
      ! Number of directions 
      !

      ! If dynamic
      if (dyn) then

        ! Get quadrature number of directions
        njdir = Geom%njdir

      ! Static
      else

        ! One direction is enough
        njdir = 1

      end if ! Dynamic

      !
      ! Number of atomic transitions
      !

      ! If polarization problem
      if (pol) then

        ! Term-term transitions
        ntran = nxtran

      ! If intensity problem
      else

        ! Level-level transitions
        ntran = nxt

      end if ! Polarization/intensity

      ! If LTE lines
      if (allocated(LTElines)) then

        ! Get size
        nl = size(LTElines)

      ! No LTE lines
      else

        ! Zero size
        nl = 0

      end if ! Allocated LTE lines

      ! Allocate indexing space
      allocate(Red%idzao(ntran+nl,Rz0:Rz1,njdir))

      ! Initialize counter
      Red%ndzao = 0

      ! Reset RAM estimation for normalization
      DRAMc = 0d0

      ! For each direction
      do idir=1,njdir

        ! For each considered height
        do iz=Rz0,Rz1

          ! If dynamic
          if (dyn) then

            ! Check local velocity
            vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))
            lvel = vel.gt.TINYVEL

          end if ! Dynamic

          ! Magnetic?
          field = Bfield%Bstrength(iz).gt.TINYB

          ! For each atom
          do ia=1,nA

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Skip absent
              if (Atom(ia)%fflag(jtran)%absent) cycle

              !
              ! If polarized
              !
              if (pol) then

                ! Rolling index
                ktran = jtran + Atom(ia)%tshift

                ! If valid in terms of velocity
                if (lvel.or.idir.eq.1) then

                  ! Advance counter
                  Red%ndzao = Red%ndzao + 1

                  ! Add RAM stimation (logical)
                  DRAMc = DRAMc + 4d-6

                  ! If magnetic
                  if (field) then

                    ! Get ncom
                    ncom = Atom(ia)%trano(jtran)%ncomB
                    if (ncom.lt.1) ncom = 1

                  ! Non magnetic
                  else

                    ! Get ncom
                    ncom = Atom(ia)%trano(jtran)%ncomNB
                    if (ncom.lt.1) ncom = 1

                  end if ! Magnetic

                  ! Add RAM for norms
                  DRAMc = DRAMc + 8d-6*dble(ncom)

                  ! Store index
                  Red%idzao(ktran,iz,idir) = Red%ndzao

                ! Not valid in terms of velocity
                else

                  ! Store index
                  Red%idzao(ktran,iz,idir) = Red%idzao(ktran,iz,1)

                end if ! Valid in terms of velocity

              !
              ! Not polarized
              !
              else

                ! For each FS transition
                do fjtran=1,Atom(ia)%fst(jtran)%nt

                  ! Get rolling indexes
                  ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
                  ffktran = ffjtran + Atom(ia)%tfshift

                  ! If valid in terms of velocity
                  if (lvel.or.idir.eq.1) then

                    ! Advance counter
                    Red%ndzao = Red%ndzao + 1

                    ! Add RAM estimation (logical + double)
                    DRAMc = DRAMc + 12d-6

                    ! Store index
                    Red%idzao(ffktran,iz,idir) = Red%ndzao

                  ! Not valid in terms of velocity
                  else

                    ! Store index
                    Red%idzao(ffktran,iz,idir) = &
                                               Red%idzao(ffktran,iz,1)

                  end if ! Valid in termsof velocity

                end do ! FS transitions

              end if ! Polarized

            end do ! Transitions
          end do ! Atoms
        end do ! Heights
      end do ! Directions

      ! Save size for only active atoms (no LTE lines)
      Red%ndzaoA = Red%ndzao

      ! For each direction
      do idir=1,njdir

        ! For each considered height
        do iz=Rz0,Rz1

          ! If dynamic
          if (dyn) then

            ! Check local velocity
            vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))
            lvel = vel.gt.TINYVEL

          end if ! Dynamic

          ! Magnetic?
          field = Bfield%Bstrength(iz).gt.TINYB

          ! For each LTE line
          do ia=1,nl

            ! Skip absent
            if (LTElines(ia)%absent) cycle

            ! Skip too high z
            if (iz.lt.LTElines(ia)%Rz0) cycle

            ! If valid in terms of velocity
            if (lvel.or.idir.eq.1) then

              ! Advance counter
              Red%ndzao = Red%ndzao + 1

              ! Add RAM stimation (logical)
              DRAMc = DRAMc + 4d-6

              ! Store index
              Red%idzao(ntran+ia,iz,idir) = Red%ndzao

            ! Not valid in terms of velocity
            else

              ! Store index
              Red%idzao(ntran+ia,iz,idir) = Red%idzao(ntran+ia,iz,1)

            end if ! Valid in terms of velocity

          end do ! LTE lines
        end do ! Heights
      end do ! Directions

      ! Allocate Red array
      allocate(Red%dzao(Red%ndzao))

      end subroutine index_norm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Index array for the redistribution quantities and allocate
      !! the structures to hold the input frequency axis\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!        Red(Red_class): Structure with redistribution input
      !!                        frequency data, redistribution
      !!                        function data, and profile or
      !!                        normalization data\n
      !!  Bstrength(double(:)): Magnetic field strength\n
      !!          pol(logical): If this is for polarization
      subroutine index_red(Atom,Red,Bstrength,pol)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Red_class), intent(inout):: Red
      logical, intent(in):: pol
      double precision, dimension(:), intent(in):: Bstrength

      ! Local

      integer:: ia,jtran,it,ffjtran,ffktran,indx,itran,iz
      integer:: mina,maxa,minto,maxto,nat,nti,ncom


      !
      ! Count maximum index of PRD atom and transition

      ! Initialize atomic and transition indexes, and counter
      ! of real elements
      mina = 10000
      maxa = 0
      minto = 10000
      maxto = 0
      nat = 0

      !
      ! If polarized
      !
      if (pol) then

        ! For each atom
        do ia=1,nA

          ! For all transitions
          do jtran=1,Atom(ia)%ntran

            ! Skip CRD
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! Add to counter
            nat = nat + 1

            ! Update limits
            if (jtran.lt.minto) minto = jtran
            if (jtran.gt.maxto) maxto = jtran
            if (ia.lt.mina) mina = ia
            if (ia.gt.maxa) maxa = ia

          end do ! Transitions
        end do ! Atoms

      !
      ! If intensity
      !
      else

        ! For each atom
        do ia=1,nA

          ! For all transitions
          do jtran=1,Atom(ia)%ntran

            ! Skip CRD
            if (.not.Atom(ia)%lemiss2(jtran)) cycle

            ! Go by fine structure
            do it=1,Atom(ia)%fst(jtran)%nt

              ! Get index
              ffjtran = Atom(ia)%ifst_ij(it,jtran)

              ! Add to counter
              nat = nat + 1

              ! Update limits
              if (ffjtran.lt.minto) minto = ffjtran
              if (ffjtran.gt.maxto) maxto = ffjtran
              if (ia.lt.mina) mina = ia
              if (ia.gt.maxa) maxa = ia

            end do ! Fine structure
          end do ! Transitions
        end do ! Atoms

      end if ! Polarization or intensity

      ! Sanity check
      if (mina.gt.maxa.or.minto.gt.maxto) then
        Red%nao = 0
        Red%nzao = 0
        return
      end if

      ! Allocate indexing array and first step of Frec and Red
      allocate(Red%izao(minto:maxto,mina:maxa,Rz0:Rz1_PRD))
      Red%nao = nat
      Red%nzao = nat*(Rz1_PRD-Rz0+1)

      ! Allocate frequency data
      allocate(Red%ao(Red%nao))
      allocate(Red%zao(Red%nzao))

      ! For each index
      do indx=1,Red%nzao

        ! Nullify input transition data
        nullify(Red%zao(indx)%trani)

      end do ! Indexes


      !
      ! Build index
      !

      !
      ! If polarized
      !
      if (pol) then

        !
        ! Build index

        ! Initialize
        nat = 0

        ! For each height
        do iz=Rz0,Rz1_PRD

          ! For each atom
          do ia=mina,maxa

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Skip CRD
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! Advance and store index
              nat = nat + 1
              Red%izao(jtran,ia,iz) = nat

              ! Initialize number of input frequencies
              Red%zao(nat)%mxfreq = 0

              ! Allocate
              if (iz.eq.Rz0) then
                allocate(Red%ao(nat)%Mi0(0:nproc-1))
                allocate(Red%ao(nat)%Mi1(0:nproc-1))
                allocate(Red%ao(nat)%nn(0:nproc-1))
              end if

              !
              ! Count minimum memory for normalization
              !

              ! Magnetic
              if (Bstrength(iz).gt.TINYB) then

                ! Get ncom
                ncom = Atom(ia)%trano(jtran)%ncomB
                if (ncom.lt.1) ncom = 1

              ! Non-magnetic
              else

                ! Get ncom
                ncom = Atom(ia)%trano(jtran)%ncomNB
                if (ncom.lt.1) ncom = 1

              end if ! Magnetic field

              ! Add logical and components norm
              DRAM2c = DRAM2c + 4d-6 + 8d-6*dble(ncom)

              !
              ! Allocate input transition structures
              !

              ! Number of input transitions
              nti = Atom(ia)%trano(jtran)%nt

              ! Allocate
              allocate(Red%zao(nat)%trani(nti))

              ! For each input transition
              do itran=1,nti

                ! Initialize sizes
                Red%zao(nat)%trani(itran)%osize = 0
                Red%zao(nat)%trani(itran)%isize = 0

              end do ! Input transitions
            end do ! Transitions
          end do ! Atoms
        end do ! Heights

      !
      ! If intensity
      !
      else

        ! Initialize
        nat = 0

        ! For each height
        do iz=Rz0,Rz1_PRD

          ! For each atom
          do ia=1,nA

            ! For each transition
            do jtran=1,Atom(ia)%ntran

              ! Skip CRD
              if (.not.Atom(ia)%lemiss2(jtran)) cycle

              ! For each fine structure transition
              do it=1,Atom(ia)%fst(jtran)%nt

                ! Get index
                ffjtran = Atom(ia)%ifst_ij(it,jtran)

                ! Advance and store
                nat = nat + 1
                Red%izao(ffjtran,ia,iz) = nat

                ! Initialize input frequency size
                Red%zao(nat)%mxfreq = 0

                ! Allocate
                if (iz.eq.Rz0) then
                  allocate(Red%ao(nat)%Mi0(0:nproc-1))
                  allocate(Red%ao(nat)%Mi1(0:nproc-1))
                  allocate(Red%ao(nat)%nn(0:nproc-1))
                end if

                ! Add logical and norm necessary sizes
                DRAM2c = DRAM2c + 12d-6

                ! Rolling index
                ffktran = Atom(ia)%itrano(ffjtran)

                ! Number of input transitions
                nti = Atom(ia)%tranoI(ffktran)%nt

                ! Allocate
                allocate(Red%zao(nat)%trani(nti))

                ! For each input transition
                do itran=1,nti

                  ! Initialize
                  Red%zao(nat)%trani(itran)%osize = 0
                  Red%zao(nat)%trani(itran)%isize = 0

                end do ! Input transitions
              end do ! Fine structure transitions
            end do ! Transitions
          end do ! Atoms
        end do ! Heights

      end if ! Polarization or intensity

      end subroutine index_red

!#####################################################################
!#####################################################################
!#####################################################################

      !> Checks where the non-coherent lower term approximation can be
      !! applied\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!  JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                          over the absorption profile\n
      !!    Bfield(Bfield_class): Structure with magnetic field data
      subroutine check_nchlt(Atom,JKQ,Bfield)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Bfield_class), intent(in):: Bfield
      complex(kind=8), dimension(-2:2,0:2,nxtran,Rz0:Rz1), &
                       intent(in):: JKQ

      ! Local

      double precision, parameter:: Bsat = 10d0

      logical:: skip

      integer:: ia,iz,jtran,itermu,itermf,iterml,itran,itran0,iJl
      integer:: itranmin,itranmax,litran,iti
      integer, dimension(:), allocatable:: nmsg

      double precision:: Blu,gJ,rJ,S,rL,Bcrit,Bcrit0,efield


      ! Check that there is at least one height with non-zero
      ! magnetic field

      ! Initialize flag
      skip = .True.

      ! For each considered height
      do iz=Rz0,Rz1

        ! If there is magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Cannot skip
          skip = .False.
          exit

        end if ! There is magnetic field

      end do ! Heights

      ! If no field, not necessary to check
      if (skip) return

      ! For each Atom
      do ia=1,nA

        ! Initialize to skip
        skip = .True.

        ! For each transition
        do jtran=1,Atom(ia)%ntran

          ! If PRD line
          if (Atom(ia)%lemiss2(jtran)) then

            ! We cannot skip this atom
            skip = .False.
            exit

          end if ! PRD line

        end do ! Transitions

        ! There are no PRD lines for this atom
        if (skip) cycle


        !
        ! Allocate space for NCHLT approximation
        !

        ! Find limit in transition indexes
        itranmin = Atom(ia)%ntran + 1
        itranmax = -1

        ! For every transition
        do jtran=1,Atom(ia)%ntran

          ! Skip CRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get terms
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml

          ! Update limits
          itranmin = min(itranmin,minval(Atom(ia)%trano(jtran)%indT))
          itranmax = max(itranmin,maxval(Atom(ia)%trano(jtran)%indT))

        end do ! Output transitions

        ! Allocate
        allocate(Atom(ia)%NCHLT(Rz0:Rz1,itranmin:itranmax))
        MRAMc = MRAMc + 1d-6*sizeof(Atom(ia)%NCHLT)

        ! Initialize transition shift
        itran0 = Atom(ia)%tshift

        ! Allocate logical array to not repeat messages
        allocate(nmsg(Atom(ia)%nMulti))
        nmsg = 0

        ! For every transition
        do jtran=1,Atom(ia)%ntran

          ! Skip CRD
          if (.not.Atom(ia)%lemiss2(jtran)) cycle

          ! Get terms
          itermu = Atom(ia)%fst(jtran)%itermu
          itermf = Atom(ia)%fst(jtran)%iterml

          ! For each input transition
          do iti=1,Atom(ia)%trano(jtran)%nt

            ! Get indexes
            litran = Atom(ia)%trano(jtran)%indT(iti)
            iterml = Atom(ia)%fst(litran)%iterml

            ! Get Blu
            Blu = Atom(ia)%Ecoeff(iterml,itermu)

            ! Get L and S
            rL = Atom(ia)%rLval(iterml)
            S = Atom(ia)%Sval(iterml)

            ! Get real itran
            itran = litran + itran0

            ! For each height
            do iz=Rz0,Rz1

              ! Critical field (Blu has factor 10^8)
              Bcrit0 = 1.137d1*Blu*dble(JKQ(0,0,itran,iz))*Bsat

              ! Initialize applicability
              skip = .True.

              ! Effective magnetic field
              efield = Bfield%Bstrength(iz)*sin(Bfield%Btheta(iz))

              ! For level in the lower term
              do iJl=1,Atom(ia)%nJ(iterml)

                ! Get J
                rJ = Atom(ia)%rJval(iJl,iterml)

                !
                ! Get Lande factor

                ! Multi-level
                if (Atom(ia)%ML) then

                  ! From atom
                  gJ = Atom(ia)%gL(iterml)

                ! Multi-term
                else

                  ! Assume LS coupling
                  gJ = 1d0 + .5d0*(rJ*(rJ+1d0) + S*(S+1d0) - &
                                   rL*(rL+1d0))/rJ/(rJ+1d0)
                end if

                ! If non-zero Landé
                if (gJ.gt.0d0) then

                  ! Reduce critical limit
                  Bcrit = Bcrit0/gJ

                ! No Landé
                else

                  ! Indinite critical
                  Bcrit = 1d99

                end if ! Landé factor

                ! Check field is big enough
                if (efield.lt.Bcrit) then

                  ! Cannot skip with non-coherent
                  skip = .False.

                  ! If magnetic field, add message
                  if (Bfield%Bstrength(iz).gt.TINYB) &
                    nmsg(iterml) = nmsg(iterml) + 1

                  ! And leave
                  exit

                ! Field is big enough
                end if

              end do ! Every lower level

              ! Store if non-coherent
              Atom(ia)%NCHLT(iz,litran) = skip

            end do ! Every height
          end do ! Input transitions
        end do ! Output transitions

        !
        ! Now print the messages
        !

        ! For each lower term
        do iterml=1,Atom(ia)%nMulti

          ! If there are messages for this term
          if (nmsg(iterml).gt.0) then

            ! Multi-level atom
            if (Atom(ia)%ML) then

              ! Write message
              write(umsg,'(A,A,1x,i3,1x,A,1x,i5,1x,A)') &
                ' # WARNING: The magnetic field is smaller than ', &
                'the saturation field of the level',iterml,'for ', &
                nmsg(iterml),'heights'

            ! Multi-term atom
            else

              ! Write message
              write(umsg,'(A,A,1x,i3,1x,A,1x,i5,1x,A)') &
                ' # WARNING: The magnetic field is smaller than ', &
                'the saturation field of the term',iterml,'for ', &
                nmsg(iterml),'heights'

            end if ! Multi level or multi term

            ! Verbose
            call verbose

          end if ! Messages to print

        end do ! Lower terms

        ! Deallocate logical array of this atom
        deallocate(nmsg)

      end do ! Atoms

      return

      end subroutine check_nchlt

!#####################################################################
!#####################################################################
!#####################################################################

      end module omegabuild_mod
