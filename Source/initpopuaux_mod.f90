      !> Initializations for populations
      module initpopuaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
!  Last version:
!     18/12/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/12/2025:    V4.0.1 - Bugfix: When loading populations form a
!                             file in 1.5DS mode, missing reading of
!                             the nx and ny sizes (TdPA)
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
!  rPopu
!    Read atomic populations from a specified file
!
!  LTEPopu
!    Compute LTE populations for a given model atom and model
!  atmosphere
!
!  LTEPopu_line
!    Compute LTE populations for a given LTE line in a given model
!  atmosphere
!
!  InitHpopu
!    Set the hydrogen populations in the model atom to the ones in the
!  model atmosphere atmosphere, distributed in LTE among sublevels of
!  the given principal levels
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use math_mod
      use parameters_mod , only : hplanck , PI , me, kb , qel , &
                                  pi4eps0 , fktoJ , TINYO
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read atomic populations from a specified file\n
      !!  filename(character(:)): Name of the file to read\n
      !!        Atom(Atom_class): Structure with atomic data
      subroutine rPopu(filename,Atom)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: filename

      ! Local

      character(len=2):: fileid

      integer:: iz, ios, nZ_popu, NL_popu, offset

      double precision, dimension(Atom%nlevel):: buff


      ! Routine name
      urou = 'rPopu'

      !
      ! If inversion
      !
      if (run_mode.eq.-1) then

        ! It is not possible to read population files, issue error
        umsg = 'Reading of populations not implemented in inversion'
        call aborted
        return

      end if ! Type of run

      ! Open the file with the populations
      open (200,file=trim(filename), &
            status='unknown', iostat=ios, err=1000, access='stream', &
            action='read', form='unformatted')

      ! Get file ID
      read (200,err=1100) fileid

      ! Only 1.5D
      if (run_mode.eq.1) then

        ! Check that it has the correct population file ID for 1.5D
        if (fileid.ne.'2D') then

          ! Error due to wrong ID
          umsg = 'The ID of the population file for '// &
                 'Atom'//Atom%Element//' is not valid.'
          call aborted
          close(200)
          return

        end if ! Wrong ID

        ! Read the other half of the ID
        read (200,err=1100) fileid

        ! Read nx and ny
        read (200,err=1100) ios
        read (200,err=1100) ios

      end if ! Only 1.5D

      ! Check that it has the correct population file ID
      if (fileid.ne.'bp'.and.fileid.ne.'bd') then

        ! Error due to wrong ID
        umsg = 'The ID of the population file for '// &
               'Atom'//Atom%Element//' is not valid.'
        call aborted
        return

      end if ! Check file ID

      ! Read dimensionality of input
      read (200,err=1100) nZ_popu

      ! Check if it has the correct number of nodes
      if (nZ.ne.nZ_popu) then

        ! Error due to incompatible model sizes
        umsg = 'Node points in atmosphere and '// &
               'populations do not coincide.'
        call aborted
        close(200)
        return

      end if ! Check number of nodes

      ! Read dimensionality of model atom
      read (200,err=1100) nL_popu

      ! Check if it has the correct number of levels
      if(Atom%nlevel.ne.nL_popu)then

        ! Error due to incompatible model sizes
        write(umsg,'(A,1x,i4,1x,A,1x,i4)') &
              'Number of levels in atom and populations'// &
              ' does not coincide:',Atom%nlevel,'!=',nL_popu
        call aborted
        close(200)
        return

      end if ! Check number of levels

      ! Only 1.5D mode
      if (run_mode.eq.1) then

        ! Jump to current column
        offset = nz*8*nL_popu*(icoords(3)-1)
        call fseek(200,offset,1)

      end if ! 1.5D mode

      ! If proper populations in number density
      if (fileid.eq.'bp') then

        ! Allcoate space for populations
        allocate(Atom%popu(Atom%nlevel,nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%popu)

        ! For each height
        do iz=1,nZ

          ! Read datum
          read (200,err=1100) buff

          ! Store read data into population array
          Atom%popu(:,iz) = buff

          ! Calculate the total population from this input
          Atom%n(iz) = sum(Atom%popu(:,iz))

        end do ! Heights

      ! If departure coefficients
      else if (fileid.eq.'bd') then

        ! Allocate space for departure coefficients
        allocate(Atom%depar(Atom%nlevel,nz))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%depar)

        ! For each height
        do iz=1,nZ

        ! Read datum
          read (200,err=1100) buff

          ! Store read data into population array
          Atom%depar(:,iz) = buff

        end do ! Heights

      end if ! Population or departure coefficient

      ! Close file
      close (200)

      ! Check if everything is fine
      call control
      if (laborted) return

      ! If global master
      if (gpid.eq.0) then

        ! If read populations
        if (fileid.eq.'bp') then

          ! Prepare message
          umsg = ' - Populations '//trim(filename)//' read'

        ! If departure coefficients
        else if (fileid.eq.'bd') then

          ! Prepare message
          umsg = ' - Departure coefficients '//trim(filename)//' read'

        end if ! Population or departure coefficient

        ! Write
        call verbose

      end if ! Global master

      return

1000  umsg = 'Error opening population file'
      call aborted
      return
1100  umsg = 'Error reading population file'
      close (200)
      call aborted
      return

      end subroutine rPopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute LTE populations for a given model atom and model
      !! atmosphere\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data
      subroutine LTEPopu(Atom,Atmo)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: ilevel,iterm,iJ,ii,iim,iz,dZ

      double precision:: C0,C1,C2,zm,dby,Tin,S,dE,gi0,arg
      double precision, dimension(:), allocatable:: debey


      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      ! 1d3 sqrt(cm^-3 -> m^-3) electron density below
      C1 = sqrt(8d0*PI/kb)*((qel*qel/pi4eps0)**(1.5d0))*1d3


      !
      ! Calculate Debey correction to the ionization potential
      !

      ! Allocate and initialize the correction
      allocate(debey(Atom%nlevel))
      debey = 0d0

      ! Initialize the index
      ilevel = 1

      ! For all terms
      do iterm=1,Atom%nMulti

        ! For all level within the term
        do iJ=1,Atom%nJ(iterm)

          ! Ignore the ground level
          if (iterm.eq.1.and.iJ.eq.1) cycle

          ! Determine the level index, the charge of the
          ! shell nucleus + rest of electrons and the
          ! change of stages between this level and the
          ! ground level of the model
          ilevel = ilevel + 1
          zm = Atom%stage(iterm) - 1
          iim = Atom%stage(iterm) - Atom%stage(1)

          ! Add the contribution to the Debey correction
          do ii=1,iim

            ! Contribution for this level
            debey(ilevel) = debey(ilevel) + zm

            ! Advance charge
            zm = zm + 1

          end do ! Debey contributions
        end do ! J levels
      end do ! Term


      !
      ! Compute LTE populations
      !

      ! For each height
      do iz=1,nz

        ! Inverse of temperature
        Tin = 1d0/Atmo%T(iz)
        ! Multiplicative factor for Debey correction
        dby = C1*sqrt(Atmo%ne(iz)*Tin)
        ! Non atomic part of the Saha function
        C2 = .5d6*Atmo%ne(iz)*((C0*Tin)**(1.5d0))
        ! Cumulative factor reset
        S = 1d0

        ! Initialize index
        ilevel = 1

        ! For each term
        do iterm=1,Atom%nMulti

          ! For each level within the term
          do iJ=1,Atom%nJ(iterm)

            ! Skip ground level
            if (iterm.eq.1.and.iJ.eq.1) cycle

            ! Determine the index and the ionization potential
            ilevel = ilevel + 1
            dE = (Atom%FSfreq(iJ,iterm) - Atom%FSfreq(1,1))*fktoJ

            ! If it is the hard-coded hydrogen
            if (Atom%cust) then

              ! Degeneration from variable
              gi0 = Atom%deg(ilevel)/Atom%deg(1)

            ! If it is a normal atom
            else

              ! Degeneracy from total angular momentum
              gi0 = (2d0*Atom%rJval(iJ,iterm) + 1d0)/ &
                    (2d0*Atom%rJval(1,1) + 1d0)

            end if ! If hard-coded hydrogen

            ! Determine the difference in charge with the ground state
            dZ = Atom%stage(iterm) - Atom%stage(1)

            ! Calculate the argument of the exponential, ionization
            ! potential with Debey correction
            arg = (debey(ilevel)*dby - dE)*Tin/kb

            !
            ! Numerator of the population solution

            ! Negative exponential argument
            if (arg.lt.0d0) then

              ! Make argument positive
              arg = -arg

              ! Get inverse exponential
              Atom%populte(ilevel,iz) = gi0*diexp(arg)

            ! Positive exponential argument
            else

              ! Get exponential controling under/overflows
              Atom%populte(ilevel,iz) = gi0*ddexp(arg)

            end if ! Exponential arguments

            ! For each height
            do ii=1,dZ

              ! Apply denominator
              Atom%populte(ilevel,iz) = Atom%populte(ilevel,iz)/C2

            end do ! Heights

            ! Accumulate for the normalization equation
            S = S + Atom%populte(ilevel,iz)

          end do ! J levels
        end do ! Terms

        ! Get the ground level population from the accumulated factor
        Atom%populte(1,iz) = Atom%n(iz)/S

        ! For each level
        do ilevel=2,Atom%nlevel

          ! Determine its population from the population of the
          ! ground state
          Atom%populte(ilevel,iz) = Atom%populte(ilevel,iz)* &
                                    Atom%populte(1,iz)

        end do ! Levels
      end do ! Heights

      ! Deallocate debey array
      deallocate(debey)

      ! Check if everything is fine
      call control

      return

      end subroutine LTEPopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute LTE populations for a given LTE line in a given model
      !! atmosphere\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!                   T(double): Temperature\n
      !!                  ne(double): Electron number density\n
      !!                  C0(double): Constant for non-atomic part
      !!                              Saha function\n
      !!                  C1(double): Constant for Debey factor\n
      !!                 iz(integer): Height index\n
      !!               nstg(integer): Number of stages there is
      !!                              information on partition
      !!                              function\n
      !!               Ei(double(:)): Ionization energy data\n
      !!               pf(double(:)): Partition function data\n
      !!               debey(double): Debey correction
      subroutine LTEPopu_line(LTEline,T,ne,C0,C1,iz,nstg,Ei,pf,debey)

      ! I/O

      type(LTEline_class), intent(inout):: LTEline
      integer, intent(in):: iz,nstg
      double precision, intent(in):: T,ne,debey,C0,C1
      double precision, dimension(:), intent(in):: Ei,pf

      ! Local

      integer:: dZ

      double precision:: C2,dby,Tin,dE,gif,arg,ikT,U
      double precision, dimension(1):: frc


      !
      ! Compute LTE populations
      !

      ! Inverse of temperature
      Tin = 1d0/T
      ! Multiplicative factor for Debey correction
      dby = C1*sqrt(ne*Tin)
      ! Non atomic part of the Saha function
      C2 = 1d0/(.5d6*ne*((C0*Tin)**(1.5d0)))
      ! Scaled to Kb
      ikT = fktoJ/kb/T

      ! Determine the difference in charge with the ground state
      dZ = LTEline%stage

      ! Partition function
      U = pf(LTEline%stage)

      ! Calculate ionization fraction
      call getfrc(nstg,pf,Ei,T,ne,LTEline%stage,frc)

      !
      ! Lower level
      !

      ! Energy and degeneration
      dE = LTEline%El*fktoJ
      gif = (2d0*LTEline%Jl + 1d0)*frc(1)

      ! Calculate the argument of the exponential, ionization
      ! potential with Debey correction
      arg = (debey*dby - dE)*Tin/kb - U

      ! If negative argument
      if (arg.lt.0d0) then

        ! Make positive and apply inverse exponential
        arg = -arg
        LTEline%nl(iz) = gif*diexp(arg)

      ! If positive argument
      else

        ! Apply exponential checking for under/overflow
        LTEline%nl(iz) = gif*ddexp(arg)

      end if ! Argument sign

      !
      ! Upper level
      !

      ! Energy and degeneration
      dE = LTEline%Eu*fktoJ
      gif = (2d0*LTEline%Ju + 1d0)*frc(1)

      ! Calculate the argument of the exponential, ionization
      ! potential with Debey correction
      arg = (debey*dby - dE)*Tin/kb - U

      ! If negative argument
      if (arg.lt.0d0) then

        ! Make positive and apply inverse exponential
        arg = -arg
        LTEline%nu(iz) = gif*diexp(arg)

      ! If positive argument
      else

        ! Apply exponential checking for under/overflow
        LTEline%nu(iz) = gif*ddexp(arg)

      end if ! Argument sign

      return

      end subroutine LTEPopu_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage hydrogen populations\n
      !> Set the hydrogen populations in the model atom to the ones in
      !! the model atmosphere atmosphere, distributed in LTE among
      !! sublevels of the given principal levels\n
      !!   Atom(Atom_class): Structure with atomic data\n
      !!   Atmo(Atmo_class): Structure with atmospheric data\n
      !!  LJstruct(logical): If the atomic model has LJ structure
      subroutine InitHpopu(Atom,Atmo,LJstruct)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), intent(inout):: Atom
      logical, intent(out):: LJstruct

      ! Local

      integer, parameter:: nn = 6

      integer:: iz,it,iJ,il,jl
      integer, dimension(Atom%nlevel):: loc

      double precision:: ikkb,ikT,deg,exu
      double precision, dimension(nn):: Emin,Emax,W,TEmin,ig0,SS,nlv


      ! Constant
      ikkb = fktoJ/kb

      !
      ! Definition of a H atom energies
      !              0     1     2      3      4        p+
      Emin = (/ -1.0d0,0.5d0,0.9d0,1.00d0,1.04d0,1.09677d0 /)
      Emax = (/  0.5d0,0.9d0,1.0d0,1.04d0,1.06d0,1.00000d9 /)

      ! Initialize populations
      Atom%popu = 0d0

      ! Initialize location, weight, TEmin, number of levels, and
      ! LJ structure
      loc = 0
      W = 0d0
      TEmin = 1d99
      ig0 = 1d0
      nlv = 0
      LJstruct = .False.

      ! Initialize level index
      il = 0

      !
      ! Check to which region each transition is part of
      !

      ! For each term
      do it=1,Atom%nmulti

        ! For each level within the term
        do iJ=1,Atom%nJ(it)

          ! Advance level index
          il = il + 1

          ! Look in all ranges
          do jl=1,nn

            ! Check if energy in the range
            if (Atom%FSfreq(iJ,it).gt.Emin(jl).and. &
                Atom%FSfreq(iJ,it).lt.Emax(jl)) then

              ! Get degeneration of the level
              deg = 2d0*Atom%rJval(iJ,it) + 1d0

              ! Check true minimum energy for this range
              if (Atom%FSfreq(iJ,it).lt.TEmin(jl)) then

                ! New minimum energy and degeneration factor
                TEmin(jl) = Atom%FSfreq(iJ,it)
                ig0(jl) = 1d0/deg

              end if ! New minimum energy

              ! Set location, add to weight, and exit
              loc(il) = jl
              nlv(jl) = nlv(jl) + 1
              W = W + deg
              exit

            end if

          end do ! Atmos ranges
        end do ! J levels
      end do ! Terms

      ! For each height node
      do iz=1,nZ

        ! 1/kT
        ikT = ikkb/Atmo%T(iz)

        ! Initialize level index
        il = 0

        ! Initialize sum
        SS = 0d0

        ! For each term
        do it=1,Atom%nmulti

          ! For each level within the term
          do iJ=1,Atom%nJ(it)

            ! Advance level index
            il = il + 1

            ! Get location
            jl = loc(il)

            ! If no location, skip
            if (jl.le.0) cycle

            ! Exponential argument
            exu = Atom%FSfreq(iJ,it) - TEmin(jl)

            ! Skip first level in the term
            if (exu.lt.TINYO) cycle

            ! Complete exponential
            exu = exu*ikT
            exu = diexp(exu)*ig0(jl)

            ! Add to population and to the sum
            Atom%popu(il,iz) = (2d0*Atom%rJval(iJ,it) + 1d0)*exu
            SS(jl) = SS(jl) + Atom%popu(il,iz)

          end do ! J level
        end do ! Atomic term

        ! Initialize level index
        il = 0

        ! For each term
        do it=1,Atom%nmulti

          ! For each level within the term
          do iJ=1,Atom%nJ(it)

            ! Advance level index
            il = il + 1

            ! Get location
            jl = loc(il)

            ! If no location, skip
            if (jl.le.0) cycle

            ! If only one level, get the whole population
            if (nlv(jl).le.1) then

              Atom%popu(il,iz) = Atmo%nH(iz,jl)

            ! If several levels
            else

              ! Set structure to True
              if (.not.LJstruct) LJstruct = .True.

              ! Exponential
              exu = Atom%FSfreq(iJ,it) - TEmin(jl)

              ! If first level
              if (exu.lt.TINYO) then

                ! Get population as fraction from total and update
                ! sum
                Atom%popu(il,iz) = Atmo%nH(iz,jl)/SS(jl)
                SS(jl) = Atom%popu(il,iz)

              ! Not first level
              else

                ! Get population from first level
                Atom%popu(il,iz) = Atom%popu(il,iz)*SS(jl)

              end if ! First level or not
            end if ! One/many levels

          end do ! J level
        end do ! Atomic term

      end do ! Heights

      ! Check if everything is fine
      call control

      return

      end subroutine InitHpopu

!#####################################################################
!#####################################################################
!#####################################################################

      end module initpopuaux_mod
