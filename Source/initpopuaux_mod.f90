      !> Initializations for populations
      module initpopuaux_mod
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
!     08/07/2023 V3.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.1 - Added LTEPopu_line (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o The reading of populations needs to
!                                be branched depending on the type of
!                                run. Added the needed logic for the
!                                1.5D case.
!                              o Added return after abortion calls.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     09/11/2020:    V1.2.6 - Added option to read departure
!                             coefficients (TdPA)
!
!     03/05/2020:    V1.2.5 - Added InitHpopu to manage hydrogen
!                             populations. Takes the populations from
!                             the atmospheric file and distribute
!                             them in LTE among subleves of a given
!                             principal level (TdPA)
!
!     11/19/2019:    V1.2.4 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     10/31/2019:    V1.2.3 - Improved debug message in rpopu (TdPA)
!
!     07/30/2019:    V1.2.2 - The Debey shielding was wrong by a
!                             factor 1d3, because m^-3 was expected
!                             for electron density, and I was giving
!                             cm^-3 (TdPA)
!
!     03/18/2019:    V1.2.1 - Introduced exponentials with argument
!                             control (TdPA)
!
!     02/20/2019:    V1.2.0 - New verbosity (TdPA)
!                           - Specified used parameters (TdPA)
!                           - Now uses unit 200 (TdPA)
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
!    Auxiliar for initpopu_mod
!
!  rPopu:
!    Reads populations from a file
!
!  LTEPopu
!    Computes the LTE population balance
!
!  LTEPopu_line
!    Computes the LTE populations for an LTE line
!
!  InitHpopu
!    Sets H atom populations to the ones in the atmosphere,
!  distributed in LTE among sublevels of the given principal levels
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

      !> Reads a file with population data.\n
      !!   filename(character(:)): Name of the file to read\n
      !!         Atom(Atom_class): Structure with the atomic data
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

        umsg = 'Reading of populations not implemented in inversion'
        call aborted
        return

      end if ! Type of run

      !
      ! Open the file with the populations
      !
      open (200,file=trim(filename), &
            status='unknown', iostat=ios, err=1000, access='stream', &
            action='read', form='unformatted')

      read (200,err=1100) fileid

      ! Only 1.5D
      if (run_mode.eq.1) then

        ! Check that it has the correct population file ID
        if (fileid.ne.'2D') then

          umsg = 'The ID of the population file for '// &
                 'Atom'//Atom%Element//' is not valid.'
          call aborted
          close(200)
          return

        end if

        read (200,err=1100) fileid

      end if

      ! Check that it has the correct population file ID
      if (fileid.ne.'bp'.and.fileid.ne.'bd') then

        umsg = 'The ID of the population file for '// &
               'Atom'//Atom%Element//' is not valid.'
        call aborted
        return

      end if

      read (200,err=1100) nZ_popu

      ! Check if it has the correct number of nodes
      if (nZ.ne.nZ_popu) then

        umsg = 'Node points in atmosphere and '// &
               'populations do not coincide.'
        call aborted
        close(200)
        return

      end if

      read (200,err=1100) nL_popu

      ! Check if it has the correct number of levels
      if(Atom%nlevel .ne. nL_popu)then

        write(umsg,'(A,1x,i4,1x,A,1x,i4)') &
              'Number of levels in atom and populations'// &
              ' does not coincide:',Atom%nlevel,'!=',nL_popu
        call aborted
        close(200)
        return

      end if

      ! Only 1.5D
      if (run_mode.eq.1) then

        ! Jump to current column
        offset = nz*8*nL_popu*(icoords(3)-1)
        call fseek(200,offset,1)

      end if

      ! If proper populations
      if (fileid.eq.'bp') then

        allocate(Atom%popu(Atom%nlevel,nz))

        do iz=1,nZ

          read (200,err=1100) buff

          ! Store read data into population array
          Atom%popu(:,iz) = buff

          ! Calculate the total population from this input
          Atom%n(iz) = sum(Atom%popu(:,iz))

        end do

      ! If departure coefficients
      else if (fileid.eq.'bd') then

        allocate(Atom%depar(Atom%nlevel,nz))

        do iz=1,nZ

          read (200,err=1100) buff

          ! Store read data into population array
          Atom%depar(:,iz) = buff

        end do

      end if

      close (200)

      ! Check if everything is fine
      call control
      if (laborted) return

      if (gpid.eq.0) then
        if (fileid.eq.'bp') then
          umsg = ' - Populations '//trim(filename)//' read'
        else if (fileid.eq.'bd') then
          umsg = ' - Departure coefficients '//trim(filename)//' read'
        end if
        call verbose
      end if

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

      !> Computes LTE populations for a given atom and atmosphere\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!    Atmo(Atmo_class): Structure with atmospheric data
      subroutine LTEPopu(Atom,Atmo)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo

      ! Local

      integer:: ilevel,iterm,iJ,ii,iim,iz,dZ

      double precision:: C0, C1, C2, zm, dby, Tin, S, dE, gi0, arg
      double precision, dimension(:), allocatable:: debey


      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      ! 1d3 sqrt(cm^-3 -> m^-3) electron density below
      C1 = sqrt(8d0*PI/kb)*((qel*qel/pi4eps0)**(1.5d0))*1d3


      !
      ! Calculate Debey correction to the ionization potential
      !

      ! Allocate the correction
      allocate(debey(Atom%nlevel))
      debey = 0d0

      ! Initialize the index
      ilevel = 1

      ! Run over all the levels using the term-J indexes
      do iterm=1,Atom%nMulti
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

            debey(ilevel) = debey(ilevel) + zm
            zm = zm + 1

          end do

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

        ! For each level running through term and J indexes
        do iterm=1,Atom%nMulti
          do iJ=1,Atom%nJ(iterm)

            ! Skip ground level
            if (iterm.eq.1.and.iJ.eq.1) cycle

            ! Determine the index and the ionization potential
            ilevel = ilevel + 1
            dE = (Atom%FSfreq(iJ,iterm) - Atom%FSfreq(1,1))*fktoJ

            ! If it is the hardwired hydrogen, use the degeneration
            ! variable
            if (Atom%cust) then

              gi0 = Atom%deg(ilevel)/Atom%deg(1)

            ! If it is a normal atom, use the J values
            else

              gi0 = (2d0*Atom%rJval(iJ,iterm) + 1d0)/ &
                    (2d0*Atom%rJval(1,1) + 1d0)

            end if

            ! Determine the difference in charge with the ground state
            dZ = Atom%stage(iterm) - Atom%stage(1)

            ! Calculate the argument of the exponential, ionization
            ! potential with Debey correction
            arg = (debey(ilevel)*dby - dE)*Tin/kb

            ! Numerator of the population solution
            if (arg.lt.0d0) then
              arg = -arg
              Atom%populte(ilevel,iz) = gi0*diexp(arg)
            else
              Atom%populte(ilevel,iz) = gi0*ddexp(arg)
            end if

            ! Denominator
            do ii=1,dZ
              Atom%populte(ilevel,iz) = Atom%populte(ilevel,iz)/C2
            end do

            ! Accumulate for the normalization equation
            S = S + Atom%populte(ilevel,iz)

          end do ! J levels
        end do ! Terms

        ! Get the ground level population from the accumulated
        ! factor
        Atom%populte(1,iz) = Atom%n(iz)/S

        ! Determine the population of each level from the
        ! population of the ground state
        do ilevel=2,Atom%nlevel
          Atom%populte(ilevel,iz) = Atom%populte(ilevel,iz)* &
                                    Atom%populte(1,iz)
        end do

      end do ! Heights

      ! Deallocate debey vector
      deallocate(debey)

      ! Check if everything is fine
      call control

      return

      end subroutine LTEPopu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes LTE populations for the pair of levels of a LTE
      !! transition\n
      !!  LTEline(Input_class): Structure with LTE line data\n
      !!             T(dfloat): Temperature\n
      !!            ne(dfloat): Electron number density\n
      !!            C0(dfloat): Constant for non-atomic part Saha
      !!                        function\n
      !!            C1(dfloat): Constant for Debey factor\n
      !!           iz(integer): Height index of the call\n
      !!         nstg(integer): Number of stages there is
      !!                        information on partition function\n
      !!         Ei(dfloat(:)): Ionization energy\n
      !!         pf(dfloat(:)): Partition function\n
      !!         debey(dfloat): Debey correction
      subroutine LTEPopu_line(LTEline,T,ne,C0,C1,iz,nstg,Ei,pf,debey)

      ! I/O
      type(LTEline_class):: LTEline
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

      !
      ! Calculate ionization fraction
      !
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

      ! Population relative to total
      if (arg.lt.0d0) then
        arg = -arg
        LTEline%nl(iz) = gif*diexp(arg)
      else
        LTEline%nl(iz) = gif*ddexp(arg)
      end if

      !
      ! Upper level
      !

      ! Energy and degeneration
      dE = LTEline%Eu*fktoJ
      gif = (2d0*LTEline%Ju + 1d0)*frc(1)

      ! Calculate the argument of the exponential, ionization
      ! potential with Debey correction
      arg = (debey*dby - dE)*Tin/kb - U

      ! Population relative to total
      if (arg.lt.0d0) then
        arg = -arg
        LTEline%nu(iz) = gif*diexp(arg)
      else
        LTEline%nu(iz) = gif*ddexp(arg)
      end if

      return

      end subroutine LTEPopu_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Manage hydrogen populations\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!        LJstruct(logical): If atomic model had LJ structure
      subroutine InitHpopu(Atom,Atmo,LJstruct)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), intent(inout):: Atom
      logical, intent(out):: LJstruct

      ! Parameter
      integer, parameter:: nn = 6

      ! Local

      integer:: iz,it,iJ,il,jl
      integer, dimension(Atom%nlevel):: loc

      double precision:: ikkb, ikT, deg, exu
      double precision, dimension(nn):: Emin, Emax, W, TEmin, ig0, SS
      double precision, dimension(nn):: nlv

      ! Constant
      ikkb = fktoJ/kb

      !
      ! Definitions of an H atom
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

      ! Initialize level counter
      il = 0

      !
      ! Check to which region each transition is part of
      !
      do it=1,Atom%nmulti
        do iJ=1,Atom%nJ(it)

          ! Advance level index
          il = il + 1

          ! Look for the range
          do jl=1,nn

            ! Check energies
            if (Atom%FSfreq(iJ,it).gt.Emin(jl).and. &
                Atom%FSfreq(iJ,it).lt.Emax(jl)) then

              ! Get degeneration
              deg = 2d0*Atom%rJval(iJ,it) + 1d0

              ! Check true minimum energy
              if (Atom%FSfreq(iJ,it).lt.TEmin(jl)) then
                TEmin(jl) = Atom%FSfreq(iJ,it)
                ig0(jl) = 1d0/deg
              end if

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

        ! For each atomic level
        do it=1,Atom%nmulti
          do iJ=1,Atom%nJ(it)

            ! Advance level index
            il = il + 1

            ! Get location
            jl = loc(il)

            ! If no location, skip
            if (jl.le.0) cycle

            ! Exponential
            exu = Atom%FSfreq(iJ,it) - TEmin(jl)

            ! Skip first level
            if (exu.lt.TINYO) cycle

            exu = exu*ikT
            exu = diexp(exu)*ig0(jl)

            ! Add to population and sum
            Atom%popu(il,iz) = (2d0*Atom%rJval(iJ,it) + 1d0)*exu
            SS(jl) = SS(jl) + Atom%popu(il,iz)

          end do ! J level
        end do ! Atomic term

        ! Initialize level index
        il = 0

        ! For each atomic level
        do it=1,Atom%nmulti
          do iJ=1,Atom%nJ(it)

            ! Advance level index
            il = il + 1

            ! Get location
            jl = loc(il)

            ! If no location, skip
            if (jl.le.0) cycle

            ! If only one level
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

                Atom%popu(il,iz) = Atmo%nH(iz,jl)/SS(jl)
                SS(jl) = Atom%popu(il,iz)

              ! Not first level
              else

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
