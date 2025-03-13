      !> Reading of partition function and abundance
      module rpfa_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     12/07/2022
!  Last version:
!     19/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     19/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  rParfunAbund
!    Read data on the partition function from the specified file.
!  Read data on abundances from the specified file or initialize them
!  to the hard-coded tabulation
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read data on the partition function from the specified file.
      !! Read data on abundances from the specified file or
      !! initialize them to the hard-coded tabulation\n
      !!  Input(Input_class): Structure with configuration data\n
      !!    Atmo(Atmo_class): Structure with atmospheric data
      subroutine rParfunAbund(Input,Atmo)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Atmo_class), intent(inout):: Atmo

      ! Local

      character(len=1):: cdump
      character(len=2):: cbuff

      integer:: iele,iele0,nele,ios


      ! Routine name
      urou = 'rParfunAbund'

      ! Number of elements
      Atmo%nele = recallnumber(0)

      !
      ! Read the partition function data from the file
      !

      ! If specific input
      if (trim(Input%pf).ne.'NONE') then

        ! Open file
        open (200,file=trim(Input%pf), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        ! Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = ' - Read partition functions from '// &
                 trim(Input%pf)
          call verbose

        end if ! Master

      ! Default input
      else

        ! Open default file
        open (200,file=trim(Input%resource)//'partfunc', &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        ! Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = ' - Read partition functions from '// &
                 trim(Input%resource)//'partfunc'
          call verbose

        end if ! Master

      end if ! Specified file

      ! Number of temperatures for the tabulated data
      read(200,err=1100) Atmo%NT

      ! Temperature inputs
      allocate(Atmo%pT(Atmo%NT))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%pT)

      ! Read the temperature
      read(200,err=1100) Atmo%pT

      ! Allocate elements
      allocate(Atmo%ele(Atmo%nele))


      !
      ! For each element get the partition function tabulation
      !

      ! For each element
      do iele0=1,Atmo%nele

        ! Read first character
        read(200,err=1100) cdump
        cbuff(1:1) = cdump

        ! Read second character
        read(200,err=1100) cdump

        ! If the second character was a space, shift the characters in
        ! the string
        if (cdump.eq.' ') then
          cbuff(2:2) = cbuff(1:1)
          cbuff(1:1) = cdump
        else
          cbuff(2:2) = cdump
        end if

        ! Get index for this atom
        iele = atom_char2index(cbuff)

        ! Memory count
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(iele))

        ! Save element name
        Atmo%ele(iele)%Element = cbuff

        ! Unnecessary character
        read(200,err=1100) cdump

        ! Read the number of stages that follows
        read(200,err=1100) Atmo%ele(iele)%nstg

        ! Allocate the partition function and ionization energies
        allocate(Atmo%ele(iele)%pf(Atmo%nT,Atmo%ele(iele)%nstg))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(iele)%pf)
        allocate(Atmo%ele(iele)%Ei(Atmo%ele(iele)%nstg))
        MRAMc = MRAMc + 1d-6*sizeof(Atmo%ele(iele)%Ei)

        ! Read the pf data
        read(200,err=1100) Atmo%ele(iele)%pf
        Atmo%ele(iele)%pf = log(Atmo%ele(iele)%pf)

        ! Get ionization energies
        read(200,err=1100) Atmo%ele(iele)%Ei
        Atmo%ele(iele)%Ei = Atmo%ele(iele)%Ei*1d-5

      end do ! Elements in the partition function file

      ! Close partition function file
      close(200)

      !
      ! Abundances
      !

      ! Allocate abundance
      allocate(Atmo%abund(99))
      MRAMc = MRAMc + 1d-6*sizeof(Atmo%abund)

      ! Set default abundance for each element
      do iele=1,99
        Atmo%abund(iele) = recallabund_ind(iele)
      end do

      ! If there is an abundance file
      if (trim(Input%abund).ne.'NONE') then

        ! Open file
        open (200,file=trim(Input%abund), &
              status='unknown', iostat=ios, err=1001, &
              access='stream', action='read', &
              form='unformatted')

        ! Read number of elements in file
        read(200,err=1101) nele

        ! For each element in input
        do iele0=1,nele

          ! Read index and abundance
          read(200,err=1101) iele
          read(200,err=1101) Atmo%abund(iele)

          ! Transform to ratio
          Atmo%abund(iele) = 1d1**(Atmo%abund(iele) - 12d0)

        end do ! Elements in input

        ! Close file
        close(200)

        ! Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = ' - Read abundances from '// &
                 trim(Input%abund)
          call verbose

        end if ! Master
      end if ! Specified an abundance file

      return

1000  umsg = 'Error opening pf file'
      call aborted
      return
1100  umsg = 'Error reading pf file'
      close(200)
      call aborted
      return
1001  umsg = 'Error opening abundance file'
      call aborted
      return
1101  umsg = 'Error reading abuncance file'
      close(200)
      call aborted
      return

      end subroutine rParfunAbund

!#####################################################################
!#####################################################################
!#####################################################################

      end module rpfa_mod
