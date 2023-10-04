      !> Reading of partition function and abundance
      module rpfa_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     07/12/2022
!  Last version:
!     07/12/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     07/12/2022:    V3.0.0 - First version (TdPA)
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
!  This subroutine reads the partition function data and sets up
!  abundances
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

      !> Reads the partition function file and sets up abundances.\n
      !!  Input(Input_class): Structure with settings data\n
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

        open (200,file=trim(Input%pf), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        if (gpid.eq.0) then
          umsg = ' - Read partition functions from '// &
                 trim(Input%pf)
          call verbose
        end if

      ! Default input
      else

        open (200,file=trim(Input%resource)//'partfunc', &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        if (gpid.eq.0) then
          umsg = ' - Read partition functions from '// &
                 trim(Input%resource)//'partfunc'
          call verbose
        end if

      end if

      ! Number of temperature input for the data table
      read(200,err=1100) Atmo%NT

      ! Temperature inputs
      allocate(Atmo%pT(Atmo%NT))

      ! Read the Temperature
      read(200,err=1100) Atmo%pT

      !
      ! Allocate elements
      !
      allocate(Atmo%ele(Atmo%nele))


      !
      ! For each element get the partition function table
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

        ! Get index
        iele = atom_char2index(cbuff)

        ! Save name
        Atmo%ele(iele)%Element = cbuff

        ! Unnecessary character
        read(200,err=1100) cdump
        ! Read the number of stages that follows
        read(200,err=1100) Atmo%ele(iele)%nstg

        ! Allocate the partition function and ionization energies
        allocate(Atmo%ele(iele)%pf(Atmo%nT,Atmo%ele(iele)%nstg))
        allocate(Atmo%ele(iele)%Ei(Atmo%ele(iele)%nstg))

        ! Read the pf data
        read(200,err=1100) Atmo%ele(iele)%pf
        Atmo%ele(iele)%pf = log(Atmo%ele(iele)%pf)

        ! Get ionization energies
        read(200,err=1100) Atmo%ele(iele)%Ei
        Atmo%ele(iele)%Ei = Atmo%ele(iele)%Ei*1d-5

      end do ! Elements in the partition function file

      close(200)

      !
      ! Abundances
      !

      ! Allocate abundance
      allocate(Atmo%abund(99))

      ! Set default for each element
      do iele=1,99
        Atmo%abund(iele) = recallabund_ind(iele)
      end do

      ! Now read from file if any
      if (trim(Input%abund).ne.'NONE') then

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
          Atmo%abund(iele) = 1d1**(Atmo%abund(iele) - 12d0)

        end do

        close(200)

        if (gpid.eq.0) then
          umsg = ' - Read abundances from '// &
                 trim(Input%abund)
          call verbose
        end if

      end if

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
