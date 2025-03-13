      !> Diagonalization manager
      module diagon_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     18/04/2017
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
!  Data:
!
!  diagon
!    Call the actual diagonalization for each height in the model
!  atmosphere
!
!  diagon_B0:
!    Initialize diagonalization related variables when there is no
!  magnetic field
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use parameters_mod , only : B2LK
      use pb_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Call the actual diagonalization for each height in the model
      !! atmosphere\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data\n
      !!         mode(integer): Type of Zeeman effect\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols
      subroutine diagon(Atom,Bfield,mode,Flgsg)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: mode

      ! Local

      integer:: iz, iterm

      double precision:: larmork


      ! Allocate eigenvalues and eigenvectors
      allocate(Atom%eval(Atom%nJmax,Atom%nMmax, &
                         Atom%nMulti,Rz0:Rz1))
      allocate(Atom%evec(Atom%nJmax,Atom%nJmax, &
                         Atom%nMmax, &
                         Atom%nMulti,Rz0:Rz1))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Atom%eval)
      MRAMc = MRAMc + 1d-6*sizeof(Atom%evec)

      ! Initialize evec and eval
      Atom%eval = 0d0
      Atom%evec = 0d0

      ! For each height
      do iz=Rz0,Rz1

        ! Calculate Larmor frequency in k units
        larmork = B2LK*Bfield%Bstrength(iz)

        ! For each term
        do iterm=1,Atom%nMulti

          ! Diagonalize
          call PB(iz,iterm,larmork,mode,Flgsg,Atom)

        enddo ! Term
      end do ! Height

      end subroutine diagon

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize diagonalization related variables when there is no
      !! magnetic field
      !!  Atom(Atom_class): Structure with atomic data
      subroutine diagon_B0(Atom)

      ! I/O

      type(Atom_class),intent(inout):: Atom

      ! Local

      integer:: iterm


      ! For each term
      do iterm=1,Atom%nMulti

        ! Call initialization of diagonalization variables without
        ! magnetic field
        call PB0(iterm,Atom)

      enddo ! Term

      end subroutine diagon_B0

!#####################################################################
!#####################################################################
!#####################################################################

      end module diagon_mod
