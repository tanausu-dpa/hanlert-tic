      !> Diagonalization manager
      module diagon_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     12/12/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     12/12/2023:    V3.0.2 - Added diagon_B0 subroutine (TdPA)
!
!     10/25/2022:    V3.0.1 - The relevant atomic variables are
!                             allocated here now (TdPA)
!                           - Limited the ranges of the height do
!                             loop (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     12/17/2019:    V1.0.2 - Passing another parameter to PB (TdPA)
!
!     10/27/2017:    V1.0.1 - Initializing evec and eval (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!  diagon:
!    This subroutine calls the routine that calculates the energy
!    eigenvalues and eigenvectors for a multiplet, in the
!    magnetic-field regime of the incomplete Paschen-Back effect.
!    Diagonality with respect to M is exploited for block-
!    diagonalization
!
!  diagon_B0:
!    This subroutine calls the routine that initializes some
!    multi-term quantities in absence of magnetic fields.
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

      !> Runs over all the heights and terms to call the
      !! diagonalization routine\n
      !!       Atom(Atom_class): Structure with the atomic data\n
      !!   Bfield(Bfield_class): Structure with magnetic field data\n
      !!          mode(integer): Type of Zeeman effect\n
      !!     Flgsg(Fctsg_class): Structure with factorials and
      !!                         signs
      subroutine diagon(Atom,Bfield,mode,Flgsg)

      ! I/O

      type(Atom_class),intent(inout):: Atom
      type(Bfield_class),intent(in):: Bfield
      type(Fctsg_class),intent(in):: Flgsg
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

      ! Initialize evec and eval
      Atom%eval = 0d0
      Atom%evec = 0d0

      ! For each height
      do iz=Rz0,Rz1

        ! Calculate Larmor frequency in k units
        larmork = B2LK*Bfield%Bstrength(iz)

        ! And term
        do iterm=1,Atom%nMulti

          call PB(iz,iterm,larmork,mode,Flgsg,Atom)

        enddo ! Term
      end do ! Height

      end subroutine diagon

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize diagonalization variables for B=0\n
      !!       Atom(Atom_class): Structure with the atomic data
      subroutine diagon_B0(Atom)

      ! I/O

      type(Atom_class),intent(inout):: Atom

      ! Local

      integer:: iterm

      ! For each term
      do iterm=1,Atom%nMulti

       call PB0(iterm,Atom)

      enddo ! Term

      end subroutine diagon_B0

!#####################################################################
!#####################################################################
!#####################################################################

      end module diagon_mod
