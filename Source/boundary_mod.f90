      !> Boundary conditions
      module boundary_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Significantly simplified what is done
!                             in the subroutines, now expecting more
!                             processed arguments (TdPA)
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
!  top
!    Define Stokes parameters at the upper boundary
!
!  topI
!    Define intensity at the upper boundary
!
!  bottom
!    Define Stokes parameters at the bottom boundary
!
!  bottomI
!    Define intensity at the bottom boundary
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use planck_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Define Stokes parameters at the upper boundary\n
      !!      if0(integer): Lower limit index for frequency\n
      !!      if1(integer): Upper limit index for frequency\n
      !!  Stk(double(:,:)): Stokes parameters at the upper boundary
      subroutine top(if0,if1,Stk)

      ! I/O

      integer, intent(in):: if0,if1
      double precision, dimension(0:3,if0:if1), intent(out):: Stk

      ! Vacuum on top
      Stk = 0d0

      end subroutine top

!#####################################################################
!#####################################################################
!#####################################################################

      !> Define intensity at the upper boundary\n
      !!    if0(integer): Lower limit index for frequency\n
      !!    if1(integer): Upper limit index for frequency\n
      !!  Stk(double(:)): Intensity at upper boundary
      subroutine topI(if0,if1,Stk)

      ! I/O

      integer, intent(in):: if0,if1
      double precision, dimension(if0:if1),intent(out):: Stk

      ! Vacuum on top
      Stk = 0d0

      end subroutine topI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Define Stokes parameters at the bottom boundary\n
      !!  omega(double(:)): Frequency array\n
      !!         T(double): Temperature\n
      !!      vfac(double): Doppler shift factor\n
      !!      if0(integer): Lower limit index for frequency\n
      !!      if1(integer): Upper limit index for frequency\n
      !!  Stk(double(:,:)): Stokes parameters at bottom boundary
      subroutine bottom(omega,T,vfac,if0,if1,Stk)

      ! I/O

      integer, intent(in):: if0,if1
      double precision, intent(in):: T,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,if0:if1), intent(out):: Stk

      ! Local

      integer:: ifreq


      ! Initialize
      Stk = 0d0

      ! For every frequency
      do ifreq=if0,if1

        ! Make the boundary planckian [k-units]
        Stk(0,ifreq) = planck(omega(ifreq)*vfac,T)

      end do ! Frequencies

      end subroutine bottom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Define intensity at the bottom boundary\n
      !!  omega(double(:)): Frequency array\n
      !!         T(double): Temperature\n
      !!      vfac(double): Doppler shift\n
      !!      if0(integer): Lower limit index for frequency\n
      !!      if1(integer): Upper limit index for frequency\n
      !!    Stk(double(:)): Intensity at bottom boundary
      subroutine bottomI(omega,T,vfac,if0,if1,Stk)

      ! I/O

      integer, intent(in):: if0,if1
      double precision, intent(in):: T,vfac
      double precision, dimension(:), intent(in)::  omega
      double precision, dimension(if0:if1), intent(out):: Stk

      ! Local

      integer:: ifreq


      ! For each frequency
      do ifreq=if0,if1

        ! Make the boundary planckian [k-units]
        Stk(ifreq) = planck(omega(ifreq)*vfac,T)

      end do ! Frequencies

      end subroutine bottomI

!#####################################################################
!#####################################################################
!#####################################################################

      end module boundary_mod
