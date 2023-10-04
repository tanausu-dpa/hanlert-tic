      !> Boundary conditions
      module boundary_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     11/24/2022 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/24/2022:    V3.0.2 - Skip computing the velocity shift if
!                             the dynamic flag is off (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case Atmo%v has
!                             changed to Atmo%vx,%vy, and %vz (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  top:
!    This subroutine initializes the top boundary
!
!  bottom:
!    This subroutine initializes the bottom boundary
!
!  topI:
!    This subroutine initializes the top boundary for only intensity
!
!  bottomI:
!    This subroutine initializes the bottom boundary for only
!  intensity
!
!  Planckian function (k-units)
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

      !> Returns Stokes parameters at the upper boundary\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !!     Stk(dfloat(:,:)): Stokes parameters at upper boundary
      subroutine top(MPID,Stk)

      ! I/O

      type(MPI_class), intent(in):: MPID
      double precision, dimension(0:3,MPID%if0(pid):MPID%if1(pid)), &
                        intent(out):: Stk

      Stk = 0d0

      end subroutine top

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns intensity at the upper boundary\n
      !!    MPID(MPI_class): Structure with MPI data\n
      !!     Stk(dfloat(:)): Intensity at upper boundary
      subroutine topI(MPID,Stk)

      ! I/O

      type(MPI_class), intent(in):: MPID
      double precision, dimension(MPID%if0(pid):MPID%if1(pid)), &
                        intent(out):: Stk

      Stk = 0d0

      end subroutine topI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns Stokes parameters at the upper boundary\n
      !!   omega(dfloat(:)): Frequency array\n
      !!          T(dfloat): Temperature\n
      !!         vx(dfloat): Velocity along X\n
      !!         vy(dfloat): Velocity along Y\n
      !!         vz(dfloat): Velocity along Z\n
      !!         mu(dfloat): Cosine of polar angle\n
      !!        mux(dfloat): Cosine of azimuth\n
      !!     muy_in(dfloat): Sign of sin of azimuth\n
      !!    MPID(MPI_class): Structure with MPI data\n
      !!   Stk(dfloat(:,:)): Stokes parameters at bottom boundary
      subroutine bottom(omega,T,vx,vy,vz,mu,mux,muy_in,MPID,Stk)

      ! I/O

      type(MPI_class), intent(in):: MPID
      double precision, intent(in):: T,mu,mux,muy_in,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(0:3,MPID%if0(pid):MPID%if1(pid)), &
                        intent(out):: Stk

      ! Local

      integer:: ifreq

      double precision:: muy,ct,st,cc,sc,vfac


      ! Initialize
      Stk = 0d0

      ! Calculate Doppler shift factor
      if (dyn) then

        ct = mu
        st = sqrt(1d0 - ct*ct)
        cc = mux
        muy = muy_in
        if (abs(muy).lt.1d-8) muy = 1d0
        sc = muy*sqrt(1d0 - cc*cc)/abs(muy)

        vfac = 1d0 - vx*st*cc - vy*st*sc - vz*ct

      ! Static
      else

        vfac = 1d0

      end if

      ! Make the boundary planckian
      do ifreq=MPID%if0(pid),MPID%if1(pid)

        Stk(0,ifreq) = planck(omega(ifreq)*vfac,T)

      end do

      end subroutine bottom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Returns intensity at the upper boundary\n
      !!  omega(dfloat(:)): Frequency array\n
      !!         T(dfloat): Temperature\n
      !!        vx(dfloat): Velocity along X\n
      !!        vy(dfloat): Velocity along Y\n
      !!        vz(dfloat): Velocity along Z\n
      !!        mu(dfloat): Cosine of polar angle\n
      !!       mux(dfloat): Cosine of azimuth\n
      !!    muy_in(dfloat): Sign of sin of azimuth\n
      !!   MPID(MPI_class): Structure with MPI data\n
      !!    Stk(dfloat(:)): Intensity at bottom boundary
      subroutine bottomI(omega,T,vx,vy,vz,mu,mux,muy_in,MPID,Stk)

      ! I/O

      type(MPI_class), intent(in):: MPID
      double precision, intent(in):: T,mu,mux,muy_in,vx,vy,vz
      double precision, dimension(:), intent(in)::  omega
      double precision, dimension(MPID%if0(pid):MPID%if1(pid)), &
                        intent(out):: Stk

      ! Local

      integer:: ifreq

      double precision:: muy,ct,st,cc,sc,vfac


      ! Calculate Doppler shift factor
      if (dyn) then

        ct = mu
        st = sqrt(1d0 - ct*ct)
        cc = mux
        muy = muy_in
        if (abs(muy).lt.1d-8) muy = 1d0
        sc = muy*sqrt(1d0 - cc*cc)/abs(muy)

        vfac = 1d0 - vx*st*cc - vy*st*sc - vz*ct

      ! Static
      else

        vfac = 1d0

      end if

      ! Make the boundary planckian
      do ifreq=MPID%if0(pid),MPID%if1(pid)

        Stk(ifreq) = planck(omega(ifreq)*vfac,T)

      end do

      end subroutine bottomI

!#####################################################################
!#####################################################################
!#####################################################################

      end module boundary_mod
