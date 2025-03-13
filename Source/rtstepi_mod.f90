      !> Short characteristics intensity
      module rtstepi_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     19/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     19/12/2024:    V4.0.0 - Removed OpenMP directives (TdPA)
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
!  RTStepi
!    Perform one short-characteristics step in the radiation transfer
!  equation for the intensity with BESSER interpolation
!
!  RTTauI
!    Calculate the optical depth at a given point in a given direction
!  and store where tau becomes unity
!
!  RTContrI
!    Calculate the contribution function to the intensity at a given
!  point in a given direction
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use math_mod
      use parameters_mod , only: TINYI,TINYCI,TINYT
      use rtstepaux_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Perform one short-characteristics step in the radiation
      !! transfer equation for the intensity with BESSER
      !! interpolation\n
      !!        iz(integer): Height index\n
      !!       ith(integer): Output direction polar index\n
      !!       iph(integer): Output direction azimuth index\n
      !!     nfreq(integer): Frequency arrays size\n
      !!        ds1(double): Geometrical distance backward\n
      !!        ds2(double): Geometrical distance forward\n
      !!      KM(double(:)): Absorptivity previous point\n
      !!      SM(double(:)): Emissivity previous point\n
      !!      KO(double(:)): Absorptivity current point\n
      !!      SO(double(:)): Emissivity current point\n
      !!      KP(double(:)): Absorptivity next point\n
      !!      SP(double(:)): Emissivity next point\n
      !!    StkM(double(:)): Intensity previous point\n
      !!    StkO(double(:)): Intensity current point\n
      !!      LO(double(:)): Lambda operator current point\n
      !!       ALI(logical): If computing Lambda operator\n
      !!  nolinear(logical): If the interpolation can be non-linear
      subroutine RTStepI(iz,ith,iph,nfreq,ds1,ds2,KM,SM,KO,SO,KP,SP, &
                         StkM,StkO,LO,ALI,nolinear)

      ! I/O

      logical, intent(in):: ALI,nolinear
      integer, intent(in):: iz,ith,iph,nfreq
      double precision, intent(in):: ds1,ds2
      double precision, dimension(:), intent(in):: StkM
      double precision, dimension(:), intent(in):: KM,SM,KO,SO,KP,SP
      double precision, dimension(:), intent(out):: StkO,LO

      ! Local

      integer:: ifreq

      double precision:: exu,dt1,dt2,psim,psio,om_m,om_o,om_c,C0


      ! Routine name
      urou = 'RTStepI'

      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate optical distance between the short characteristic
        ! points
        call ftau(KM(ifreq),KO(ifreq),ds1,dt1)
        if (nolinear) call ftau(KO(ifreq),KP(ifreq),ds2,dt2)

        ! Get exponential
        exu = diexp(dt1)

        !
        ! Formal solution
        !

        ! Attenuation of incoming intensity
        StkO(ifreq) = StkM(ifreq)*exu

        ! If Besser with three points
        if (nolinear) then

          ! Besser interpolation of source function
          call RTomega(exu,dt1,om_m,om_o,om_c)
          call QBezierC0(dt1,dt2,SM(ifreq),SO(ifreq),SP(ifreq),C0)

          ! Stokes calculation
          StkO(ifreq) = StkO(ifreq) + &
                        om_m*SM(ifreq) + om_o*SO(ifreq) + om_c*C0

          ! Lambda operator
          if (ALI) LO(ifreq) = om_o + om_c

        ! If linear with two points
        else

          ! Linear interpolation of source function
          call psi_lin(exu,dt1,psim,psio)

          ! Stokes calculation
          StkO(ifreq) = StkO(ifreq) + psim*SM(ifreq) + psio*SO(ifreq)

          ! Lambda operator
          if (ALI) LO(ifreq) = psio

        end if


        !
        ! Control the physicallity of the result
        !

        ! If NaN
        if(isnan(StkO(ifreq)))then

          ! Issue error
          write(umsg,'(A,3(i4,","),i5,A,'// &
                     '2(A,1x,es9.2),'// &
                     '3(A,2(1x,es9.2)))') &
          'Error in RTStepI: NaN in intensity'// &
          new_line('A')//'(ith,iph,iz,ifreq)=(', &
          ith,iph,iz,ifreq,')', &
          new_line('A')//'Stokes(M):',StkM(ifreq), &
          new_line('A')//'Stokes(O):',StkO(ifreq), &
          new_line('A')//'eta(M),S(M):',KM(ifreq),SM(ifreq), &
          new_line('A')//'eta(O),S(O):',KO(ifreq),SO(ifreq), &
          new_line('A')//'eta(P),S(P):',KP(ifreq),SP(ifreq)

          call abortedS(umsg,urou,.True.,.True.)

        endif ! if NaN

        ! Control underflow and negativity
        if (StkO(ifreq).lt.TINYI) StkO(ifreq) = 0d0

      end do ! Frequencies

      return

      end subroutine RTStepI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the optical depth at a given point in a given
      !! direction and store where tau becomes unity\n
      !!        ds1(double): Geometrical distance backward\n
      !!     nfreq(integer): Frequency arrays size\n
      !!         z1(double): Height previous point\n
      !!         z2(double): Height current point\n
      !!    etaM(double(:)): Absorptivity for intensity previous
      !!                     point\n
      !!    etaO(double(:)): Absorptivity for intensity current
      !!                     point\n
      !!    tauM(double(:)): Optical depth previous point\n
      !!    tauO(double(:)): Optical depth current point\n
      !!  tau1(double(:,:)): Height where optical depth is equal
      !!                     to one
      subroutine RTTauI(ds1,nfreq,z1,z2,etaM,etaO,tauM,tauO,tau1)

      ! I/O

      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,z1,z2
      double precision, dimension(:), intent(in):: etaM,etaO,tauM
      double precision, dimension(:), intent(inout):: tauO
      double precision, dimension(:,:), intent(inout):: tau1

      ! Local

      integer:: ifreq

      double precision:: a,b,D,dt1


      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate optical distance between the two points
        call ftau(etaM(ifreq),etaO(ifreq),ds1,dt1)

        ! Add this to the optical depth
        tauO(ifreq) = tauM(ifreq) + dt1
        tau1(1,ifreq) = tauO(ifreq)

        ! If the optical was lower than one before
        if (tauM(ifreq).lt.1d0) then

          ! And now is larger than one, interpolate heights
          if (tauO(ifreq).ge.1d0) then

            ! If it is exactly one
            if (abs(tauO(ifreq) - 1d0).lt.TINYT) then

              ! We know the height
              tau1(2,ifreq) = z2

            ! If it is not exact
            else

              ! Interpolate the height
              D = tauM(ifreq) - tauO(ifreq)
              a = z1 - z2
              b = tauM(ifreq)*z2 - tauO(ifreq)*z1
              tau1(2,ifreq) = (a + b)/D

            end if ! exactly one

          ! If it is not larger than one yet
          else

            ! Update the variables
            tau1(2,ifreq) = z2

          end if ! New tau larger than one
        end if ! Previous tau larger than one

      end do ! Every frequency

      end subroutine RTTauI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the contribution function to the intensity at a
      !! given point in a given direction\n
      !!        ds1(double): Geometrical distance backward\n
      !!        ds2(double): Geometrical distance forward\n
      !!     nfreq(integer): Frequency arrays size\n
      !!        dz1(double): Geometrical vertical distance backward\n
      !!        dz2(double): Geometrical vertical distance forward\n
      !!     K0M(double(:)): Absorptivity previous point\n
      !!     K0O(double(:)): Absorptivity current point\n
      !!      SO(double(:)): Emissivity current point\n
      !!     K0P(double(:)): Absorptivity next point\n
      !!     tau(double(:)): Optical depth at current point\n
      !!   Contr(double(:)): Contribution function at current point\n
      !!  nolinear(logical): If the interpolation can be non-linear
      subroutine RTContrI(ds1,ds2,nfreq,dz1,dz2,K0M,K0O,SO,K0P, &
                          tau,Contr,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,ds2,dz1,dz2
      double precision, dimension(:), intent(in):: tau,K0M,K0O,SO,K0P
      double precision, dimension(:), intent(out):: Contr

      ! Local

      integer:: ifreq

      double precision:: dtau,exu,dt1,dt2


      ! Routine name
      urou = 'RTContrI'

      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate optical distance with the previous point
        call ftau(K0M(ifreq),K0O(ifreq),ds1,dt1)
        dtau = dt1/dz1

        ! Calculate optical distance with the next point and the
        ! optical depth derivative
        if(nolinear)then
          call ftau(K0O(ifreq),K0P(ifreq),ds2,dt2)
          dtau = .5d0*(dtau + dt2/dz2)
        end if

        ! Exponential
        exu = diexp(tau(ifreq))

        ! Compute contribution function
        Contr(ifreq) = SO(ifreq)*exu*dtau

        !
        ! Control the physicallity of the result
        !

        ! Control underflow
        if (abs(Contr(ifreq)).lt.TINYCI) Contr(ifreq) = 0d0

        ! Control NaN
        if (isnan(Contr(ifreq))) then

          ! Issue error
          write(umsg,'(2(A,1x,es9.2),'// &
                     '5(A,1x,es9.2))') &
          'Error in RTContrI: NaN in Contribution Function'// &
          new_line('A')//'Contribution:',Contr(ifreq), &
          new_line('A')//'eta(M):',K0M(ifreq), &
          new_line('A')//'eta(O):',K0O(ifreq), &
          new_line('A')//'S(O):',SO(ifreq), &
          new_line('A')//'eta(P):',K0P(ifreq)

          call abortedS(umsg,urou,.True.,.True.)

        end if ! If NaN

      end do ! All frequencies

      return

      end subroutine RTContrI

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstepi_mod
