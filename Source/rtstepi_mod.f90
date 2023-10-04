      !> Short characteristics intensity
      module rtstepi_mod
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
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS and generation
!                             of the error message (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Moved the frequency loop inside and
!                             added OpenMP (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Now uses diexp (TdPA)
!                           - Now uses especific TINY variables (TdPA)
!
!     07/11/2017:    V1.0.4 - Bugfix:One of the conditional
!                             branchings in RTTauI was wrong (TdPA)
!
!     06/09/2017:    V1.0.3 - Now the negative intensities are forced
!                             to 0, instead of aborting the whole
!                             calculation (TdPA)
!                           - The underflow is checked together with
!                             the sign (TdPA)
!
!     05/05/2017:    V1.0.2 - RTStep knows which point is computing
!                             during the call (TdPA)
!
!     05/01/2017:    V1.0.1 - RTTauI now stores the current tau in
!                             tau1(1) too (TdPA)
!
!     04/20/2017:    V1.0.0 - Started module (TdPA)
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
!  RTStepi:
!    Applies short characteristics to the points in input (only
!    intensity)
!
!  RTTauI:
!    Calculates the tau and stores where tau=1 (only intensity)
!
!  RTContrI:
!    Applies short characteristics to the points in input to calculate
!    the contribution function (only intensity)
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use math_mod
      use omp_mod
      use parameters_mod , only: TINYI,TINYCI,TINYT
      use rtstepaux_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Advances the radiation transfer equation propagation for
      !! intensity using short characteristics with the BESSER
      !! interpolation algorithm.\n
      !!           iz(integer): Height index\n
      !!          ith(integer): Output direction polar index\n
      !!          iph(integer): Output direction azimuth index\n
      !!        nfreq(integer): Frequency size\n
      !!           ds1(dfloat): Geometrical distance backward\n
      !!           ds2(dfloat): Geometrical distance forward\n
      !!         KM(dfloat(:)): Absorptivity previous point\n
      !!         SM(dfloat(:)): Emissivity previous point\n
      !!         KO(dfloat(:)): Absorptivity current point\n
      !!         SO(dfloat(:)): Emissivity current point\n
      !!         KP(dfloat(:)): Absorptivity next point\n
      !!         SP(dfloat(:)): Emissivity next point\n
      !!       StkM(dfloat(:)): Intensity previous point\n
      !!       StkO(dfloat(:)): Intensity current point\n
      !!         LO(dfloat(:)): Lambda operator current point\n
      !!          ALI(logical): If computing Lambda operator\n
      !!     nolinear(logical): Bool to specify that the interpolation
      !!                        must not be linear
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

      integer:: ifreq,tid=-1

      double precision:: exu,dt1,dt2
      double precision:: psim,psio
      double precision:: om_m,om_o,om_c,C0

      ! Routine name
      urou = 'RTStepI'

!$omp parallel do default(none) &
!$omp private(ifreq,exu,dt1,dt2,psim,psio,om_m,om_o,om_c,C0,tid) &
!$omp private(umsg) &
!$omp shared(ALI,nolinear,iz,ith,iph,nfreq,ds1,ds2,laborted) &
!$omp shared(StkM,KM,SM,KO,SO,KP,SP,StkO,LO,urou)

      ! For each frequency
      do ifreq=1,nfreq

        !
        ! Calculate optical distance between the short characteristic
        ! points
        !
        call ftau(KM(ifreq),KO(ifreq),ds1,dt1)
        if(nolinear)call ftau(KO(ifreq),KP(ifreq),ds2,dt2)

        exu = diexp(dt1)

        ! Attenuation of incoming intensity
        StkO(ifreq) = StkM(ifreq)*exu


        !
        ! Formal solution
        !

        ! If Besser with three points
        if(nolinear)then

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

        ! Control NaN
        if(isnan(StkO(ifreq)))then
#ifdef _OPENMP
          tid = omp_get_thread_num() + 1
#endif
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

          call abortedS(umsg,urou,tid,.True.,.True.)

        endif ! if NaN

        ! Control underflow and negativity
        if(StkO(ifreq).lt.TINYI) StkO(ifreq) = 0d0

      end do
!$omp end parallel do

      return

      end subroutine RTStepI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds to the optical depth the contribution of a section
      !! between two points and checks if the optical depth becomes
      !! one.\n
      !!           ds1(dfloat): Geometrical distance backward\n
      !!        nfreq(integer): Frequency size\n
      !!            z1(dfloat): Height previous point\n
      !!            z2(dfloat): Height current point\n
      !!       etaM(dfloat(:)): Absorptivity for intensity previous
      !!                        point\n
      !!       etaO(dfloat(:)): Absorptivity for intensity current
      !!                        point\n
      !!       tauM(dfloat(:)): Optical depth previous point\n
      !!       tauO(dfloat(:)): Optical depth current point\n
      !!     tau1(dfloat(:,:)): Height where optical depth is equal
      !!                        to one
      subroutine RTTauI(ds1,nfreq,z1,z2,etaM,etaO,tauM,tauO,tau1)

      ! I/O
      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,z1,z2
      double precision, dimension(:), intent(in):: etaM,etaO
      double precision, dimension(:), intent(inout):: tauM,tauO
      double precision, dimension(:,:), intent(inout):: tau1

      ! Local

      integer:: ifreq

      double precision:: a,b,D,dt1


!$omp parallel do default(none) &
!$omp private(ifreq,a,b,D,dt1) &
!$omp shared(ds1,nfreq,z1,z2,etaM,etaO,tauM,tauO,tau1)
      do ifreq=1,nfreq

        !
        ! Calculate optical distance between the two points
        !
        call ftau(etaM(ifreq),etaO(ifreq),ds1,dt1)

        ! Add this to the optical depth
        tauO(ifreq) = tauM(ifreq) + dt1
        tau1(1,ifreq) = tauO(ifreq)

        ! If the optical was lower than one before
        if (tauM(ifreq).lt.1d0) then

          ! And now is larger than one, interpolate heights
          if (tauO(ifreq).ge.1d0) then

            ! If it is exactly one, we know the height
            if (abs(tauO(ifreq) - 1d0).lt.TINYT) then

              tau1(2,ifreq) = z2

            ! If it is not exact, interpolate the height
            else

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

!$omp end parallel do

      end subroutine RTTauI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the contribution function to the intensity in a step
      !! of the short characteristics.\n
      !!           ds1(dfloat): Geometrical distance backward\n
      !!           ds2(dfloat): Geometrical distance forward\n
      !!        nfreq(integer): Frequency size\n
      !!           dz1(dfloat): Geometrical vertical distance
      !!                        backward\n
      !!           dz2(dfloat): Geometrical vertical distance
      !!                        forward\n
      !!        K0M(dfloat(:)): Absorptivity previous point\n
      !!        K0O(dfloat(:)): Absorptivity current point\n
      !!         SO(dfloat(:)): Emissivity current point\n
      !!        K0P(dfloat(:)): Absorptivity next point\n
      !!        tau(dfloat(:)): Optical depth at current point\n
      !!      Contr(dfloat(:)): Contribution function at current
      !!                        point\n
      !!     nolinear(logical): Bool to specify that the interpolation
      !!                        must not be linear
      subroutine RTContrI(ds1,ds2,nfreq,dz1,dz2,K0M,K0O,SO,K0P, &
                          tau,Contr,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,ds2,dz1,dz2
      double precision, dimension(:), intent(in):: tau,K0M,K0O,SO,K0P
      double precision, dimension(:), intent(out):: Contr

      ! Local

      integer:: ifreq,tid=-1
      double precision:: dtau,exu,dt1,dt2


      ! Routine name
      urou = 'RTContrI'

!$omp parallel do default(none) &
!$omp private(ifreq,exu,dt1,dt2,dtau,tid,umsg) &
!$omp shared(nfreq,ds1,ds2,dz1,dz2,laborted) &
!$omp shared(K0M,K0O,SO,K0P,tau,Contr,nolinear,urou)
      do ifreq=1,nfreq

        !
        ! Calculate optical distance with the previous point
        !
        call ftau(K0M(ifreq),K0O(ifreq),ds1,dt1)
        dtau = dt1/dz1

        !
        ! Calculate optical distance with the next point and the
        ! optical depth derivative
        !
        if(nolinear)then
          call ftau(K0O(ifreq),K0P(ifreq),ds2,dt2)
          dtau = .5d0*(dtau + dt2/dz2)
        end if

        exu = diexp(tau(ifreq))

        !
        ! Compute contribution function
        !
        Contr(ifreq) = SO(ifreq)*exu*dtau


        !
        ! Control the physicallity of the result
        !

        ! Control underflow
        if (abs(Contr(ifreq)).lt.TINYCI) Contr(ifreq)=.0D0

        ! Control NaN
        if (isnan(Contr(ifreq))) then
#ifdef _OPENMP
          tid = omp_get_thread_num() + 1
#endif
          write(umsg,'(2(A,1x,es9.2),'// &
                     '5(A,1x,es9.2))') &
          'Error in RTContrI: NaN in Contribution Function'// &
          new_line('A')//'Contribution:',Contr(ifreq), &
          new_line('A')//'eta(M):',K0M(ifreq), &
          new_line('A')//'eta(O):',K0O(ifreq), &
          new_line('A')//'S(O):',SO(ifreq), &
          new_line('A')//'eta(P):',K0P(ifreq)

          call abortedS(umsg,urou,tid,.True.,.True.)

        end if ! If NaN

      end do ! All frequencies
!$omp end parallel do

      return

      end subroutine RTContrI

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstepi_mod
