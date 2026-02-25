      !> Short characteristics algorithm
      module rtstep_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     25/02/2026 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     25/02/2026:    V4.0.2 - Added the length of the path along the
!                             los to the debugging information when
!                             the output Stokes parameters are non
!                             physical (TdPA)
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
!  RTStep
!    Perform one short-characteristics step in the radiation transfer
!  equation with BESSER interpolation
!
!  RTContr
!    Calculate the contribution function at a given point in a given
!  direction
!
!  RTCStep:
!    Solve the radiative transfer equation with short-characteristics
!  in a constant properties slab
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use math_mod
      use parameters_mod , only: vacuum , TINYS , TINYT , TINYCS
      use rtstepaux_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Perform one short-characteristics step in the radiation
      !! transfer equation with BESSER interpolation\n
      !!        iz(integer): Height index\n
      !!       ith(integer): Output direction polar index\n
      !!       iph(integer): Output direction azimuth index\n
      !!     nfreq(integer): Number of frequencies\n
      !!        ds1(double): Geometrical distance backward\n
      !!        ds2(double): Geometrical distance forward\n
      !!   K0M(double(:,:)): First row absorption matrix previous
      !!                     point\n
      !!   K1M(double(:,:)): Second row absorption matrix previous
      !!                     point\n
      !!   K2M(double(:,:)): Third row absorption matrix previous
      !!                     point\n
      !!    SM(double(:,:)): Emissivity previous point\n
      !!   K0O(double(:,:)): First row absorption matrix current
      !!                     point\n
      !!   K1O(double(:,:)): Second row absorption matrix current
      !!                     point\n
      !!   K2O(double(:,:)): Third row absorption matrix current
      !!                     point\n
      !!    SO(double(:,:)): Emissivity current point\n
      !!   K0P(double(:,:)): First row absorption matrix next
      !!                     point\n
      !!    SP(double(:,:)): Emissivity next point\n
      !!  StkM(double(:,:)): Stokes parameters previous point\n
      !!  StkO(double(:,:)): Stokes parameters current point\n
      !!  nolinear(logical): If the interpolation can be non-linear
      subroutine RTStep(iz,ith,iph,nfreq,ds1,ds2,K0M,K1M,K2M,SM, &
                        K0O,K1O,K2O,SO,K0P,SP,StkM,StkO,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: iz,nfreq,ith,iph
      double precision, intent(in):: ds1,ds2
      double precision, dimension(0:3,nfreq), intent(in):: StkM
      double precision, dimension(0:3,nfreq), intent(in):: K0M,K1M
      double precision, dimension(0:3,nfreq), intent(in):: K2M,SM
      double precision, dimension(0:3,nfreq), intent(in):: K0O,K1O
      double precision, dimension(0:3,nfreq), intent(in):: K2O,SO
      double precision, dimension(0:3,nfreq), intent(in):: K0P,SP
      double precision, dimension(0:3,nfreq), intent(out):: StkO

      ! Local

      integer:: ifreq,iS,i1

      double precision:: exu,dt1,dt2,psim,psio,om_m,om_o,om_c
      double precision, dimension(0:1):: CA,vector1A,vector2A,vector3A
      double precision, dimension(0:3):: C0,vector1,vector2,vector3
      double precision, dimension(0:1,0:1):: kappaA,matrixA
      double precision, dimension(0:3,0:3):: kappa,matrix


      ! Routine name
      urou = 'RTStep'

      ! Initialize Stokes to calculate
      StkO = 0d0

      !
      ! Choose the method of solution given the symmetry
      !

      !
      ! If the problem is axial (only I and Q)
      !
      if (RTaxial) then !I+Q/I-Q

        ! For each frequency
        do ifreq=1,nfreq

          ! Calculate optical distance
          call ftau(K0M(0,ifreq),K0O(0,ifreq),ds1,dt1)
          if(nolinear)call ftau(K0O(0,ifreq),K0P(0,ifreq),ds2,dt2)

          ! Exponential
          exu = diexp(dt1)

          ! Calculate linear coefficients
          call psi_lin(exu,dt1,psim,psio)

          ! If there is a P point
          if (nolinear) then

            ! Calculate BESSER coefficients
            call RTomega(exu,dt1,om_m,om_o,om_c)

            ! For each Stokes parameter
            do iS=0,1

              ! Calculate correction
              call QBezierC0(dt1,dt2,SM(iS,ifreq),SO(iS,ifreq), &
                                     SP(iS,ifreq),CA(iS))

            end do ! Stokes parameter

          endif ! There is P point

          !
          ! Compute the matrices of the RT solution
          !

          ! First row and column
          kappaA(1,0) = psio*K0O(1,ifreq)
          matrixA(1,0) = -psim*K0M(1,ifreq)
          kappaA(0,1) = kappaA(1,0)
          matrixA(0,1) = matrixA(1,0)

          ! Diagonal elements
          do i1=0,1
            kappaA(i1,i1) = 1d0
            matrixA(i1,i1) = exu
          end do

          ! Invert kappa
          call Matinv_2d(kappaA)

          ! Compute the Stokes solution
          call MatVec_2d(matrixA,StkM(0:1,ifreq),vector1A)
          call MatVec_2d(kappaA,vector1A,vector2A)

          ! If three points short characteristics
          if (nolinear) then

            ! Compute source function contribution
            call MatVec_2d(kappaA,om_m*SM(0:1,ifreq) + &
                                  om_o*SO(0:1,ifreq) + &
                                  om_c*CA,vector3A)

          ! If two points short characteristics
          else

            ! Compute source function contribution
            call MatVec_2d(kappaA,psim*SM(0:1,ifreq) + &
                                  psio*SO(0:1,ifreq), &
                           vector3A)

          end if ! Three or two points

          ! Compute final Stokes parameters
          StkO(0:1,ifreq) = vector2A + vector3A


          !
          ! Control the physicallity of the result
          !

          ! For each Stokes parameter
          do iS=0,1

            ! Control underflow
            if (abs(StkO(iS,ifreq)).lt.TINYS) StkO(iS,ifreq) = 0d0

            ! Control NaN
            if (ieee_is_nan(StkO(iS,ifreq))) then

              ! Issue error
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es9.2)),'// &
                         '6(A,2(1x,es9.2)))') &
              'Error in RTStep: NaN in Stokes'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'ds(M):',ds1, &
              new_line('A')//'ds(P):',ds2, &
              new_line('A')//'Stk(M):',StkM(:,ifreq), &
              new_line('A')//'Stk(O):',StkO(:,ifreq), &
              new_line('A')//'etaI(M),SI(M):', &
                          K0M(0,ifreq),SM(0,ifreq), &
              new_line('A')//'etaI(O),SI(O):', &
                          K0O(0,ifreq),SO(0,ifreq), &
              new_line('A')//'etaI(P),SI(P):', &
                          K0P(0,ifreq),SP(0,ifreq), &
              new_line('A')//'etaQ(M),SQ(M):', &
                          K0M(1,ifreq),SM(1,ifreq), &
              new_line('A')//'etaQ(O),SQ(O):', &
                          K0O(1,ifreq),SO(1,ifreq), &
              new_line('A')//'etaQ(P),SQ(P):', &
                          K0P(1,ifreq),SP(1,ifreq)
              call abortedS(umsg,urou,.True.,.True.)

            endif ! if NaN

          end do ! Stokes parameters

          ! Control Physics
          if (StkO(0,ifreq).lt.0d0.or. &
              StkO(0,ifreq)*StkO(0,ifreq).lt. &
              (StkO(1,ifreq)*StkO(1,ifreq) + &
               StkO(2,ifreq)*StkO(2,ifreq) + &
               StkO(3,ifreq)*StkO(3,ifreq))) then

            ! Write message
            write(umsg,'(A,3(i4,","),i5,A,'// &
                       '2(A,1x,es8.1),'// &
                       '2(A,4(1x,es9.2)),'// &
                       '6(A,2(1x,es9.2)),A)') &
            'Non-phys. Stokes'// &
            new_line('A')//'(ith,iph,iz,ifreq)=(', &
            ith,iph,iz,ifreq,')', &
            new_line('A')//'ds(M):',ds1, &
            new_line('A')//'ds(P):',ds2, &
            new_line('A')//'Stk(M):',StkM(:,ifreq), &
            new_line('A')//'Stk(O):',StkO(:,ifreq), &
            new_line('A')//'etaI(M),SI(M):', &
                        K0M(0,ifreq),SM(0,ifreq), &
            new_line('A')//'etaI(O),SI(O):', &
                        K0O(0,ifreq),SO(0,ifreq), &
            new_line('A')//'etaI(P),SI(P):', &
                        K0P(0,ifreq),SP(0,ifreq), &
            new_line('A')//'etaQ(M),SQ(M):', &
                        K0M(1,ifreq),SM(1,ifreq), &
            new_line('A')//'etaQ(O),SQ(O):', &
                        K0O(1,ifreq),SO(1,ifreq), &
            new_line('A')//'etaQ(P),SQ(P):', &
                        K0P(1,ifreq),SP(1,ifreq), &
            new_line('A')//' Do lin.'

            ! If no linear
            if (nolinear) then

              ! Issue warning
              call abortedS(umsg,urou,.False.,.True.)

            ! If already linear
            else

              ! Abort
              call abortedS(umsg,urou,.not.nphysS,.True.)

            end if ! Linear or not

            ! Try linear
            call MatVec_2d(kappaA,psim*SM(0:1,ifreq) + &
                                  psio*SO(0:1,ifreq), &
                           vector3A)

            ! Compute final Stokes parameters
            StkO(0:1,ifreq) = vector2A + vector3A

            ! Control Physics linear
            if (StkO(0,ifreq).lt.0d0.or. &
                StkO(0,ifreq)*StkO(0,ifreq).lt. &
                (StkO(1,ifreq)*StkO(1,ifreq) + &
                 StkO(2,ifreq)*StkO(2,ifreq) + &
                 StkO(3,ifreq)*StkO(3,ifreq))) then

              ! Issue error
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es9.2)),'// &
                         '6(A,2(1x,es9.2)))') &
              'Non-phys. Stokes'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'ds(M):',ds1, &
              new_line('A')//'ds(P):',ds2, &
              new_line('A')//'Stk(M):',StkM(:,ifreq), &
              new_line('A')//'Stk(O):',StkO(:,ifreq), &
              new_line('A')//'etaI(M),SI(M):', &
                          K0M(0,ifreq),SM(0,ifreq), &
              new_line('A')//'etaI(O),SI(O):', &
                          K0O(0,ifreq),SO(0,ifreq), &
              new_line('A')//'etaI(P),SI(P):', &
                          K0P(0,ifreq),SP(0,ifreq), &
              new_line('A')//'etaQ(M),SQ(M):', &
                          K0M(1,ifreq),SM(1,ifreq), &
              new_line('A')//'etaQ(O),SQ(O):', &
                          K0O(1,ifreq),SO(1,ifreq), &
              new_line('A')//'etaQ(P),SQ(P):', &
                          K0P(1,ifreq),SP(1,ifreq)

              call abortedS(umsg,urou,.not.nphysS,.True.)

            end if ! Physical Stokes parameters
          end if ! Physical Stokes parameters

        end do ! For each frequency

      !
      ! Full Stokes problem
      !
      else

        ! For each frequency
        do ifreq=1,nfreq

          ! Calculate optical distance
          call ftau(K0M(0,ifreq),K0O(0,ifreq),ds1,dt1)
          if(nolinear)call ftau(K0O(0,ifreq),K0P(0,ifreq),ds2,dt2)

          ! Exponential
          exu = diexp(dt1)

          ! Calculate linear coefficients
          call psi_lin(exu,dt1,psim,psio)

          ! Calculate BESSER coefficients
          if(nolinear)then

            ! Calculate BESSER coefficients
            call RTomega(exu,dt1,om_m,om_o,om_c)

            ! For each Stokes parameter
            do iS=0,3

              ! Calculate correction
              call QBezierC0(dt1,dt2,SM(iS,ifreq),SO(iS,ifreq), &
                                     SP(iS,ifreq),C0(iS))

            end do ! Stokes parameter

          endif ! There is P point

          !
          ! Compute the matrices of the RT solution
          !

          ! First row and column
          kappa(1:3,0) = psio*K0O(1:3,ifreq)
          matrix(1:3,0) = -psim*K0M(1:3,ifreq)
          kappa(0,1:3) = kappa(1:3,0)
          matrix(0,1:3) = matrix(1:3,0)

          ! Second row
          kappa(1,2) = K1O(2,ifreq)*psio
          kappa(1,3) = K1O(3,ifreq)*psio
          matrix(1,2) = -psim*K1M(2,ifreq)
          matrix(1,3) = -psim*K1M(3,ifreq)

          ! Third row
          kappa(2,1) = -kappa(1,2)
          kappa(2,3) = K2O(3,ifreq)*psio
          matrix(2,1) = -matrix(1,2)
          matrix(2,3) = -psim*K2M(3,ifreq)

          ! Fourth row
          kappa(3,1) = -kappa(1,3)
          kappa(3,2) = -kappa(2,3)
          matrix(3,1) = -matrix(1,3)
          matrix(3,2) = -matrix(2,3)

          ! Diagonal elements
          do i1=0,3
            kappa(i1,i1) = 1d0
            matrix(i1,i1) = exu
          end do

          ! Invert kappa
          call Matinv(kappa)

          ! Compute the Stokes solution
          call MatVec(matrix,StkM(:,ifreq),vector1)
          call MatVec(kappa,vector1,vector2)

          ! If three points short characteristics
          if (nolinear) then

            ! Compute source function contribution
            call MatVec(kappa,om_m*SM(:,ifreq) + om_o*SO(:,ifreq) + &
                              om_c*C0,vector3)

          ! If two points short characteristics
          else

            ! Compute source function contribution
            call MatVec(kappa,psim*SM(:,ifreq) + psio*SO(:,ifreq), &
                        vector3)

          end if ! Three or two points

          ! Compute final Stokes parameters
          StkO(:,ifreq) = vector2 + vector3

          !
          ! Control the physicallity of the result
          !

          ! For each Stokes parameter
          do iS=0,3

            ! Control underflow
            if (abs(StkO(iS,ifreq)).lt.TINYS) StkO(iS,ifreq) = 0d0

            ! Control NaN
            if (ieee_is_nan(StkO(iS,ifreq))) then

              ! Issue error
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es8.1)),'// &
                         '3(A,5(1x,es8.1)))') &
              'NaN in Stokes'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'ds(M):',ds1, &
              new_line('A')//'ds(P):',ds2, &
              new_line('A')//'Stk(M):',StkM(:,ifreq), &
              new_line('A')//'Stk(O):',StkO(:,ifreq), &
              new_line('A')//'etaI(M),S(M):', &
                          K0M(0,ifreq),SM(:,ifreq), &
              new_line('A')//'etaI(O),S(O):', &
                          K0O(0,ifreq),SO(:,ifreq), &
              new_line('A')//'etaI(P),S(P):', &
                          K0P(0,ifreq),SP(:,ifreq)
              call abortedS(umsg,urou,.True.,.True.)
              write(umsg,'(2(A,3(4(1x,es8.1),A)),'// &
                         'A,4(1x,es8.1))') &
              'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(:,ifreq),new_line('A'), &
                           K2M(:,ifreq),new_line('A'), &
              'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(:,ifreq),new_line('A'), &
                           K2O(:,ifreq),new_line('A'), &
              'Kabs(P)'//new_line('A'), &
                           K0P(:,ifreq)
              call abortedS(umsg,urou,.True.,.True.)

            endif ! if NaN

          end do ! Stokes parameters

          ! Control Physics
          if (StkO(0,ifreq).lt.0d0.or. &
              StkO(0,ifreq)*StkO(0,ifreq).lt. &
              (StkO(1,ifreq)*StkO(1,ifreq) + &
               StkO(2,ifreq)*StkO(2,ifreq) + &
               StkO(3,ifreq)*StkO(3,ifreq))) then

            ! If was doing non-linear
            if (nolinear) then

              ! Issue warning
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es8.1)),'// &
                         '3(A,5(1x,es8.1)))') &
                'Non-phys. Stokes'// &
                new_line('A')//'(ith,iph,iz,ifreq)=(', &
                ith,iph,iz,ifreq,')', &
                new_line('A')//'ds(M):',ds1, &
                new_line('A')//'ds(P):',ds2, &
                new_line('A')//'Stk(M):',StkM(:,ifreq), &
                new_line('A')//'Stk(O):',StkO(:,ifreq), &
                new_line('A')//'etaI(M),S(M):', &
                            K0M(0,ifreq),SM(:,ifreq), &
                new_line('A')//'etaI(O),S(O):', &
                            K0O(0,ifreq),SO(:,ifreq), &
                new_line('A')//'etaI(P),S(P):', &
                            K0P(0,ifreq),SP(:,ifreq)
              call abortedS(umsg,urou,.False.,.True.)
              write(umsg,'(2(A,4(1x,es8.1),A,'// &
                            '18x,2(1x,es8.1),A,'// &
                            '28x,es8.1,A),A,4(1x,es8.1),A)') &
                new_line('A')//'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(2:3,ifreq),new_line('A'), &
                           K2M(3,ifreq),new_line('A'), &
                'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(2:3,ifreq),new_line('A'), &
                           K2O(3,ifreq),new_line('A'), &
                'Kabs(P)'//new_line('A'), &
                             K0P(:,ifreq), &
                new_line('A')//' do lin.'
              call abortedS(umsg,urou,.False.,.True.)

            ! If already linear
            else

              ! Issue error
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es8.1)),'// &
                         '2(A,5(1x,es8.1)))') &
                'Non-phys. Stokes'// &
                new_line('A')//'(ith,iph,iz,ifreq)=(', &
                ith,iph,iz,ifreq,')', &
                new_line('A')//'ds(M):',ds1, &
                new_line('A')//'ds(P):',ds2, &
                new_line('A')//'Stk(M):',StkM(:,ifreq), &
                new_line('A')//'Stk(O):',StkO(:,ifreq), &
                new_line('A')//'etaI(M),S(M):', &
                            K0M(0,ifreq),SM(:,ifreq), &
                new_line('A')//'etaI(O),S(O):', &
                            K0O(0,ifreq),SO(:,ifreq)
              call abortedS(umsg,urou,.False.,.True.)
              write(umsg,'(2(A,4(1x,es8.1),A,'// &
                            '18x,2(1x,es8.1),A,'// &
                            '28x,es8.1,A))') &
                new_line('A')//'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(2:3,ifreq),new_line('A'), &
                           K2M(3,ifreq),new_line('A'), &
                'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(2:3,ifreq),new_line('A'), &
                           K2O(3,ifreq),''
              call abortedS(umsg,urou,.not.nphysS,.True.)

            end if ! Non-linear

            ! Try linear
            call MatVec(kappa,psim*SM(:,ifreq) + psio*SO(:,ifreq), &
                        vector3)

            ! Compute final Stokes parameters
            StkO(:,ifreq) = vector2 + vector3

            ! Control Physics linear
            if (StkO(0,ifreq).lt.0d0.or. &
                StkO(0,ifreq)*StkO(0,ifreq).lt. &
                (StkO(1,ifreq)*StkO(1,ifreq) + &
                 StkO(2,ifreq)*StkO(2,ifreq) + &
                 StkO(3,ifreq)*StkO(3,ifreq))) then

              ! Issue error
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,1x,es8.1),'// &
                         '2(A,4(1x,es8.1)),'// &
                         '2(A,5(1x,es8.1)))') &
                'Non-phys. Stokes'// &
                new_line('A')//'(ith,iph,iz,ifreq)=(', &
                ith,iph,iz,ifreq,')', &
                new_line('A')//'ds(M):',ds1, &
                new_line('A')//'ds(P):',ds2, &
                new_line('A')//'Stk(M):',StkM(:,ifreq), &
                new_line('A')//'Stk(O):',StkO(:,ifreq), &
                new_line('A')//'etaI(M),S(M):', &
                            K0M(0,ifreq),SM(:,ifreq), &
                new_line('A')//'etaI(O),S(O):', &
                            K0O(0,ifreq),SO(:,ifreq)
              call abortedS(umsg,urou,.False.,.True.)
              write(umsg,'(2(A,4(1x,es8.1),A,'// &
                            '18x,2(1x,es8.1),A,'// &
                            '28x,es8.1,A))') &
                new_line('A')//'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(2:3,ifreq),new_line('A'), &
                           K2M(3,ifreq),new_line('A'), &
                'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(2:3,ifreq),new_line('A'), &
                           K2O(3,ifreq),''
              call abortedS(umsg,urou,.not.nphysS,.True.)

            end if ! Physical Stokes parameters
          end if ! Physical Stokes parameters

        end do ! For each frequency

      end if ! Method of solution

      return

      end subroutine RTStep

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the contribution function at a given point in a
      !! given direction\n
      !!      nfreq(integer): Number of frequencies\n
      !!         ds1(double): Geometrical distance backward\n
      !!         ds2(double): Geometrical distance forward\n
      !!         dz1(double): Geometrical vertical distance
      !!                      backward\n
      !!         dz2(double): Geometrical vertical distance
      !!                      forward\n
      !!    K0M(double(:,:)): First row absorption matrix previous
      !!                      point\n
      !!    K0O(double(:,:)): First row absorption matrix current
      !!                      point\n
      !!    K1O(double(:,:)): Second row absorption matrix current
      !!                      point\n
      !!    K2O(double(:,:)): Third row absorption matrix current
      !!                      point\n
      !!     SO(double(:,:)): Emissivity current point\n
      !!    K0P(double(:,:)): First row absorption matrix next point\n
      !!   StkO(double(:,:)): Stokes parameters current point\n
      !!      tau(double(:)): Optical depth at current point\n
      !!  Contr(double(:,:)): Contribution function at current point\n
      !!   nolinear(logical): If the interpolation can be non-linear
      subroutine RTContr(nfreq,ds1,ds2,dz1,dz2,K0M, &
                         K0O,K1O,K2O,SO,K0P, &
                         StkO,tau,Contr,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,ds2,dz1,dz2
      double precision, dimension(:), intent(in):: tau
      double precision, dimension(0:3,nfreq), intent(in):: StkO,K0M
      double precision, dimension(0:3,nfreq), intent(in):: K0O,K1O
      double precision, dimension(0:3,nfreq), intent(in):: K2O,SO
      double precision, dimension(0:3,nfreq), intent(in):: K0P
      double precision, dimension(0:3,nfreq), intent(out):: Contr

      ! Local

      integer:: iS,ifreq

      double precision:: dtau,exu,dt1,dt2
      double precision, dimension(0:3):: vector1
      double precision, dimension(0:3,0:3):: kappa


      ! Routine name
      urou = 'RTContr'

      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate optical distance with the previous point
        call ftau(K0M(0,ifreq),K0O(0,ifreq),ds1,dt1)
        dtau = dt1/dz1

        ! Calculate optical distance with the next point and the
        ! optical depth derivative
        if(nolinear)then
          call ftau(K0O(0,ifreq),K0P(0,ifreq),ds2,dt2)
          dtau = .5d0*(dtau + dt2/dz2)
        end if

        ! Exponential
        exu = diexp(tau(ifreq))

        ! If axial
        if (RTaxial) then

          !
          ! Compute the polarization matrix part of the effective
          ! source function
          !
          vector1(0) = K0O(1,ifreq)*StkO(1,ifreq)
          vector1(1) = K0O(1,ifreq)*StkO(0,ifreq)


          !
          ! Compute contribution function
          !
          Contr(0,ifreq) = (SO(0,ifreq) - vector1(0))*exu*dtau
          Contr(1,ifreq) = (SO(1,ifreq) - vector1(1))*exu*dtau
          Contr(2,ifreq) = 0d0
          Contr(3,ifreq) = 0d0

        ! If not axial
        else

          !
          ! Compute the matrices
          !
          kappa(0,0) = 0d0
          kappa(0,1:3) = K0O(1:3,ifreq)
          kappa(1:3,0) = K0O(1:3,ifreq)
          kappa(1,1) = 0d0
          kappa(1,2:3) = K1O(2:3,ifreq)
          kappa(2,1) = -kappa(1,2)
          kappa(2,2) = 0d0
          kappa(2,3) = K2O(3,ifreq)
          kappa(3,1) = -kappa(1,3)
          kappa(3,2) = -kappa(2,3)
          kappa(3,3) = 0d0

          !
          ! Compute the polarization matrix part of the effective
          ! source function
          !
          call MatVec(kappa,StkO(:,ifreq),vector1)


          !
          ! Compute contribution function
          !
          Contr(:,ifreq) = (SO(:,ifreq) - vector1)*exu*dtau

        end if

        !
        ! Control the physicallity of the result
        !

        ! For each Stokes parameter
        do iS=0,3

          ! Control underflow
          if (abs(Contr(iS,ifreq)).lt.TINYCS) Contr(iS,ifreq)=.0D0

          ! Control NaN
          if (ieee_is_nan(Contr(iS,ifreq))) then

            ! Issue error
            write(umsg,'(2(A,1x,es9.2),'// &
                       'A,2(1x,es9.2),'// &
                       'A,1x,es9.2)') &
            'Error in RTContr: NaN in Contribution Function'// &
            new_line('A')//'ContrI:',Contr(0,ifreq), &
            new_line('A')//'etaI(M):',K0M(0,ifreq), &
            new_line('A')//'etaI(O),SI(O):', &
                        K0O(0,ifreq),SO(0,ifreq), &
            new_line('A')//'etaI(P):',K0P(0,ifreq)

            call abortedS(umsg,urou,.True.,.True.)

          end if ! If NaN

        end do ! Stokes parameters
      end do ! For each frequency

      end subroutine RTContr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the radiative transfer equation with
      !! short-characteristics in a constant properties slab\n
      !!    nfreq(integer): Number of frequencies\n
      !!      fact(double): Optical depth scale factor\n
      !!   K0(double(:,:)): First row absorption matrix\n
      !!   K1(double(:,:)): Second row absorption matrix\n
      !!   K2(double(:,:)): Third row absorption matrix\n
      !!    S(double(:,:)): Emissivity current point\n
      !!  Stk(double(:,:)): Stokes parameters
      subroutine RTCStep(nfreq,fact,K0,K1,K2,S,Stk)

      ! I/O

      integer, intent(in):: nfreq
      double precision, intent(in):: fact
      double precision, dimension(0:3,nfreq), intent(in):: K0,K1
      double precision, dimension(0:3,nfreq), intent(in):: K2,S
      double precision, dimension(0:3,nfreq), intent(inout):: Stk

      ! Local

      integer:: ifreq,iS

      double precision:: dt,exu,psim,psio
      double precision, dimension(0:3):: vector1,vector2,vector3
      double precision, dimension(0:3,0:3):: kappa,matrix


      ! Routine name
      urou = 'RTCStep'

      ! For each frequency
      do ifreq=1,nfreq

        ! Optical depth
        dt = K0(0,ifreq)*fact

        ! Exponential
        exu = diexp(dt)

        ! Calculate linear coefficients
        call psi_lin(exu,dt,psim,psio)

        !
        ! Compute the matrices of the RT solution
        !

        ! First row and column
        kappa(1:3,0) = psio*K0(1:3,ifreq)
        matrix(1:3,0) = -psim*K0(1:3,ifreq)
        kappa(0,1:3) = kappa(1:3,0)
        matrix(0,1:3) = matrix(1:3,0)

        ! Second row
        kappa(1,2) = K1(2,ifreq)*psio
        kappa(1,3) = K1(3,ifreq)*psio
        matrix(1,2) = -psim*K1(2,ifreq)
        matrix(1,3) = -psim*K1(3,ifreq)

        ! Third row
        kappa(2,1) = -kappa(1,2)
        kappa(2,3) = K2(3,ifreq)*psio
        matrix(2,1) = -matrix(1,2)
        matrix(2,3) = -psim*K2(3,ifreq)

        ! Fourth row
        kappa(3,1) = -kappa(1,3)
        kappa(3,2) = -kappa(2,3)
        matrix(3,1) = -matrix(1,3)
        matrix(3,2) = -matrix(2,3)

        ! Diagonal elements
        do is=0,3
          kappa(is,is) = 1d0
          matrix(is,is) = exu
        end do

        ! Invert kappa
        call Matinv(kappa)

        ! Compute the Stokes solution
        call MatVec(matrix,Stk(:,ifreq),vector1)
        call MatVec(kappa,vector1,vector2)
        call MatVec(kappa,(psim+psio)*S(:,ifreq),vector3)

        ! Compute final Stokes parameters
        Stk(:,ifreq) = vector2 + vector3

        !
        ! Control the physicallity of the result
        !

        ! For each Stokes parameter
        do iS=0,3

          ! Control underflow
          if (abs(Stk(iS,ifreq)).lt.TINYS) Stk(iS,ifreq) = 0d0

          ! Control NaN
          if (ieee_is_nan(Stk(iS,ifreq))) then

            ! Issue error
            write(umsg,'(A,i4,i5,A,'// &
                       'A,4(1x,es8.1),'// &
                       'A,5(1x,es8.1),'// &
                       'A,3(4(1x,es8.1),A))') &
                  'Error in RTCStep: NaN in Stokes parameters'// &
                  new_line('A')//'(ifreq)=(',ifreq,')', &
                  new_line('A')//'Stokes:',Stk(:,ifreq), &
                  new_line('A')//'etaI,S:',K0(0,ifreq),S(:,ifreq), &
                  new_line('A')//'Kabs'//new_line('A'), &
                               K0(:,ifreq),new_line('A'), &
                               K1(:,ifreq),new_line('A'), &
                               K2(:,ifreq),' '

            call abortedS(umsg,urou,.True.,.True.)

          endif ! if NaN

        end do ! Stokes parameters

        ! Control Physics
        if (Stk(0,ifreq).lt.0d0.or. &
            Stk(0,ifreq)*Stk(0,ifreq).lt. &
            (Stk(1,ifreq)*Stk(1,ifreq) + &
             Stk(2,ifreq)*Stk(2,ifreq) + &
             Stk(3,ifreq)*Stk(3,ifreq))) then

          ! Issue error
          write(umsg,'(A,i4,i5,A,'// &
                     'A,4(1x,es8.1),'// &
                     'A,5(1x,es8.1),'// &
                     'A,3(4(1x,es8.1),A))') &
                'Error in RTCStep: Non physical Stokes '// &
                'parameters'// &
                new_line('A')//'(ifreq)=(',ifreq,')', &
                new_line('A')//'Stokes:',Stk(:,ifreq), &
                new_line('A')//'etaI,S:',K0(0,ifreq),S(:,ifreq), &
                new_line('A')//'Kabs'//new_line('A'), &
                             K0(:,ifreq),new_line('A'), &
                             K1(:,ifreq),new_line('A'), &
                             K2(:,ifreq),' '
          call abortedS(umsg,urou,.not.nphysS,.True.)

        end if ! Physical Stokes parameters

      end do ! For each frequency

      return

      end subroutine RTCStep

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstep_mod
