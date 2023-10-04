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
!     04/20/2017
!  Last version:
!     09/21/2023 V3.1.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/21/2023:    V3.1.0 - Removed I+Q/I-Q solution algorithm, now
!                             solve the coupled equations. This solves
!                             an issue in which the Q source function
!                             was not really forced to be monotonic
!                             because S+ and S- were dominated by the
!                             well behaved I contribution (TdPA)
!
!     05/25/2023:    V3.0.4 - Removed P variables from the linear
!                             SC case error verbosity (TdPA)
!
!     05/16/2023:    V3.0.3 - When bezier fails, the formar solver
!                             tries linear interpolation (TdPA)
!
!     03/23/2023:    V3.0.2 - Bugfix: Ensured the OpenMP version
!                             compiles after some years of changes,
!                             albeit did not test it works (TdPA)
!
!     11/24/2022:    V3.0.1 - Added RTCStep (TdPA)
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
!     05/31/2019:    V1.1.2 - Typo in error message (TdPA)
!
!     05/08/2019:    V1.1.1 - There is a parameter that allows for
!                             wrong physics (TdPA)
!                           - Added more verbosity to error (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!                           - Now uses diexp (TdPA)
!                           - Now uses especific TINY variables (TdPA)
!
!     01/30/2019:    V1.0.4 - Changed axial by RTaxial (TdPA)
!
!     07/11/2017:    V1.0.3 - Bugfix:One of the conditional
!                             branchings in RTTau was wrong (TdPA)
!
!     07/07/2017:    V1.0.2 - Bugfix: The contribution function was
!                             simply wrong (TdPA)
!
!     05/05/2017:    V1.0.1 - RTStep know which point is computing
!                             during the call (TdPA)
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
!  RTStep:
!    Applies short characteristics to the points in input
!
!  RTTau:
!    Calculates the tau and stores where tau=1
!
!  RTContr:
!    Applies short characteristics to the points in input to calculate
!    the contribution function
!
!  RTCStep:
!    Applies short characteristics to solve a constant property
!  slab
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use math_mod
      use omp_mod
      use parameters_mod , only: vacuum,TINYS,TINYT,TINYCS
      use rtstepaux_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Advances the radiation transfer equation propagation using
      !! short characteristics with the BESSER interpolation
      !! algorithm.\n
      !!           iz(integer): Height index\n
      !!          ith(integer): Output direction polar index\n
      !!          iph(integer): Output direction azimuth index\n
      !!        nfreq(integer): Number of frequencies\n
      !!           ds1(dfloat): Geometrical distance backward\n
      !!           ds2(dfloat): Geometrical distance forward\n
      !!      K0M(dfloat(:,:)): First row absorption matrix previous
      !!                        point\n
      !!      K1M(dfloat(:,:)): Second row absorption matrix previous
      !!                        point\n
      !!      K2M(dfloat(:,:)): Third row absorption matrix previous
      !!                        point\n
      !!       SM(dfloat(:,:)): Emissivity previous point\n
      !!      K0O(dfloat(:,:)): First row absorption matrix current
      !!                        point\n
      !!      K1O(dfloat(:,:)): Second row absorption matrix current
      !!                        point\n
      !!      K2O(dfloat(:,:)): Third row absorption matrix current
      !!                        point\n
      !!       SO(dfloat(:,:)): Emissivity current point\n
      !!      K0P(dfloat(:,:)): First row absorption matrix next
      !!                        point\n
      !!       SP(dfloat(:,:)): Emissivity next point\n
      !!     StkM(dfloat(:,:)): Stokes parameters previous point\n
      !!     StkO(dfloat(:,:)): Stokes parameters current point\n
      !!     nolinear(logical): Bool to specify that the interpolation
      !!                        must not be linear
      subroutine RTStep(iz,ith,iph,nfreq,ds1,ds2,K0M,K1M,K2M,SM, &
                        K0O,K1O,K2O,SO,K0P,SP,StkM,StkO,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: iz,nfreq,ith,iph
      double precision, intent(in):: ds1,ds2
      double precision, dimension(0:3,nfreq), intent(in):: StkM, &
                                                           K0M,K1M, &
                                                           K2M,SM, &
                                                           K0O,K1O, &
                                                           K2O,SO, &
                                                           K0P,SP
      double precision, dimension(0:3,nfreq), intent(out):: StkO

      ! Local

      integer:: ifreq,iS,i1,tid=-1
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

!$omp parallel do default(none) &
!$omp private(iS,i1,ifreq,exu,dt1,dt2,psim,psio,om_m,om_o,om_c,C0) &
!$omp private(CA,vector1A,vector2A,vector3A,kappaA,matrixA) &
!$omp private(tid,umsg) &
!$omp shared(nolinear,iz,ith,iph,nfreq,ds1,ds2,nphysS,laborted) &
!$omp shared(K0M,K1M,K2M,SM,K0O,K1O,K2O,SO,K0P,SP,StkM,StkO,urou)
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

            call RTomega(exu,dt1,om_m,om_o,om_c)

            do iS=0,1

              call QBezierC0(dt1,dt2,SM(iS,ifreq),SO(iS,ifreq), &
                                     SP(iS,ifreq),CA(iS))

            end do

          endif

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
          if(nolinear)then

            call MatVec_2d(kappaA,om_m*SM(0:1,ifreq) + &
                                  om_o*SO(0:1,ifreq) + &
                                  om_c*CA,vector3A)

          ! If two points short characteristics
          else

            call MatVec_2d(kappaA,psim*SM(0:1,ifreq) + &
                                  psio*SO(0:1,ifreq), &
                           vector3A)

          end if

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
            if (isnan(StkO(iS,ifreq))) then
#ifdef _OPENMP
              tid = omp_get_thread_num() + 1
#endif
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,4(1x,es9.2)),'// &
                         '6(A,2(1x,es9.2)))') &
              'Error in RTStep: NaN in Stokes parameters'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'Stokes(M):',StkM(:,ifreq), &
              new_line('A')//'Stokes(O):',StkO(:,ifreq), &
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
              call abortedS(umsg,urou,tid,.True.,.True.)

            endif ! if NaN

          end do ! Stokes parameters

          ! Control Physics
          if (StkO(0,ifreq).lt.0d0.or. &
              StkO(0,ifreq)*StkO(0,ifreq).lt. &
              (StkO(1,ifreq)*StkO(1,ifreq) + &
               StkO(2,ifreq)*StkO(2,ifreq) + &
               StkO(3,ifreq)*StkO(3,ifreq))) then
#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
#endif
            write(umsg,'(A,3(i4,","),i5,A,'// &
                       '2(A,4(1x,es9.2)),'// &
                       '6(A,2(1x,es9.2)),A)') &
            'Non-phys. Stokes par.'// &
            new_line('A')//'(ith,iph,iz,ifreq)=(', &
            ith,iph,iz,ifreq,')', &
            new_line('A')//'Stokes(M):',StkM(:,ifreq), &
            new_line('A')//'Stokes(O):',StkO(:,ifreq), &
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

            ! Do not quit if non-linear
            if (nolinear) then
              call abortedS(umsg,urou,tid,.False.,.True.)
            ! Abort if already linear
            else
              call abortedS(umsg,urou,tid,.not.nphysS,.True.)
            end if

            ! Try linear
            call MatVec(kappaA,psim*SM(0:1,ifreq) + &
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

              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,4(1x,es9.2)),'// &
                         '6(A,2(1x,es9.2)))') &
              'Non-phys. Stokes par.'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'Stokes(M):',StkM(:,ifreq), &
              new_line('A')//'Stokes(O):',StkO(:,ifreq), &
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

              call abortedS(umsg,urou,tid,.not.nphysS,.True.)

            end if ! Physical Stokes parameters
          end if ! Physical Stokes parameters

        end do ! For each frequency
!$omp end parallel do

      !
      ! Full Stokes problem
      !
      else

!$omp parallel do default(none) &
!$omp private(iS,i1,ifreq,exu,dt1,dt2,psim,psio,om_m,om_o,om_c) &
!$omp private(C0,vector1,vector2,vector3,kappa,matrix,umsg,tid) &
!$omp shared(nolinear,iz,ith,iph,nfreq,ds1,ds2,nphysS,laborted) &
!$omp shared(K0M,K1M,K2M,SM,K0O,K1O,K2O,SO,K0P,SP,StkM,StkO,urou)
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

            call RTomega(exu,dt1,om_m,om_o,om_c)

            do iS=0,3

              call QBezierC0(dt1,dt2,SM(iS,ifreq),SO(iS,ifreq), &
                                     SP(iS,ifreq),C0(iS))

            end do

          endif

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
          if(nolinear)then

            call MatVec(kappa,om_m*SM(:,ifreq) + om_o*SO(:,ifreq) + &
                              om_c*C0,vector3)

          ! If two points short characteristics
          else

            call MatVec(kappa,psim*SM(:,ifreq) + psio*SO(:,ifreq), &
                        vector3)

          end if

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
            if (isnan(StkO(iS,ifreq))) then
#ifdef _OPENMP
              tid = omp_get_thread_num() + 1
#endif
              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,4(1x,es8.1)),'// &
                         '3(A,5(1x,es8.1)),'// &
                         '2(A,3(4(1x,es8.1),A)),'// &
                         'A,4(1x,es8.1))') &
              'NaN in Stokes parameters'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'Stokes(M):',StkM(:,ifreq), &
              new_line('A')//'Stokes(O):',StkO(:,ifreq), &
              new_line('A')//'etaI(M),S(M):', &
                          K0M(0,ifreq),SM(:,ifreq), &
              new_line('A')//'etaI(O),S(O):', &
                          K0O(0,ifreq),SO(:,ifreq), &
              new_line('A')//'etaI(P),S(P):', &
                          K0P(0,ifreq),SP(:,ifreq), &
              new_line('A')//'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(:,ifreq),new_line('A'), &
                           K2M(:,ifreq),new_line('A'), &
              'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(:,ifreq),new_line('A'), &
                           K2O(:,ifreq),new_line('A'), &
              'Kabs(P)'//new_line('A'), &
                           K0P(:,ifreq)

              call abortedS(umsg,urou,tid,.True.,.True.)

            endif ! if NaN

          end do ! Stokes parameters

          ! Control Physics
          if (StkO(0,ifreq).lt.0d0.or. &
              StkO(0,ifreq)*StkO(0,ifreq).lt. &
              (StkO(1,ifreq)*StkO(1,ifreq) + &
               StkO(2,ifreq)*StkO(2,ifreq) + &
               StkO(3,ifreq)*StkO(3,ifreq))) then
#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
#endif
            write(umsg,'(A,3(i4,","),i5,A,'// &
                       '2(A,4(1x,es8.1)),'// &
                       '3(A,5(1x,es8.1)),'// &
                       '2(A,3(4(1x,es8.1),A)),'// &
                       'A,4(1x,es8.1),A)') &
              'Non-phys. Stokes par.'// &
              new_line('A')//'(ith,iph,iz,ifreq)=(', &
              ith,iph,iz,ifreq,')', &
              new_line('A')//'Stokes(M):',StkM(:,ifreq), &
              new_line('A')//'Stokes(O):',StkO(:,ifreq), &
              new_line('A')//'etaI(M),S(M):', &
                          K0M(0,ifreq),SM(:,ifreq), &
              new_line('A')//'etaI(O),S(O):', &
                          K0O(0,ifreq),SO(:,ifreq), &
              new_line('A')//'etaI(P),S(P):', &
                          K0P(0,ifreq),SP(:,ifreq), &
              new_line('A')//'Kabs(M)'//new_line('A'), &
                           K0M(:,ifreq),new_line('A'), &
                           K1M(:,ifreq),new_line('A'), &
                           K2M(:,ifreq),new_line('A'), &
              'Kabs(O)'//new_line('A'), &
                           K0O(:,ifreq),new_line('A'), &
                           K1O(:,ifreq),new_line('A'), &
                           K2O(:,ifreq),new_line('A'), &
              'Kabs(P)'//new_line('A'), &
                           K0P(:,ifreq), &
              new_line('A')//' do lin.'

            ! Do not quit if non-linear
            if (nolinear) then
              call abortedS(umsg,urou,tid,.False.,.True.)
            else
              call abortedS(umsg,urou,tid,.not.nphysS,.True.)
            end if

            ! Linear
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

              write(umsg,'(A,3(i4,","),i5,A,'// &
                         '2(A,4(1x,es8.1)),'// &
                         '2(A,5(1x,es8.1)),'// &
                         '2(A,3(4(1x,es8.1),A)))') &
                'Nonphys. Stokes par.'// &
                new_line('A')//'(ith,iph,iz,ifreq)=(', &
                ith,iph,iz,ifreq,')', &
                new_line('A')//'Stokes(M):',StkM(:,ifreq), &
                new_line('A')//'Stokes(O):',StkO(:,ifreq), &
                new_line('A')//'etaI(M),S(M):', &
                            K0M(0,ifreq),SM(:,ifreq), &
                new_line('A')//'etaI(O),S(O):', &
                            K0O(0,ifreq),SO(:,ifreq), &
                new_line('A')//'Kabs(M)'//new_line('A'), &
                             K0M(:,ifreq),new_line('A'), &
                             K1M(:,ifreq),new_line('A'), &
                             K2M(:,ifreq),new_line('A'), &
                'Kabs(O)'//new_line('A'), &
                             K0O(:,ifreq),new_line('A'), &
                             K1O(:,ifreq),new_line('A'), &
                             K2O(:,ifreq),new_line('A')

              call abortedS(umsg,urou,tid,.not.nphysS,.True.)

            end if ! Physical Stokes parameters
          end if ! Physical Stokes parameters

        end do ! For each frequency
!$omp end parallel do

      end if ! Method of solution

      return

      end subroutine RTStep

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds to the optical depth the contribution of a section
      !! between two points and checks if the optical depth becomes
      !! one.\n
      !!        nfreq(integer): Number of frequencies\n
      !!           ds1(dfloat): Geometrical distance backward\n
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
      subroutine RTTau(nfreq,ds1,z1,z2,etaM,etaO,tauM,tauO,tau1)

      ! I/O

      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,z1,z2
      double precision, dimension(:), intent(in):: etaM,etaO,tauM
      double precision, dimension(:), intent(inout):: tauO
      double precision, dimension(:,:), intent(inout):: tau1

      ! Local

      integer:: ifreq
      double precision:: a,b,D,dt1

!$omp parallel do default(none)  &
!$omp private(ifreq,a,b,D,dt1) &
!$omp shared(nfreq,ds1,z1,z2,etaM,etaO,tauM,tauO,tau1)
      do ifreq=1,nfreq

        !
        ! Calculate optical distance between the two points
        !
        call ftau(etaM(ifreq),etaO(ifreq),ds1,dt1)

        ! Add this to the optical depth
        tauO(ifreq) = tauM(ifreq) + dt1
        tau1(1,ifreq) = tauO(ifreq)

        ! If the optical was lower than one before
        if (tauM(ifreq).lt.1) then

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

       end do ! Frequencies
!$omp end parallel do

      end subroutine RTTau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the contribution function in a step of the short
      !! characteristics.\n
      !!        nfreq(integer): Number of frequencies\n
      !!           ds1(dfloat): Geometrical distance backward\n
      !!           ds2(dfloat): Geometrical distance forward\n
      !!           dz1(dfloat): Geometrical vertical distance
      !!                        backward\n
      !!           dz2(dfloat): Geometrical vertical distance
      !!                        forward\n
      !!      K0M(dfloat(:,:)): First row absorption matrix previous
      !!                        point\n
      !!      K0O(dfloat(:,:)): First row absorption matrix current
      !!                        point\n
      !!      K1O(dfloat(:,:)): Second row absorption matrix current
      !!                        point\n
      !!      K2O(dfloat(:,:)): Third row absorption matrix current
      !!                        point\n
      !!       SO(dfloat(:,:)): Emissivity current point\n
      !!      K0P(dfloat(:,:)): First row absorption matrix next
      !!                        point\n
      !!     StkO(dfloat(:,:)): Stokes parameters current point\n
      !!        tau(dfloat(:)): Optical depth at current point\n
      !!    Contr(dfloat(:,:)): Contribution function at current
      !!                        point\n
      !!     nolinear(logical): Bool to specify that the interpolation
      !!                        must not be linear
      subroutine RTContr(nfreq,ds1,ds2,dz1,dz2,K0M, &
                         K0O,K1O,K2O,SO,K0P, &
                         StkO,tau,Contr,nolinear)

      ! I/O

      logical, intent(in):: nolinear
      integer, intent(in):: nfreq
      double precision, intent(in):: ds1,ds2,dz1,dz2
      double precision, dimension(:), intent(in):: tau
      double precision, dimension(0:3,nfreq), intent(in):: StkO,K0M, &
                                                           K0O,K1O, &
                                                           K2O,SO, &
                                                           K0P
      double precision, dimension(0:3,nfreq), intent(out):: Contr

      ! Local

      integer:: iS,ifreq,tid=-1

      double precision:: dtau,exu,dt1,dt2
      double precision, dimension(0:3):: vector1
      double precision, dimension(0:3,0:3):: kappa


      ! Routine name
      urou = 'RTContr'

!$omp parallel do default(none) &
!$omp private(iS,ifreq,exu,dt1,dt2,dtau,vector1,kappa,tid,umsg) &
!$omp shared(nolinear,nfreq,ds1,ds2,dz1,dz2,RTaxial,laborted) &
!$omp shared(K0M,K0O,K1O,K2O,SO,K0P,StkO,tau,Contr,urou)
      do ifreq=1,nfreq

        !
        ! Calculate optical distance with the previous point
        !
        call ftau(K0M(0,ifreq),K0O(0,ifreq),ds1,dt1)
        dtau = dt1/dz1

        !
        ! Calculate optical distance with the next point and the
        ! optical depth derivative
        !
        if(nolinear)then
          call ftau(K0O(0,ifreq),K0P(0,ifreq),ds2,dt2)
          dtau = .5d0*(dtau + dt2/dz2)
        end if

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
          Contr(0,ifreq) = (SO(0,ifreq) - &
                            vector1(0))*exu*dtau/K0O(0,ifreq)
          Contr(1,ifreq) = (SO(1,ifreq) - &
                            vector1(1))*exu*dtau/K0O(0,ifreq)
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
          if (isnan(Contr(iS,ifreq))) then
#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
#endif
            write(umsg,'(2(A,1x,es9.2),'// &
                       'A,2(1x,es9.2),'// &
                       'A,1x,es9.2)') &
            'Error in RTContr: NaN in Contribution Function'// &
            new_line('A')//'ContrI:',Contr(0,ifreq), &
            new_line('A')//'etaI(M):',K0M(0,ifreq), &
            new_line('A')//'etaI(O),SI(O):', &
                        K0O(0,ifreq),SO(0,ifreq), &
            new_line('A')//'etaI(P):',K0P(0,ifreq)

            call abortedS(umsg,urou,tid,.True.,.True.)

          end if ! If NaN

        end do ! Stokes parameters
      end do ! For each frequency
!$omp end parallel do

      end subroutine RTContr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the slab RT equation with short characteristics
      !1 assuming constant properties.\n
      !!        nfreq(integer): Number of frequencies\n
      !!          fact(dfloat): Optical depth scale factor\n
      !!       K0(dfloat(:,:)): First row absorption matrix\n
      !!       K1(dfloat(:,:)): Second row absorption matrix\n
      !!       K2(dfloat(:,:)): Third row absorption matrix\n
      !!        S(dfloat(:,:)): Emissivity current point\n
      !!      Stk(dfloat(:,:)): Stokes parameters
      subroutine RTCStep(nfreq,fact,K0,K1,K2,S,Stk)

      ! I/O

      integer, intent(in):: nfreq
      double precision, intent(in):: fact
      double precision, dimension(0:3,nfreq), intent(in):: K0,K1, &
                                                           K2,S
      double precision, dimension(0:3,nfreq), intent(inout):: Stk

      ! Local

      integer:: ifreq,iS,tid=-1

      double precision:: dt,exu,psim,psio
      double precision, dimension(0:3):: vector1,vector2,vector3
      double precision, dimension(0:3,0:3):: kappa,matrix


      ! Routine name
      urou = 'RTCStep'

!$omp parallel do default(none) &
!$omp private(ifreq,dt,psim,psio,kappa,matrix,exu) &
!$omp private(vector1,vector2,vector3,iS,tid,umsg) &
!$omp shared(nfreq,fact,laborted,K0,K1,K2,S,Stk,urou,nphysS)
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
          if (isnan(Stk(iS,ifreq))) then
#ifdef _OPENMP
            tid = omp_get_thread_num() + 1
#endif
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

            call abortedS(umsg,urou,tid,.True.,.True.)

          endif ! if NaN

        end do ! Stokes parameters

        ! Control Physics
        if (Stk(0,ifreq).lt.0d0.or. &
            Stk(0,ifreq)*Stk(0,ifreq).lt. &
            (Stk(1,ifreq)*Stk(1,ifreq) + &
             Stk(2,ifreq)*Stk(2,ifreq) + &
             Stk(3,ifreq)*Stk(3,ifreq))) then
#ifdef _OPENMP
          tid = omp_get_thread_num() + 1
#endif
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

          call abortedS(umsg,urou,tid,.not.nphysS,.True.)

        end if ! Physical Stokes parameters

      end do ! For each frequency
!$omp end parallel do

      return

      end subroutine RTCStep

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtstep_mod
