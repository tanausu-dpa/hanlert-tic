      !> Statistical equilibrium equations
      module see_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/26/2017
!  Last version:
!     09/29/2023 V3.0.7
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.7 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!
!     07/03/2023:    V3.0.6 - Added an alternative trace equation
!                             which fixes the lower term population
!                             but not ensures the particle
!                             conservation (TdPA)
!
!     03/08/2023:    V3.0.5 - Bugfix: iJJ1 cannot be limited to
!                             use with zconj1 if iJJ is also limited
!                             to use with zconj. Removed zconj1 and
!                             extended to the full iJJ1 loop (TdPA)
!                           - The normalization condition for the
!                             atomic trace is now imposed at the
!                             end of SEbuild (TdPA)
!
!     02/01/2023:    V3.0.4 - Added an argument with the number of
!                             transitions to rA and rS (TdPA)
!
!     11/24/2022:    V3.0.3 - Removed non-used variables (TdPA)
!
!     11/10/2022:    V3.0.2 - Added the option to zero out the
!                             density matrix of the last ion of the
!                             atom. This is useful when using
!                             two-level atoms while keeping the
!                             ionization stage to compute the
!                             broadening parameters (TdPA)
!
!     10/26/2022:    V3.0.1 - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS and generation
!                             of the error message (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Stopped using Atom%rho and created a
!                             local variable to be able to call
!                             SEE with threads (TdPA)
!
!     02/17/2021:    V1.3.5 - Added argument to forbidden collisional
!                             transfer calls (TdPA)
!                           - Added back the factorials and signs
!                             argument to forbidden collisional
!                             transfer calls (TdPA)
!
!     10/18/2019:    V1.3.4 - Possibility to keep populations fixed
!                             if specified in the input (TdPA)
!                           - Added initrho routine (TdPA)
!
!     05/08/2019:    V1.3.3 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!                           - There is a parameter that allows for
!                             wrong physics (TdPA)
!
!     04/12/2019:    V1.3.2 - Changed ifill (integer) to lfill
!                             (logical). They have the same use (TdPA)
!
!     03/18/2019:    V1.3.1 - Physical errors do not make rhosol to
!                             return (TdPA)
!
!     02/20/2019:    V1.3.0 - New verbosity (TdPA)
!                           - Now uses especific TINY variables (TdPA)
!
!     11/14/2017:    V1.2.4 - Bugfix: Relaxation rates were called
!                             for photoionizations even when there
!                             were not in the model atom (TdPA)
!
!     10/30/2017:    V1.2.3 - Using the multilevel flag instead of
!                             the comparison with S (TdPA)
!                           - By popular demand, stored 1/4 into a
!                             parameter (what is used to determine
!                             if two quantum numbers are the
!                             same (TdPA)
!
!     10/25/2017:    V1.2.2 - Bugfix: For the multilevel case, forgot
!                             to check K diagonality (TdPA)
!
!     10/13/2017:    V1.2.1 - Improved information given when see
!                             solution is not physical when it is not
!                             diagonal in J (TdPA)
!
!     10/11/2017:    V1.2.0 - Able to consider multi-level if spin
!                             is zero (TdPA)
!
!     09/22/2017:    V1.1.0 - Possibility to limit number K (TdPA)
!
!     09/14/2017:    V1.0.5 - Bugfix: The way to force rhoKQ=0 in
!                             isolted ground states was wrong because
!                             of the application of conjugation
!                             properties (TdPA)
!
!     09/08/2017:    V1.0.4 - Added checks for physicallity of the
!                             solution of the SEE (TdPA)
!                           - The ionizing collisional rates are
!                             only taken into account for K=0 (TdPA)
!                           - Solved a memory problem with the
!                             forbidden collisions when there were
!                             no collisions (TdPA)
!
!     06/15/2017:    V1.0.3 - Added check L=L1=0 in SEbuild (TdPA)
!                           - Added check S=SS in SEbuild (TdPA)
!                           - Bugfix: rLL used before defining (TdPA)
!
!     06/14/2017:    V1.0.2 - In densmatrI, rho is intent(inout), to
!                             avoid valgrind complains (TdPA)
!
!     05/05/2017:    V1.0.1 - Fixed multilevel part, the conditional
!                             logic was wrong (TdPA)
!
!     04/26/2017:    V1.0.0 - First version (TdPA)
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
!  SEE:
!    Manage the solution of the statistical equilibrium equations
!
!  SEbuild:
!    This subroutine builds the coefficient matrix for the system
!  of S.E. equations, for the multi-term atom in the L-S coupling
!  approximation (real version)
!
!  initrho:
!    Initializes rho vector or considers fixed populations and zero
!  ion
!
!  densmatr:
!    This subroutine calculates the density matrix solution
!  of the SE equations
!
!  rhosol:
!    Restore the solution density matrix
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : k2f , TINYR
      use seeaux_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calls subroutines to solve the statistical equilibrium
      !! equations.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       JRad(dcmplx(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!      JRadS(dcmplx(:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!        JPhot(dfloat(:,:)): Mean intensity integrated with
      !!                            photoionizations\n
      !!            larmor(dfloat): Magnetic field in larmor frequency
      !!                            units\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!              tid(integer): thread index
      subroutine SEE(Atom,JRad,JRadS,Jphot,larmor,Flgsg,iz,tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,tid
      double precision, intent(in):: larmor
      double precision, dimension(:,:),intent(in):: Jphot
      complex(kind=8), dimension(-2:2,0:2,Atom%ntran), &
                                             intent(in):: JRad,JRadS

      ! Local
      double precision:: rho(Atom%ndim)
      double precision:: STcoeff(Atom%ndim,Atom%ndim)


      !
      ! Build SEE equations for multi-term atom
      !
      call SEbuild(Atom,JRad(:,:,1:Atom%ntran), &
                   JRadS(:,:,1:Atom%ntran),Jphot, &
                   larmor,STcoeff,Flgsg,iz,tid)


      !
      ! Initialize rho and fix populations if requested
      !
      call initrho(Atom,STcoeff,rho,iz)


      !
      ! Solve the SEE
      !
      call densmatr(rho,Atom%ndim,STcoeff)

      !
      ! Rearrange the solution into the rhoKQ matrices
      !
      call rhosol(Atom,rho,Flgsg,iz,tid)

      end subroutine SEE


!#####################################################################
!#####################################################################
!#####################################################################

      !> Builds the statistical equilibrium equations system\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!       RadJ(dcmplx(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!      RadJS(dcmplx(:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!           JP(dfloat(:,:)): Mean intensity integrated with
      !!                            photoionizations\n
      !!            larmor(dfloat): Magnetic field in larmor frequency
      !!                            units\n
      !!      STcoeff(dfloat(:,:)): Statistical equilibrium equations
      !!                            system\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!              tid(integer): Thread index
      subroutine SEbuild(Atom,RadJ,RadJS,JP,larmor,STcoeff, &
                         Flgsg,iz,tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,tid
      double precision, intent(in):: larmor
      double precision, dimension(:,:),intent(in):: JP
      double precision, dimension(:,:), intent(out):: STcoeff
      complex(kind=8),dimension(-2:2,0:2,Atom%ntran), &
                                             intent(in):: RadJ,RadJS

      ! Local

      logical:: zJ,zJ1,zK,zQ,zJ_JJ,zJ1_JJ1,zK1KK,zK2KK,zJlt,zJgt,zJeq
      logical:: zpermit,zpermitJ,zRpermit,zCpermit,zRpermitJ,zCpermitJ
      logical:: zrelax,zupper,zlower,zrelaxJ,zupperJ,zlowerJ,zion
      logical:: zflag,zreal,zzreal,zTS,zTA,zR,zdiag,zconj
      logical:: zMK_J,zMK
      logical, dimension(Atom%ndim,Atom%ndim):: lfill

      integer:: itran,icol,iphot,iterm,itterm,iJ,iJ1,iJJ,iJJ1
      integer:: K,iQ,KK,iQQ,i,ii,i1,ii1,ir,iir,mKr,mKc

      double precision, parameter:: dSLJ = .25d0

      double precision:: rL,S,rJ,rJ1,rLL,SS,rJJ,rJJ1,rK,Q,rKK,QQ
      double precision:: rJsgn,tJsgn,rJJsgn,tJJsgn,dFS
      double precision:: rMKcoeff,Dcoeff,rSCcoeff,rACcoeff
      double precision:: tSCcoeff,tACcoeff,rEcoeff,tEcoeff
      double precision:: rSPcoeff,rAPcoeff,tSPcoeff,tAPcoeff
      double precision:: rEPcoeff,tEPcoeff
      double precision, dimension(Atom%ndim):: rho
      double precision, dimension(Atom%ndim,Atom%ndim):: Scoeff
      double precision, dimension(Atom%ndim,Atom%ndim):: Tcoeff

      complex(kind=8):: rScoeff,rAcoeff
      complex(kind=8):: tScoeff,tAcoeff


      ! Routine name
      urou = 'SEbuild'


      !
      ! Initialize flag matrix and SEE coefficients
      !
      lfill = .True.
      Scoeff = 0d0
      Tcoeff = 0d0
      STcoeff = 0d0
      rho = 0d0

      ! Reset the element indexes
      i = 0
      ii = 0


      !
      ! Build SEE
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! Get term quantities
        rL = Atom%rLval(iterm)
        S = Atom%Sval(iterm)

        ! Spontaneous relaxation b-b rate
        call rE(iterm,Atom%Ecoeff,rEcoeff)

        ! For each level (row_a)
        do iJ=1,Atom%nJ(iterm)

          ! Get level momentum
          rJ = Atom%rJval(iJ,iterm)

          ! Get level index
          ir = Atom%irho(iterm)%irho_ij(iJ)

          ! We use the conjugation properties of the SE matrix
          ! (see zconj), so we can restrict the range of iJ1
          ! For each level (row_b)
          do iJ1=iJ,Atom%nJ(iterm)  ! use with zconj
         !do iJ1=1,Atom%nJ(iterm)
         !do iJ1=iJ,iJ              ! no J,J'

            ! Get level momentum
            rJ1 = Atom%rJval(iJ1,iterm)

            ! Get energy diference in the appropriate unnits
            dFS = k2f*(Atom%FSfreq(iJ,iterm) - &
                       Atom%FSfreq(iJ1,iterm))

            ! Get sign from the difference of J,J'
            rJsgn = Flgsg%sg(nint(rJ-rJ1))

            ! For each K
            do K=nint(abs(rJ-rJ1)), &
                 min(nint(rJ+rJ1),Atom%Kcut(iterm))

              ! Get real number
              rK = dble(K)


              !
              ! Depolarizing collisions
              !

              ! If not populations
              if (K.ne.0) then

                ! Get rate from atomic data
                Dcoeff = Atom%gk(iterm,iJ,iJ1,K,iz)

              ! There are not depolarizing collisions for populations
              else

                Dcoeff = .0D0

              end if ! K!=0

              ! For each Q
              do iQ=-K,K

                ! Get real number
                Q = dble(iQ)

                ! Get sign from the difference of J,J' and Q
                tJsgn = rJsgn*Flgsg%sg(iQ)

                ! Get the indexes of the elements J,J' and J',J
                i = Atom%irho(iterm)%Jrho(iJ1,iJ)%kq(iQ,K)
                i1 = Atom%irho(iterm)%Jrho(iJ,iJ1)%kq(-iQ,K)

                ! Flag of filled line
                zflag = .True.

                ! If J!=J' or Q!=0, use conjugation properties
                zconj = (iJ.ne.iJ1).or.(iQ.ne.0)

                ! For each term (column)
                do itterm=1,Atom%nMulti

                  ! Get the term quantities
                  rLL = Atom%rLval(itterm)
                  SS = Atom%Sval(itterm)

                  ! Check if there is a permitted transition between
                  ! these terms
                  itran = Atom%irad(itterm,iterm)
                  zRpermit = Atom%irad(itterm,iterm).ne.0.and. &
                             nint(abs(rL-rLL)).le.1.and. &
                             nint(rL+rLL).gt.0.and. &
                             nint(abs(S-SS)).lt.1

                  ! Check if there is a permitted collision between
                  ! these terms
                  icol = Atom%icol(itterm,iterm)
                  zCpermit = Atom%icol(itterm,iterm).ne.0.and. &
                             nint(abs(rL-rLL)).le.1.and. &
                             nint(rL+rLL).gt.0.and. &
                             nint(abs(S-SS)).lt.1

                  ! Check how the column relates to the row
                  zrelax = itterm.eq.iterm
                  zupper = itterm.gt.iterm
                  zlower = itterm.lt.iterm

                  ! Check ionization
                  zion = abs(Atom%stage(iterm) - &
                             Atom%stage(itterm)).gt.0

                  ! Check if there is any permitted transition between
                  ! these terms
                  zpermit = zRpermit.or.zCpermit

                  ! For each level (column_a)
                  do iJJ=1,Atom%nJ(itterm)

                    ! Get level momentum
                    rJJ = Atom%rJval(iJJ,itterm)

                    ! Get level index
                    iir = Atom%irho(itterm)%irho_ij(iJJ)

                    ! Check if J = J''
                    zJ = abs(rJ-rJJ).lt.dSLJ

                    ! Check if J and J'' differ in less than 1 unit
                    zJ_JJ = nint(abs(rJ-rJJ)).le.1

                    ! For each level (column_b)
                   !do iJJ1=iJJ,Atom%nJ(itterm) ! use with zconj1
                    do iJJ1=1,Atom%nJ(itterm)
                   !do iJJ1=iJJ,iJJ             ! no J,J'

                      ! Get level momentun
                      rJJ1 = Atom%rJval(iJJ1,itterm)

                      ! Get sign given J'' and J'''
                      rJJsgn = Flgsg%sg(nint(rJJ-rJJ1))

                      ! Check if J' = J'''
                      zJ1 = abs(rJ1-rJJ1).lt.dSLJ

                      ! Check if J' and J''' differ in less than 1
                      zJ1_JJ1 = nint(abs(rJ1-rJJ1)).le.1

                      ! Determine if there is a permitted FS trans.
                      zRpermitJ = zRpermit.and.zJ_JJ.and.zJ1_JJ1

                      ! Determine if there is a permitted FS coll.
                      zCpermitJ = zCpermit.and.zJ_JJ.and.zJ1_JJ1

                      ! Determine if there is any kind of permitted
                      ! transition
                      zpermitJ = zRpermitJ.or.zCpermitJ

                      ! Check if one of the pairs of J differs in one,
                      ! but the other is the same
                      zMK_J = (zJ_JJ.and.zJ1).or.(zJ1_JJ1.and.zJ)

                      if (zupper.and.zRpermitJ) &
                      call tE(rL,rJ,rJ1,rK,rLL,rJJ,rJJ1,S, &
                              Atom%Ecoeff(itterm,iterm),Flgsg,tEcoeff)

                      ! For each K'
                      do KK=nint(abs(rJJ-rJJ1)), &
                            min(nint(rJJ+rJJ1),Atom%Kcut(itterm))

                        ! Get a real number
                        rKK = dble(KK)

                        ! Check if K = K'
                        zK = K.eq.KK

                        ! Check if K and K' differ in less than 2
                        zK1KK = abs(K-KK).le.1

                        ! Check if K and K' differ in less than 3
                        zK2KK = abs(K-KK).le.2

                        ! Check if it is a diagonal element
                        zdiag = (zJ.and.zJ1).and.zK

                        ! Check if the conditions of the magnetic
                        ! kernel are satisfied
                        zMK = zMK_J.and.zK1KK

                        ! Check if stimulated transfer rate is to be
                        ! calculated from the flag and K values
                        zTS = stm.and.zK2KK

                        ! Check if absorption transfer rate is to be
                        ! calculated from the term order and K values
                        zTA = zlower.and.zK2KK

                        ! Check if a relaxation rate is to be
                        ! calculated from the J and K values
                        zR = (zJ.or.zJ1).and.zK2KK

                        ! For each Q'
                        do iQQ=-KK,KK

                          ! Get a real number
                          QQ = dble(iQQ)

                          ! Get the sign given all the J and Q
                          tJJsgn = tJsgn*rJJsgn*Flgsg%sg(iQQ)

                          ! Check if Q = Q'
                          zQ = iQ.eq.iQQ

                          ! Get the indexes of the elements J'',J'''
                          !  and J''',J''
                          ii = Atom%irho(itterm)% &
                                    Jrho(iJJ1,iJJ)%kq(iQQ,KK)
                          ii1 = Atom%irho(itterm)% &
                                     Jrho(iJJ,iJJ1)%kq(-iQQ,KK)

                          ! If J!=J' or Q!=0, use conjugation
                          ! properties
                         !zconj1 = (iJJ.ne.iJJ1).or.(iQQ.ne.0)

      !
      ! Reset the Identation
      !

      !
      ! If we have not touched this element yet, proceed
      !
      if (lfill(i,ii)) then

        ! Transition rates are nonvanishing only for radiatively
        ! permitted transitions between the pairs (iJ,iJ1) and
        ! (iJJ,iJJ1)
        if (zpermitJ) then

          ! Contribute the transition rate for emission (spontaneous
          ! and stimulated) from the upper terms
          if (zupper) then

            ! If there is a permitted radiative transition
            if(zRpermitJ)then

              ! If diagonal in K and Q, spontaneous transfer rate
              if (zK.and.zQ) Scoeff(i,ii) = tEcoeff

              ! If stimulated emission transfer rate to calculate
              if (zTS) then

                ! Maximum Kr
                mKr = Atom%Krad(itran)

                ! Stimulated emission transfer rate
                call tS(rL,rJ,rJ1,rK,Q, &
                        rLL,rJJ,rJJ1,rKK,QQ,S, &
                        Atom%Ecoeff(iterm,itterm), &
                        RadJS(:,:,itran),mKr,Flgsg,tScoeff)

                ! Add the contribution to the SEE
                Scoeff(i,ii) = dble(tScoeff) + Scoeff(i,ii)
                Tcoeff(i,ii) = dimag(tScoeff) + Tcoeff(i,ii)

              end if ! stimulated emission
            end if ! There is a permitted radiative transition

            ! If K and K' differ in 2 as maximum and there is a
            ! permitted collisional transition
            if (zK2KK.and.zCpermitJ) then

              ! De-excitation collisional transfer rate
              call tSC(rL,rJ,rJ1,rK,Q,rLL,rJJ,rJJ1,rKK,QQ,S, &
                       Atom%Ccoeff(itterm,iterm,iz),Flgsg,tSCcoeff)

              ! Add the contribution to the SEE
              Scoeff(i,ii) = tSCcoeff + Scoeff(i,ii)

            end if ! If K-K' <= 2 and there is a collisional trans.
          end if ! If column > row

          ! If absorption transfer rate to be calculated
          if (zTA) then

            ! If there is a radiative permitted transition
            if (zRpermitJ) then

              ! Maximum Kr
              mKr = Atom%Krad(itran)

              ! Maximum Kcut
              mKc = max(Atom%Kcut(iterm),Atom%Kcut(itterm))

              ! Absorption transfer rate
              call tA(rL,rJ,rJ1,rK,Q, &
                      rLL,rJJ,rJJ1,rKK,QQ,S, &
                      Atom%Ecoeff(itterm,iterm), &
                      RadJ(:,:,itran),mKr,mKc, &
                      Flgsg,tAcoeff)

              ! Add contribution to SEE
              Scoeff(i,ii) = dble(tAcoeff)
              Tcoeff(i,ii) = dimag(tAcoeff)

            end if ! permitted radiative transition

            ! If there is a permitted collisional transition
            if (zCpermitJ) then

              ! Excitation collisional transfer rate
              call tAC(rL,rJ,rJ1,rK,Q,rLL,rJJ,rJJ1,rKK,QQ,S, &
                       Atom%Ccoeff(itterm,iterm,iz),Flgsg,tACcoeff)

              ! Add the contribution to the SEE
              Scoeff(i,ii)= Scoeff(i,ii) + tACcoeff

            end if ! If there is a permitted collisional transition
          end if ! If absorption transfer rate to be calculated
        end if ! If there is any permitted transition

        ! If it is a diagonal element (in terms)
        if (zrelax) then

          ! If Q=Q'
          if (zQ) then

            ! If magnetic kernel to be calculated
            if (zMK) then

              ! Multilevel
              if (Atom%ML) then

                ! If diagonal in K
                if (zK) then

                    rMKcoeff = Atom%gL(iterm)*larmor*Q

                ! Non-diagonal in K
                else

                    rMKcoeff = 0d0

                end if

              ! Multiterm
              else

                ! Calculate magnetic kernel
                call MK(rL,S,rJ,rJ1,rK,rJJ,rJJ1, &
                        rKK,Q,zJ,zJ1,zK,dFS,larmor,Flgsg, &
                        rMKcoeff)

              end if ! Multilevel or multiterm

              ! Add contribution to SEE
              Tcoeff(i,ii) = -rMKcoeff

            end if ! magnetic kernel to be calculated

            ! If completely diagonal, add spontaneous relaxation
            ! and depolarizing collisions
            if (zdiag) Scoeff(i,ii) = -rEcoeff - Dcoeff

          end if ! Q=Q'

          ! If diagonal in one of the J pairs and K-K'<=2
          if (zR) then

            ! If stimulated emission
            if (stm) then

              ! Calculate stimulated emission relaxation rate
              call rS(iterm,Atom%irad(:,iterm),Atom%ntran, &
                      Atom%rLval,rL,S,rJ,rJ1,rK,Q,rJJ,rJJ1,rKK,QQ, &
                      zJ,zJ1,Atom%Ecoeff(:,iterm),RadJS,Atom%Krad, &
                      Flgsg,rScoeff)

            ! If not stimulated emission
            else

              ! Force 0
              rScoeff = dcmplx(0d0,0d0)

            end if ! stimulated emission

            ! Calculate absorption relaxation rate
            call rA(iterm,Atom%irad(:,iterm),Atom%ntran, &
                    Atom%rLval,Atom%nMulti,rL,S,rJ,rJ1,rK,Q,rJJ, &
                    rJJ1,rKK,QQ,zJ,zJ1,Atom%Ecoeff(iterm,:),RadJ, &
                    Atom%Krad,Atom%Kcut,Flgsg,rAcoeff)

            ! Add contribution to SEE
            Scoeff(i,ii) = -dble(rScoeff) - dble(rAcoeff) + &
                           Scoeff(i,ii)
            Tcoeff(i,ii) = -dimag(rScoeff) - dimag(rAcoeff) + &
                           Tcoeff(i,ii)

            ! Calculate de-excitation collisional transfer rate
            call rSC(iterm,Atom%icol(:,iterm),Atom%rLval,rL,S,rJ, &
                     rJ1,rK,Q,rJJ,rJJ1,rKK,QQ,zJ,zJ1, &
                     Atom%Ccoeff(iterm,:,iz),Flgsg,rSCcoeff)

            ! Calculate excitation collisional transfer rate
            call rAC(iterm,Atom%icol(:,iterm),Atom%rLval, &
                     Atom%nMulti,rL,S,rJ,rJ1,rK,Q,rJJ,rJJ1,rKK, &
                     QQ,zJ,zJ1,Atom%Ccoeff(iterm,:,iz),Flgsg, &
                     rACcoeff)

            ! Add contribution to SEE
            Scoeff(i,ii)= -rSCcoeff - rACcoeff + Scoeff(i,ii)

          end if ! if diagonal in one pair of J and K-K'<=2
        end if ! If it is a diagonal element (in terms)

        ! If diagonal in K, Q, and in rhoKQ coherences
        if (zK.and.zQ.and.(iJ.eq.iJ1).and.(iJJ.eq.iJJ1)) then

          ! Check relative positions of levels
          zrelaxJ = iir.eq.ir

          !
          ! Multi-level contributions
          !

          ! If column = row
          if (zrelaxJ) then

            ! If there are collisions
            if (maxval(Atom%CcoeffJ(ir,:,iz)).gt.0d0) then

              ! De-excitation collisional relaxation rate
              call rSFC(ir,K,Atom%CcoeffJ(ir,:,iz), &
                        Atom%fcflag(:,ir),rSCcoeff)

              ! Excitation collisional relaxation rate
              call rAFC(ir,K,Atom%nlevel,Atom%CcoeffJ(ir,:,iz), &
                        Atom%fcflag(:,ir),rACcoeff)

              ! Add contribution to SEE
              Scoeff(i,ii) = -rSCcoeff - rACcoeff + Scoeff(i,ii)

            end if ! Collisions

            ! If dealing with populations
            if (K.eq.0) then

              ! If there are photoionizations
              if (maxval(Atom%iphot(:,ir)).gt.0) then

                ! Spontaneous relaxation b-f rate
                call rEP(ir,Atom%iphot(:,ir),Atom%phot,iz,rEPcoeff)

                ! Stimulated b-f relaxation rate
                call rSP(ir,Atom%iphot(:,ir),JP(:,2),rSPcoeff)

                ! Absorption b-f relaxation rate
                call rAP(ir,Atom%iphot(:,ir),Atom%nlevel, &
                         JP(:,1),rAPcoeff)

                Scoeff(i,ii)= -rSPcoeff - rAPcoeff - rEPcoeff + &
                               Scoeff(i,ii)

              end if

            end if ! populations

          ! If not relaxation rate
          else

            ! Check relative positions of levels
            zupperJ = iir.gt.ir
            zlowerJ = iir.lt.ir

            ! Check there are collisions
            if (Atom%CcoeffJ(iir,ir,iz).gt.0d0) then

              ! Check there is a forbidden collision
              if (Atom%fcflag(iir,ir).gt.0.and. &
                  (K.eq.0.or..not.zion)) then

                ! If column > row
                if (zupperJ) then

                  ! De-excitation collisional transfer rate
                  call tSFC(rJ,rJJ,rK,Atom%CcoeffJ(iir,ir,iz), &
                            Flgsg,fcol_transfer,tSCcoeff)

                  ! Add contribution to SEE
                  Scoeff(i,ii) = tSCcoeff + Scoeff(i,ii)

                ! If column < row
                else if (zlowerJ) then

                  ! Excitation collisional transfer rate
                  call tAFC(rJ,rJJ,rK,Atom%CcoeffJ(iir,ir,iz), &
                            Flgsg,fcol_transfer,tACcoeff)

                  ! Add contribution to SEE
                  Scoeff(i,ii) = Scoeff(i,ii) + tACcoeff

                end if ! Direction of transition

              end if ! Forbidden transfer rate

            end if ! Collisions

            ! If dealing with populations
            if (K.eq.0) then

              ! Check if there is a b-f transition between levels
              iphot = Atom%iphot(iir,ir)

              ! If there is a b-f transition
              if (iphot.gt.0) then

                ! If column > row
                if (zupperJ) then

                  ! Spontaneous b-f transfer rate
                  call tEP(rJ,rJJ,Atom%phot(iphot)%TEI(iz),tEPcoeff)

                  ! Stimulated b-f transfer rate
                  call tSP(rJ,rJJ,JP(iphot,2),tSPcoeff)

                  ! Up-down b-f transfer rate
                  Scoeff(i,ii) = tEPcoeff + tSPcoeff + Scoeff(i,ii)

                ! If row > column
                else if (zlowerJ) then

                  ! Absorption b-f transfer rate
                  call tAP(rJ,rJJ,JP(iphot,1),tAPcoeff)

                  ! Down-up b-f transfer rate
                  Scoeff(i,ii)= Scoeff(i,ii) + tAPcoeff

                end if ! direction of transition
              end if ! There is a b-f transition
            end if ! populations
          end if ! Relaxation rate
        end if ! Level diagonal element
      end if ! Filled element

      ! Check if the row is not empty
      if (zflag.and.(abs(Scoeff(i,ii)).gt..0D0)) zflag=.False.

      ! Uncomment to avoid conjugation
     !zconj=.false.
     !zconj1=.false.

     !! If any of the two J pairs requires the use of conjugation
     !if (zconj.or.zconj1) then
      ! If the row J pairs requires the use of conjugation
      if (zconj) then

        Scoeff(i1,ii1) =  tJJsgn*Scoeff(i,ii)
        Tcoeff(i1,ii1) = -tJJsgn*Tcoeff(i,ii)

        ! Flag this element so we do not go over it again
        lfill(i1,ii1) = .False.

      end if ! If conjugation


                          !
                          ! Recover the identation
                          !

                        end do ! Q'
                      end do ! K'
                    end do ! J''' (column_b)
                  end do ! J'' (column_a)
                end do ! terms (column)

                ! If the row is empty, oportunity to correct it by
                ! checking that is the last ionization stage
                if (zflag) then

                  ! Check that it is the last one
                  if (iterm.eq.Atom%nMulti) then

                    ! Check that it is not a population
                    if (K.ne.0) then

                      ! Check that is a lone term in the ion
                      if (Atom%stage(iterm).ne. &
                          Atom%stage(iterm-1)) then

                        ! It is a lone ion and not rho00, flag this
                        ! row to receive the correct equation
                        ! rhoKQ = 0
                        rho(i) = -1d0

                        ! And unflag it
                        zflag = .False.

                      end if
                    end if
                  end if

                end if ! Empty row

                ! If the row is empty
                if (zflag) then

                  ! If we are zeroing out the ion
                  if (Atom%zero_ion) then
                    ! And this is from the last ion
                    if (Atom%stage(iterm).eq.Atom%stage(Atom%nMulti)) then

                      ! Then do not worry and just continue
                      cycle

                    end if ! Last stage
                  end if ! zeroing out the ion

                  write(umsg,*) ' # Element',iterm,real(rJ), &
                                real(rJ1),K,iQ,' isolated'
                  call abortedS(umsg,urou,tid,.True.,.True.)

                  return

                end if ! Empty row

              end do ! Q
            end do ! K
          end do ! J' (row_b)
        end do ! J (row_a)
      end do ! Terms (row)


      !
      ! Check that we have an square SE
      !
      if (i.ne.ii) then
        umsg = 'STcoeff is not square'
        call abortedS(umsg,urou,tid,.True.,.True.)
        return
      end if


      !
      ! Finish building the SEE with conjugation properties
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! For each level (row_a)
        do iJ=1,Atom%nJ(iterm)

          ! Get angular momentum
          rJ = Atom%rJval(iJ,iterm)

          ! Get level index
          ir = Atom%irho(iterm)%irho_ij(iJ)

          ! For each level (row_b)
          do iJ1=1,Atom%nJ(iterm)

            ! Get angular momentum
            rJ1 = Atom%rJval(iJ1,iterm)

            ! For each K
            do K=nint(abs(rJ-rJ1)), &
                 min(nint(rJ+rJ1),Atom%Kcut(iterm))

              ! For each Q
              do iQ=-K,K

                ! Get SEE index
                i = Atom%irho(iterm)%Jrho(iJ1,iJ)%kq(iQ,K)

                ! Check if special equation
                if (rho(i).lt.0d0) then

                  STcoeff(i,i) = 1d0
                  rho(i) = 0d0

                else

                  ! Check if Q<0 OR (Q==0 and J!=J')
                  zreal = (iQ.lt.0).or.((iQ.eq.0).and.(iJ.le.iJ1))

                  ! For each term (column)
                  do itterm=1,Atom%nMulti

                    ! For each level (column_a)
                    do iJJ=1,Atom%nJ(itterm)

                      ! Get angular momentum
                      rJJ = Atom%rJval(iJJ,itterm)

                      ! Get level index
                      iir = Atom%irho(itterm)%irho_ij(iJJ)

                      ! For each level (column_b)
                      do iJJ1=1,Atom%nJ(itterm)

                        ! Get angular momentum
                        rJJ1 = Atom%rJval(iJJ1,itterm)

                        ! Get sign given J
                        rJJsgn = Flgsg%sg(nint(rJJ-rJJ1))

                        ! Check relations between J indexes
                        zJlt = iJJ.lt.iJJ1
                        zJgt = iJJ.gt.iJJ1
                        zJeq = iJJ.eq.iJJ1

                        ! For each K'
                        do KK=nint(abs(rJJ-rJJ1)), &
                              min(nint(rJJ+rJJ1), Atom%Kcut(itterm))

                          ! For each Q'
                          do iQQ=-KK,KK

                            ! Get the SEE indexes
                            ii = Atom%irho(itterm)% &
                                      Jrho(iJJ1,iJJ)%kq(iQQ,KK)
                            ii1 = Atom%irho(itterm)% &
                                       Jrho(iJJ,iJJ1)%kq(-iQQ,KK)

                            ! Check if Q'<0 OR (Q'==0 and J''!=J''')
                            zzreal = (iQQ.lt.0).or.((iQQ.eq.0).and. &
                                     (iJJ.le.iJJ1))

                            ! Get sign including Q'
                            tJJsgn = rJJsgn*Flgsg%sg(iQQ)

                            ! If the flag is true for the row
                            if (zreal) then

                              ! If the flag is true for the
                              ! row AND column
                              if (zzreal) then

                                STcoeff(i,ii) = Scoeff(i,ii) + &
                                                tJJsgn*Scoeff(i,ii1)

                              ! If the flag is true for the row
                              ! and false for the column
                              else

                                STcoeff(i,ii) = tJJsgn*Tcoeff(i,ii1) &
                                                - Tcoeff(i,ii)

                              end if ! column flag

                            ! If the flag is false for the row
                            else

                              ! If the flag is false for the row and
                              ! true for the column
                              if (zzreal) then

                                STcoeff(i,ii) = tJJsgn*Tcoeff(i,ii1) &
                                                + Tcoeff(i,ii)

                              ! If the flag is false for the
                              ! row AND column
                              else

                                STcoeff(i,ii) = Scoeff(i,ii) - &
                                                tJJsgn*Scoeff(i,ii1)

                              end if ! column flag

                            end if ! row flag

                            ! If Q'==0 and J''==J''', it was accounted
                            ! twice
                            if (iQQ.eq.0.and.zJeq) &
                              STcoeff(i,ii) = .5d0*STcoeff(i,ii)

                          end do ! Q'
                        end do ! K '
                      end do ! J''' (column_b)
                    end do ! J'' (column_a)
                  end do ! terms (column)

                end if ! Flagged

              end do ! Q
            end do ! K
          end do ! J' (row_b)
        end do ! J (row_a)
      end do ! terms (row)

      end subroutine SEbuild

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes independent vector and changes the SEE matrix if
      !! requested\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !! STcoeff(dfloat(:,:)): Statistical equilibrium equations\n
      !!       rho(dfloat(:)): Array to store the solution of the
      !!                       statistical equilibrium equations\n
      !!          iz(integer): Height index
      subroutine initrho(Atom,STcoeff,rho,iz)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      double precision, dimension(:), intent(inout):: rho
      double precision, dimension(:,:), intent(inout):: STcoeff
      integer, intent(in):: iz

      ! Local
      integer:: it,iJ,iR,lstage,liR,ciR
      double precision:: rJJ


      ! Routine name
      urou = 'initrho'

      ! Initialize independent vector
      rho = 0d0

      ! If fixing populations
      if (Atom%fixp) then

        ! For each term
        do it=1,Atom%nMulti

          ! For each level (J)
          do iJ=1,Atom%nJ(it)

            ! Get SEE index of level i K=Q=0
            iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

            ! Make the row zero
            STcoeff(iR,:) = 0d0
            ! Except the diagonal
            STcoeff(iR,iR) = 1d0

            ! And put the current rho00 in the independent term
            rho(iR) = dble(Atom%crho(iR,iz))

          end do ! levels J
        end do ! terms

      ! If fixing only the lower term
      elseif (Atom%fixplt) then

        ! If multi-level atom
        if (Atom%ML) then

          ! Just fix lower level
          STcoeff(1,:) = 0d0
          STcoeff(1,1) = 1d0
          rho(1) = dble(Atom%crho(1,iz))

        ! If multi-term atom
        else

          ! For each sub-level
          do iJ=1,Atom%nJ(1)

            ! Get column for this sub-level
            iR = Atom%irho(1)%Jrho(iJ,iJ)%kq(0,0)

            ! Set ST element
            STcoeff(iR,:) = 0d0
            STcoeff(iR,iR) = 1d0

            ! Independent term
            rho(iR) = dble(Atom%crho(iR,iz))

          end do ! Sub-levels in ground-term

        end if ! Multi-level or multi-term

      ! Not fixing
      else

        ! Trace in ground level == 1
        rho(1) = 1d0

        ! Reset first row
        STcoeff(1,:) = 0d0

        ! For each term column
        do it=1,Atom%nMulti

          ! For each sub-level column
          do iJ=1,Atom%nJ(it)

            ! Get J value
            rJJ = Atom%rJval(iJ,it)

            ! Get population index
            iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

            ! Put the weight in the first row
            STcoeff(1,iR) = sqrt(2d0*rJJ+1d0)

          end do ! sub-level column
        end do ! term column

      end if ! Type of population fixing

      ! If zero_ion
      if (Atom%zero_ion) then

        ! Get last stage
        lstage = Atom%stage(Atom%nMulti)

        ! Get largest iR
        liR = Atom%ndim

        ! For every term
        do it=Atom%nMulti,1,-1

          ! If below the stage, done
          if (Atom%stage(it).lt.lstage) exit

          ! Get current iR
          ciR = Atom%irho(it)%Jrho(1,1)%kq(0,0)

          ! For every iR, nullify the ion
          do iR=ciR,lIR

            ! Make the row zero
            STcoeff(iR,:) = 0d0
            ! Except the diagonal
            STcoeff(iR,iR) = 1d0

          end do ! Term variables

          ! And update last iR
          liR = cIR - 1

        end do ! Terms

      end if

      return

      end subroutine initrho

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the statistical equilibrium equations\n
      !!        rho(dfloat(:)): Array to store the solution of the
      !!                        statistical equilibrium equations\n
      !!         ndim(integer): Size of the statistical
      !!                        equilibrium equations system\n
      !!  STcoeff(dfloat(:,:)): Statistical equilibrium equations
      subroutine densmatr(rho,ndim,STcoeff)

      ! I/O

      integer, intent(in):: ndim
      double precision, dimension(:,:), intent(in):: STcoeff
      double precision, dimension(:), intent(inout):: rho

      ! Local

      integer:: i
      integer, dimension(ndim):: indx


      !
      ! Solve SEE
      !
      call DGESV(ndim,1,STcoeff,ndim,indx,rho,ndim,i)

      end subroutine densmatr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Assigns the solution of the statistical equilibrium
      !! equations to the correct indexed variable.\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!       rho(dfloat(:)): Array to store the solution of the
      !!                       statistical equilibrium equations\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!               iz(integer): Height index\n
      !!              tid(integer): Thread index
      subroutine rhosol(Atom,rho,Flgsg,iz,tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: iz,tid
      double precision, dimension(:), intent(in):: rho

      ! Local

      logical:: zreal

      integer:: i,j,it,iJ,iJ1,K,iQ,iR,iI,iR0,iR1

      double precision:: rJ,rJ1,tsgn,rho0,rho0i


      ! Routine name
      urou = 'rhosol'


      !
      ! Initialize null flag
      !
      Atom%rhonull(:,iz) = .False.

      !
      ! Rearrange solution of SEE into rhoKQ array
      !

      ! For each term
      do it=1,Atom%nMulti

        ! For each level (J)
        do iJ=1,Atom%nJ(it)

          ! Get angular momentum
          rJ = Atom%rJval(iJ,it)

          ! Get level index
          i = Atom%irho(it)%irho_ij(iJ)

          ! For each level (J')
          do iJ1=1,Atom%nJ(it)
         !do iJ1=iJ,iJ !no J,J'

            ! Get angular momentum
            rJ1 = Atom%rJval(iJ1,it)

            ! Get level index
            j = Atom%irho(it)%irho_ij(iJ1)

            !
            ! Get the value of rho00 in order to check if rhoKQ are
            ! small quantities
            !

            ! If J!=J'
            if(i.ne.j)then

              ! Get SEE index of level i
              iR0 = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! Get SEE index of level j
              iR1 = Atom%irho(it)%Jrho(iJ1,iJ1)%kq(0,0)

              ! Check NaN
              if (isnan(rho(iR0)).or. &
                  isnan(rho(iR1))) then

                ! Not aborted
                write(umsg,'(A,i4,3(",",i4),A,2(1x,es11.4))') &
                  'NaN in SEE solution'//new_line('A')// &
                  '(iz,it,iJ,iJ1)=(',iz,it,iJ,iJ1,')'// &
                  new_line('A')//'rho00: ',rho(iR0),rho(iR1)

                call abortedS(umsg,urou,tid,.True.,.True.)

              end if

              ! Check negativity
              if (rho(iR0).lt.0d0.or.rho(iR1).lt.0d0) then

                write(umsg,'(A,i4,3(",",i4),A,2(1x,es11.4))') &
                  'Negative population in SEE solution'// &
                   new_line('A')// &
                  '(iz,it,iJ,iJ1)=(',iz,it,iJ,iJ1,')'// &
                  new_line('A')//'rho00: ',rho(iR0),rho(iR1)

                call abortedS(umsg,urou,tid,.not.nphysR,.True.)

              end if

              ! Geometric average of the two rho00
              rho0 = sqrt(abs(rho(iR0)*rho(iR1)))

              ! Inverse
              if (rho0.gt.0) then
                rho0i = 1d0/rho0
              else
                rho0i = 0d0
              end if

            ! If J==J'
            else

              ! Get SEE index
              iR0 = Atom%irho(it)%Jrho(iJ1,iJ)%kq(0,0)

              ! Check NaN
              if (isnan(rho(iR0))) then

                write(umsg,'(A,i4,2(",",i4),A,1x,es11.4)') &
                  'NaN in SEE solution'//new_line('A')// &
                  '(iz,it,iJ)=(',iz,it,iJ,')'// &
                  new_line('A')//'rho00: ',rho(iR0)

                call abortedS(umsg,urou,tid,.True.,.True.)

              end if

              ! Check negativity
              if (rho(iR0).lt.0d0) then

                write(umsg,'(A,i4,2(",",i4),A,1x,es11.4)') &
                  'Negative population in SEE solution'// &
                  new_line('A')// &
                  '(iz,it,iJ)=(',iz,it,iJ,')'// &
                  new_line('A')//'rho00: ',rho(iR0)

                call abortedS(umsg,urou,tid,.not.nphysR,.True.)

              end if

              ! rho00
              rho0 = rho(iR0)

              ! Inverse
              if (rho0.gt.0) then
                rho0i = 1d0/rho0
              else
                rho0i = 0d0
              end if

            endif ! Compare J and J'

            ! For each K
            do K=nint(abs(rJ-rJ1)), &
                 min(nint(rJ+rJ1),Atom%Kcut(it))

              ! For each Q
              do iQ=-K,K

                ! Get indexes of the real and imaginary parts in rho
                iR = Atom%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)
                iI = Atom%irho(it)%Jrho(iJ,iJ1)%kq(-iQ,K)

                ! If Q<0 or (Q'==0 and iJ<=iJ')
                zreal = (iQ.lt.0).or.((iQ.eq.0).and.(iJ.le.iJ1))

                ! Get the sign from J, J' and Q
                tsgn = Flgsg%sg(nint(rJ-rJ1)-iQ)

                ! Check NaN
                if (isnan(rho(iR)).or.isnan(rho(iI))) then

                  write(umsg,'(A,i4,3(",",i4),",",i1,",",i2,'// &
                             'A,2(1x,es11.4))') &
                    'NaN in SEE solution'//new_line('A')// &
                    '(iz,it,iJ,iJ1,K,Q)=(',iz,it,iJ,iJ1,K,iQ,')'// &
                    new_line('A')//'rhoKQ: ',rho(iR),rho(iI)

                  call abortedS(umsg,urou,tid,.True.,.True.)

                end if

                ! Check magnitude
                if (abs(rho(iR)).gt.rho0.or.abs(rho(iI)).gt.rho0) then

                  ! If J=J'
                  if (i.eq.j) then

                    write(umsg,'(A,i4,3(",",i4),",",i1,",",i2,'// &
                               'A,1x,es11.4,A,2(1x,es11.4))') &
                      'Atomic polarization larger than '// &
                      'population in SEE solution'// &
                      new_line('A')// &
                      '(iz,it,iJ,iJ1,K,Q)=(',iz,it,iJ,iJ1,K,iQ,')'// &
                      new_line('A')//'rho00: ',rho0, &
                      new_line('A')//'rhoKQ: ',rho(iR),rho(iI)

                  ! If J!=J'
                  else

                    write(umsg,'(A,i4,3(",",i4),",",i1,",",i2,'// &
                               '3(A,1x,es11.4),A,2(1x,es11.4))') &
                      'Atomic polarization larger than '// &
                      'population in SEE solution'// &
                      new_line('A')// &
                      '(iz,it,iJ,iJ1,K,Q)=(',iz,it,iJ,iJ1,K,iQ,')'// &
                      new_line('A')//'rho00: ',rho0, &
                      '= sqrt(',rho(iR0),'*',rho(iR1),')'// &
                      new_line('A')//'rhoKQ: ',rho(iR),rho(iI)

                  end if ! J diagonality

                  call abortedS(umsg,urou,tid,.not.nphysR,.True.)

                end if

                ! If Q<0 or (Q'==0 and iJ<=iJ')
                if (zreal) then

                  Atom%crho(iR,iz) = dcmplx(rho(iR), &
                                            -tsgn*rho(iI))

                ! If Q>=0 and (Q'!=0 or iJ>iJ')
                else

                  Atom%crho(iR,iz) = dcmplx(tsgn*rho(iI), &
                                            rho(iR))

                end if

                ! If Q==0 and diagonal, no imaginary part
                if (iQ.eq.0.and.iJ.eq.iJ1) &
                   Atom%crho(iR,iz) = dcmplx(rho(iR),.0d0)

                !
                ! Test if rhoKQ is significant
                !

                ! Check value relative to population
                if (abs(Atom%crho(iR,iz)*rho0i).lt.TINYR) &
                  Atom%rhonull(iR,iz) = .True.

              end do ! Q
            end do ! K
          end do ! levels J'
        end do ! levels J
      end do ! terms

      end subroutine rhosol

!#####################################################################
!#####################################################################
!#####################################################################

      end module see_mod
