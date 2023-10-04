      !> Populations statistical equilibrium equations
      module seei_mod
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
!     07/03/2023 V3.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     07/03/2023:    V3.0.4 - Renamed SEEI to SEEI_actual and added
!                             a new SEEI that calls it. Now, if
!                             populations are negative and you are
!                             running with ALI, the SEEI try again
!                             switching it off (TdPA)
!                           - Added an alternative trace equation
!                             which fixes the lower term population
!                             but not ensures the particle
!                             conservation (TdPA)
!
!     11/24/2022:    V3.0.3 - Removed non-used variables (TdPA)
!
!     11/10/2022:    V3.0.2 - Added the option to zero out the
!                             populations of the last ion of the
!                             atom. This is useful when using
!                             two-level atoms while keeping the
!                             ionization stage to compute the
!                             broadening parameters (TdPA)
!                           - Added the initrhoI routine to
!                             achieve the latter (TdPA)
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
!     11/13/2019:    V1.3.7 - The Atom%popu variable is updated with
!                             the normalized population (TdPA)
!
!     10/18/2019:    V1.3.6 - Skip if the atom ask to keep the
!                             populations fixed (TdPA)
!
!     08/19/2019:    V1.3.5 - Bugfix: In 1.3.4, when the checking was
!                             changed to level indexes, they were put
!                             after a conditional that needed them
!                             defined already (TdPA)
!
!     08/08/2019:    V1.3.4 - Removed check on L values for permitted
!                             radiative transition (TdPA)
!                           - Checking zupper, zlower, and zrelax on
!                             the level indexes directly (TdPA)
!
!     06/04/2019:    V1.3.3 - Bugfix: The Lambda contributions for
!                             the RA associated coefficients had the
!                             wrong sign (TdPA)
!
!     05/31/2019:    V1.3.2 - LamL and LamP do not need to have
!                             declared sizes (TdPA)
!
!     03/18/2019:    V1.3.1 - Physical errors do not make rhosol to
!                             return (TdPA)
!                           - Can skip ALI calculations (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     09/08/2017:    V1.0.2 - Added checks to avoid non physical
!                             solutions (TdPA)
!
!     06/14/2017:    V1.0.1 - In densmatrI, rho is intent(inout), to
!                             avoid valgrind complains (TdPA)
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
!  SEEI:
!    Manage the solution of the statistical equilibrium equations by
!  calling the actual manager with/out ALI
!
!  SEEI_actual:
!    Manage the solution of the statistical equilibrium equations
!
!  SEbuildI:
!    This subroutine builds the coefficient matrix for the system
!  of S.E. equations, for the multi-level atom
!
!  ALIbuildI:
!    This subroutine builds the lambda operator coefficient matrix
!  for the system of S.E. equations, for the multi-level atom
!
!  initrho:
!    Initializes rho vector or considers zero ion
!
!  densmatrI:
!    This subroutine calculates the density matrix solution
!  of the SE equations
!
!  rhosolI:
!    Restore the solution density matrix
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : k2f
      use seeiaux_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calls subroutines to solve the statistical equilibrium
      !! equations for the populations.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atom0(Rhoc_class): A copy of Atom\n
      !!           JRad(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!          JRadS(dfloat(:)): Mean intensity integrated over
      !!                            emission profile\n
      !!        JPhot(dfloat(:,:)): Mean intensity integrated with
      !!                            photoionizations\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!         LamL(dfloat(:,:)): Lambda operator for bound-bound
      !!                            transitions\n
      !!       LamP(dfloat(:,:,:)): Lambda operator for bound-free
      !!                            transitions\n
      !!               iz(integer): Height index\n
      !!             lALI(logical): Apply ALI\n
      !!              tid(integer): thread index
      subroutine SEEI(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP,iz,lALI, &
                      tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(inout):: Atom0
      logical, intent(in):: lALI
      integer, intent(in):: iz,tid
      double precision, dimension(:), intent(in):: JRad,JRadS
      double precision, dimension(:,:), intent(in):: LamL
      double precision, dimension(:,:),intent(in):: Jphot
      double precision, dimension(:,:,:),intent(in):: LamP

      ! Local
      logical:: try_no_ALI


      ! If fixing populations, do not do anything
      if (Atom%fixp) return

      !
      ! Call real SEEI with current ALI configuration
      !
      call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                       iz,lALI,try_no_ALI,tid)

      !
      ! If we need to try without ALI now because there was a
      ! negative population, do it now
      !
      if (try_no_ALI) then

        call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                         iz,.False.,try_no_ALI,tid)

      end if

      end subroutine SEEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calls subroutines to solve the statistical equilibrium
      !! equations for the populations.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atom0(Rhoc_class): A copy of Atom\n
      !!           JRad(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!          JRadS(dfloat(:)): Mean intensity integrated over
      !!                            emission profile\n
      !!        JPhot(dfloat(:,:)): Mean intensity integrated with
      !!                            photoionizations\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!         LamL(dfloat(:,:)): Lambda operator for bound-bound
      !!                            transitions\n
      !!       LamP(dfloat(:,:,:)): Lambda operator for bound-free
      !!                            transitions\n
      !!               iz(integer): Height index\n
      !!             lALI(logical): Apply ALI\n
      !!       try_no_ALI(logical): Output signal to try without ALI\n
      !!              tid(integer): thread index
      subroutine SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                             iz,lALI,try_no_ALI,tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(inout):: Atom0
      logical, intent(in):: lALI
      logical, intent(out):: try_no_ALI
      integer, intent(in):: iz,tid
      double precision, dimension(:), intent(in):: JRad,JRadS
      double precision, dimension(:,:), intent(in):: LamL
      double precision, dimension(:,:),intent(in):: Jphot
      double precision, dimension(:,:,:),intent(in):: LamP

      ! Local

      double precision:: rho(Atom%nlevel)
      double precision:: STcoeff(Atom%nlevel,Atom%nlevel)


      !
      ! Build SEE equations for multi-level atom
      !
      call SEbuildI(Atom,JRad,JRadS,Jphot,STcoeff,iz,tid)


      !
      ! Add the contributions due to Lambda operator
      !
      if (lALI) call ALIbuildI(Atom,Atom0,LamL,LamP,STcoeff,iz)


      !
      ! Initialize rho and fix populations if requested
      !
      call initrhoI(Atom,Atom0,STcoeff,rho,iz)


      !
      ! Solve the SEE
      !
      call densmatrI(rho,Atom%nlevel,STcoeff)


      !
      ! Rearrange the solution into the rhoKQ matrices
      !
      call rhosolI(Atom,rho,iz,lALI,try_no_ALI,tid)


      end subroutine SEEI_actual

!#####################################################################
!#####################################################################
!#####################################################################

      !> Builds the statistical equilibrium equations system for the
      !! populations\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           RadJ(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!          RadJS(dfloat(:)): Mean intensity integrated over
      !!                            emission profile\n
      !!             JP(dfloat(:)): Mean intensity integrated with
      !!                            photoionizations\n
      !!      STcoeff(dfloat(:,:)): Statistical equilibrium equations
      !!                            system\n
      !!               iz(integer): Height index\n
      !!              tid(integer): thread index
      subroutine SEbuildI(Atom,RadJ,RadJS,JP,STcoeff,iz,tid)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz,tid
      double precision,dimension(:),intent(in):: RadJ,RadJS
      double precision,dimension(:,:),intent(in):: JP
      double precision, dimension(:,:), intent(out):: STcoeff

      ! Local

      logical:: zpermit,zpermitJ,zPpermit,zRpermitJ
      logical:: zrelax,zupper,zlower,zrelaxJ,zupperJ,zlowerJ
      logical:: zJ,zflag

      integer:: itran,ftran,fftran,iphot,iterm,itterm,iJ,iJJ,i,ii

      double precision:: rL,S,rJ,rLL,SS,rJJ
      double precision:: rSCcoeff,rACcoeff
      double precision:: tSCcoeff,tACcoeff
      double precision:: rEcoeff,tEcoeff
      double precision:: rScoeff,rAcoeff
      double precision:: tScoeff,tAcoeff
      double precision:: rEPcoeff,tEPcoeff
      double precision:: rSPcoeff,rAPcoeff
      double precision:: tSPcoeff,tAPcoeff


      ! Routine name
      urou = 'SEbuildI'


      !
      ! Initialize SEE coefficients
      !
      STcoeff = 0d0


      !
      ! Build SEE
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! Get term quantities
        rL = Atom%rLval(iterm)
        S = Atom%Sval(iterm)

        ! Spontaneous relaxation b-b rate
        call rEI(iterm,Atom%Ecoeff,rEcoeff)

        ! For each level (row)
        do iJ=1,Atom%nJ(iterm)

          ! Get level momentum
          rJ = Atom%rJval(iJ,iterm)

          ! Get level index
          i = Atom%irho(iterm)%irho_ij(iJ)

          ! Spontaneous relaxation b-f rate
          call rEPI(i,Atom%iphot(:,i),Atom%phot,iz,rEPcoeff)

          ! Flag of filled line
          zflag = .True.

          ! For each term (column)
          do itterm=1,Atom%nMulti

            ! Check how the column relates to the row
            zrelax = itterm.eq.iterm
            zupper = itterm.gt.iterm
            zlower = itterm.lt.iterm

            ! Get term quantities
            rLL = Atom%rLval(itterm)
            SS = Atom%Sval(itterm)

            ! Check if there is a permitted transition between
            ! these terms
            itran = Atom%irad(itterm,iterm)
            zpermit = itran.gt.0

            ! For each level (column)
            do iJJ=1,Atom%nJ(itterm)

              ! Get level quantity
              rJJ = Atom%rJval(iJJ,itterm)

              ! Get level index
              ii = Atom%irho(itterm)%irho_ij(iJJ)

              ! If there is a permitted transition between terms,
              ! check that there is a permitted transition between
              ! the levels
              if (itran.gt.0) then

                if (zupper) then
                  ftran = Atom%fst(itran)%irad(iJJ,iJ)
                else
                  ftran = Atom%fst(itran)%irad(iJ,iJJ)
                end if
                if (ftran.gt.0) fftran = Atom%ifst_ij(ftran,itran)

              else

                ftran = 0
                fftran = 0

              end if ! permitted term b-b transition

              ! Check if there is a b-f transition between the levels
              iphot = Atom%iphot(ii,i)
              zPpermit = iphot.gt.0

              ! Flag for permitted FS transition
              zRpermitJ = zpermit.and.ftran.gt.0

              ! Flag for permitted at all
              zpermitJ = zRpermitJ.or.zPpermit

              ! Same J for the two levels
              zJ = iJ.eq.iJJ

              ! If column > row and it is a b-b transition,
              ! spontaneous b-b transfer rate
              if (zupper.and.zRpermitJ) &
                call tEI(rJ,rJJ,Atom%fst(itran)%Aul(iJJ,iJ),tEcoeff)

              ! If column > row and it is a b-f transition,
              ! spontaneous b-f transfer rate
              if (zupper.and.zPpermit) &
                call tEPI(rJ,rJJ,Atom%phot(iphot)%TEI(iz),tEPcoeff)

      !
      ! Reset the Identation
      !

      ! If there is a transition between levels
      if (zpermitJ) then

        ! If column > row
        if (zupper) then

          ! If it is a b-b transition
          if(zRpermitJ)then

            ! Contribution to spontaneous transfer rate
            STcoeff(i,ii) = tEcoeff

            ! If there is stimulated emission
            if (stm) then

              ! Stimulated b-b transfer rate
              call tSI(rJ,rJJ,Atom%fst(itran)%Blu(iJ,iJJ), &
                       RadJS(fftran),tScoeff)

              STcoeff(i,ii) = tScoeff + STcoeff(i,ii)

            end if ! Stimulated emission

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! Stimulated b-f transfer rate
            call tSPI(rJ,rJJ,JP(iphot,2),tSPcoeff)

            ! Up-down b-f transfer rate
            STcoeff(i,ii) = tEPcoeff + tSPcoeff + STcoeff(i,ii)

          end if ! b-f transition

        end if ! column > row

        ! If row > column
        if (zlower) then

          ! If it is a b-b transition
          if (zRpermitJ) then

            ! Absorption b-b transfer rate
            call tAI(rJ,rJJ,Atom%fst(itran)%Blu(iJJ,iJ), &
                     RadJ(fftran),tAcoeff)

            STcoeff(i,ii) = tAcoeff

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! Absorption b-f transfer rate
            call tAPI(rJ,rJJ,JP(iphot,1),tAPcoeff)

            STcoeff(i,ii) = STcoeff(i,ii) + tAPcoeff

          end if ! b-f transition

        end if ! row > column

      end if ! There is a transition


      ! Diagonal element
      if (zrelax) then

        if (zJ) then

          ! Spontaneous relaxation rates
          STcoeff(i,ii) = -rEcoeff - rEPcoeff

          ! If stimulated emission
          if (stm) then

            ! Stimulated b-b relaxation rate
            call rSI(iterm,iJ,Atom%irad(:,iterm),Atom%fst, &
                     Atom%ifst_ij,Atom%nJ,Atom%rJval,rJ, &
                     RadJS,rScoeff)

          else

            rScoeff = 0d0

          end if ! Stimulated emission

          ! Absorption b-b relaxation rate
          call rAI(iterm,iJ,Atom%irad(:,iterm),Atom%fst, &
                   Atom%ifst_ij,Atom%nJ,Atom%nMulti,RadJ,rAcoeff)

          STcoeff(i,ii) = -rScoeff - rAcoeff + STcoeff(i,ii)

          ! Stimulated b-f relaxation rate
          call rSPI(i,Atom%iphot(:,i),JP(:,2),rSPcoeff)

          ! Absorption b-f relaxation rate
          call rAPI(i,Atom%iphot(:,i),Atom%nlevel, &
                    JP(:,1),rAPcoeff)

          STcoeff(i,ii)= -rSPcoeff - rAPcoeff + STcoeff(i,ii)

        end if ! Same level

      end if ! Same term

      ! Collisions
      if (maxval(Atom%CcoeffJ(i,:,iz)).gt.0d0.or. &
          Atom%CcoeffJ(ii,i,iz).gt.0d0) then

        ! Check relative positions of levels
        zupperJ = ii.gt.i
        zlowerJ = ii.lt.i
        zrelaxJ = ii.eq.i

        ! column > row
        if (zupperJ) then

          ! Colisional de-excitation transfer rate
          call tSFCI(rJ,rJJ,Atom%CcoeffJ(ii,i,iz),tSCcoeff)

          STcoeff(i,ii) = tSCcoeff + STcoeff(i,ii)

        ! row > column
        else if (zlowerJ) then

          ! Colisional excitation transfer rate
          call tAFCI(rJ,rJJ,Atom%CcoeffJ(ii,i,iz),tACcoeff)

          STcoeff(i,ii) = STcoeff(i,ii) + tACcoeff

        ! diagonal element
        else if (zrelaxJ) then

          ! Colisional de-excitation relaxation rate
          call rSFCI(i,Atom%CcoeffJ(i,:,iz),rSCcoeff)

          ! Colisional excitation relaxation rate
          call rAFCI(i,Atom%nlevel,Atom%CcoeffJ(i,:,iz),rACcoeff)

          STcoeff(i,ii) = -rSCcoeff - rACcoeff + STcoeff(i,ii)

        end if ! Relation between levels
      end if ! There are collisions

      ! Check if there is any non-zero coefficients in this row
      if (zflag.and.(abs(STcoeff(i,ii)).gt..0D0)) zflag=.false.

              !
              ! Recover the identation
              !

            end do ! iJJ
          end do ! itterm

          ! If the matrix is singular
          if (zflag) then

            ! If we are zeroing out the ion
            if (Atom%zero_ion) then
              ! And this is from the last ion
              if (Atom%stage(iterm).eq.Atom%stage(Atom%nMulti)) then

                ! Then do not worry and just continue
                cycle

              end if ! Last stage
            end if ! zeroing out the ion


            ! If not aborted yet
            write(umsg,*) ' # Element',iterm,real(rJ), &
                          ' isolated'
            call abortedS(umsg,urou,tid,.True.,.True.)

            return

          end if ! Singular matrix

        end do ! iJ
      end do ! iterm

      ! Check the matrix is actually square
      if (i.ne.ii) then
        umsg = 'STcoeff is not square'
        call abortedS(umsg,urou,tid,.True.,.True.)
        return
      end if

      return

      end subroutine SEbuildI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds the ALI terms to the statistical equilibrium equations
      !! system for the populations.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atom0(Rhoc_class): A copy of Atom\n
      !!         LamL(dfloat(:,:)): Lambda operator for bound-bound
      !!                            transitions\n
      !!       LamP(dfloat(:,:,:)): Lambda operator for bound-free
      !!                            transitions\n
      !!      STcoeff(dfloat(:,:)): Statistical equilibrium equations
      !!                            system\n
      !!               iz(integer): Height index
      subroutine ALIbuildI(Atom,Atom0,LamL,LamP,STcoeff,iz)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Rhoc_class), intent(inout):: Atom0
      integer, intent(in):: iz
      double precision,dimension(:,:),intent(in):: LamL
      double precision,dimension(:,:,:),intent(in):: LamP
      double precision, dimension(:,:), intent(inout):: STcoeff

      ! Local

      logical:: zpermit,zpermitJ,zPpermit,zRpermitJ
      logical:: zrelax,zupper,zlower,zJ

      integer:: itran,ftran,fftran,iphot,iterm,itterm
      integer:: iJ,iJJ,i,ii,iR,iRR

      double precision:: rL,S,rJ,rLL,SS,rJJ
      double precision:: tAcoeff,tAPcoeff


      !
      ! Build ALI contribution to SEE
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! Get term quantities
        rL = Atom%rLval(iterm)
        S = Atom%Sval(iterm)

        ! For each level (row)
        do iJ=1,Atom%nJ(iterm)

          ! Get level quantity
          rJ = Atom%rJval(iJ,iterm)

          ! Get level index
          i = Atom%irho(iterm)%irho_ij(iJ)

          ! Get rhoKQ index
          iR = Atom%irho(iterm)%Jrho(iJ,iJ)%kq(0,0)

          ! For each term (column)
          do itterm=1,Atom%nMulti

            ! Get term quantities
            rLL = Atom%rLval(itterm)
            SS = Atom%Sval(itterm)

            itran = Atom%irad(itterm,iterm)

            ! Check if there is a permitted transition between
            ! these terms
            zpermit = itran.gt.0

            ! For each level (column)
            do iJJ=1,Atom%nJ(itterm)

              ! Get level quantity
              rJJ = Atom%rJval(iJJ,itterm)

              ! Get level index
              ii = Atom%irho(itterm)%irho_ij(iJJ)

              ! Get rhoKQ index
              iRR = Atom%irho(itterm)%Jrho(iJJ,iJJ)%kq(0,0)

              ! Check how the column relates to the row
              zrelax = i.eq.ii
              zupper = ii.gt.i
              zlower = ii.lt.i

              ! If there is a permitted transition between terms,
              ! check that there is a permitted transition between
              ! the levels
              if (itran.gt.0) then

                if (zupper) then
                  ftran = Atom%fst(itran)%irad(iJJ,iJ)
                else
                  ftran = Atom%fst(itran)%irad(iJ,iJJ)
                end if
                if(ftran.gt.0) fftran = Atom%ifst_ij(ftran,itran)

              else

                ftran = 0
                fftran = 0

              end if

              ! Check if there is a b-f transition between the levels
              iphot = Atom%iphot(ii,i)
              zPpermit = iphot.gt.0

              ! Flag for permitted FS transition
              zRpermitJ = zpermit.and.ftran.gt.0

              ! Flag for permitted at all
              zpermitJ = zRpermitJ.or.zPpermit

              ! Same J for the two levels
              zJ = iJ.eq.iJJ


      !
      ! Reset the Identation
      !

      ! First row is not modified, it is the trace
      if (i.eq.1) cycle

      ! If there is a transition between levels
      if (zpermitJ) then

        ! If column > row
        if (zupper) then

          ! If it is a b-b transition
          if (zRpermitJ) then

            ! ALI contribution is like absorption transfer rate
            call tAI(0d0,0d0,Atom%fst(itran)%Blu(iJ,iJJ), &
                     LamL(1,fftran),tAcoeff)

            STcoeff(i,ii) = STcoeff(i,ii) - &
                            tAcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! ALI contribution is like absorption transfer rate
            call tAPI(0d0,0d0,LamP(1,iphot,1),tAPcoeff)

            STcoeff(i,ii) = STcoeff(i,ii) - &
                           tAPcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAPcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-f transition

        end if ! column > row

        ! If row > column
        if (zlower) then

          ! If it is a b-b transition
          if (zRpermitJ) then

            ! ALI contribution is like absorption transfer rate
            call tAI(rJ,rJJ,Atom%fst(itran)%Blu(iJJ,iJ), &
                     LamL(1,fftran),tAcoeff)

            STcoeff(i,ii) = STcoeff(i,ii) - &
                            tAcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! ALI contribution is like absorption transfer rate
            call tAPI(rJ,rJJ,LamP(1,iphot,1),tAPcoeff)

            STcoeff(i,ii)= STcoeff(i,ii) - &
                           tAPcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAPcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-f transition
        end if ! row > column
      end if ! There is a transition

              !
              ! Recover identation
              !

            end do ! iJJ
          end do ! itterm
        end do ! iJ
      end do ! iterm

      return

      end subroutine ALIbuildI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initializes independent vector and changes the SEE matrix if
      !! requested\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atom0(Rhoc_class): A copy of Atom\n
      !! STcoeff(dfloat(:,:)): Statistical equilibrium equations\n
      !!       rho(dfloat(:)): Array to store the solution of the
      !!                       statistical equilibrium equations\n
      !!          iz(integer): Height index
      subroutine initrhoI(Atom,Atom0,STcoeff,rho,iz)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(inout):: Atom0
      integer, intent(in):: iz
      double precision, dimension(:), intent(inout):: rho
      double precision, dimension(:,:), intent(inout):: STcoeff

      ! Local
      integer:: it,iJ,i,iR,lstage
      double precision:: rJJ


      ! Routine name
      urou = 'initrhoI'

      ! Initialize independent vector
      rho = 0d0

      ! If fixing only the lower term
      if (Atom%fixplt) then

        ! If multi-level atom
        if (Atom%ML) then

          ! Just fix lower level
          STcoeff(1,:) = 0d0
          STcoeff(1,1) = 1d0
          rho(1) = dble(Atom0%crho(1,iz))

        ! If multi-term atom
        else

          ! For each sub-level
          do iJ=1,Atom%nJ(1)

            ! Get level and rhoKQ index
            i = Atom%irho(1)%irho_ij(iJ)
            iR = Atom%irho(1)%Jrho(iJ,iJ)%kq(0,0)

            ! Fix level
            STcoeff(i,:) = 0d0
            STcoeff(i,i) = 1d0
            rho(i) = dble(Atom0%crho(iR,iz))

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

            ! Get level index
            i = Atom%irho(it)%irho_ij(iJ)

            ! Get J value
            rJJ = Atom%rJval(iJ,it)

            ! Put the weight in the first row
            STcoeff(1,i) = sqrt(2d0*rJJ+1d0)

          end do ! sub-level column
        end do ! term column

      end if ! Type of population fixing


      ! If zero_ion
      if (Atom%zero_ion) then

        ! Get last stage
        lstage = Atom%stage(Atom%nMulti)

        ! For every term
        do it=Atom%nMulti,1,-1

          ! If below the stage, done
          if (Atom%stage(it).lt.lstage) exit

          ! For every level
          do iJ=1,Atom%nJ(it)

            ! Get level index
            i = Atom%irho(it)%irho_ij(iJ)
            iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

            ! Make the row zero
            STcoeff(i,:) = 0d0
            ! Except the diagonal
            STcoeff(i,i) = 1d0
            ! And keep population
            rho(i) = dble(Atom0%crho(iR,iz))

          end do ! Term-levels
        end do ! Term

      end if

      return

      end subroutine initrhoI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solves the statistical equilibrium equations for
      !! populations equations for the populations.\n
      !!         rho(dfloat(:)): Array to store the solution of the
      !!                         statistical equilibrium equations\n
      !!          ndim(integer): Size of the statistical
      !!                         equilibrium equations system\n
      !!   STcoeff(dfloat(:,:)): Statistical equilibrium equations
      !!                         system
      subroutine densmatrI(rho,ndim,STcoeff)

      ! I/O

      integer, intent(in):: ndim
      double precision, dimension(:), intent(out):: rho
      double precision, dimension(:,:), intent(in):: STcoeff

      ! Local

      integer:: i
      integer, dimension(ndim):: indx

      !
      ! Solve SEE
      !
      call DGESV(ndim,1,STcoeff,ndim,indx,rho,ndim,i)

      end subroutine densmatrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Assigns the solution of the statistical equilibrium
      !! equations to the correct indexed variable.\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!       rho(dfloat(:)): Solution of SEE\n
      !!               iz(integer): Height index\n
      !!             lALI(logical): Apply ALI\n
      !!       try_no_ALI(logical): Output signal to try without ALI\n
      !!              tid(integer): thread index
      subroutine rhosolI(Atom,rho,iz,lALI,try_no_ALI,tid)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      logical, intent(in):: lALI
      logical, intent(out):: try_no_ALI
      integer, intent(in):: iz,tid
      double precision, dimension(:), intent(in):: rho

      ! Local

      integer:: i, it, iJ, iR
      double precision:: rJ


      ! Routine name
      urou = 'rhosolI'

      ! Initialize check flag
      try_no_ALI = .False.


      !
      ! Initialize null flag
      !
      Atom%rhonull(:,iz) = .True.

      !
      ! Check negative populations
      !
      if (minval(rho).lt.0d0) then

        ! If doing ALI
        if (lALI) then

            write(umsg,'(A)') 'Negative population in SEEI '// &
                              'solution, will try without ALI'
            call abortedS(umsg,urou,tid,.False.,.True.)

            ! Flag and go back
            try_no_ALI = .True.
            return

        ! Not doing ALI
        else

          ! Find
          do i=1,Atom%nlevel

            ! Check negativity
            if(rho(i).lt.0d0)then

              write(umsg,'(A,i4,",",i4,A,1x,es11.4)') &
                'Negative population in SEE solution'// &
                new_line('A')// &
                '(iz,il)=(',iz,i,')'// &
                new_line('A')//'rho00: ',rho(i)

              call abortedS(umsg,urou,tid,.not.nphysR,.True.)

            end if ! Negative population at this heright

          end do ! Levels

          ! Go back, we are aborting
          return

        end if ! Doing ALI
      end if ! Negative population


      !
      ! Rearrange solution of SEE into rhoKQ array
      !

      ! For each term
      do it=1,Atom%nMulti

        ! For each level
        do iJ=1,Atom%nJ(it)

          ! Level index
          i = Atom%irho(it)%irho_ij(iJ)

          ! J value
          rJ = Atom%rJval(iJ,it)

          ! rhoKQ index
          iR = Atom%irho(it)%Jrho(iJ,iJ)%kq(0,0)

          ! Check NaN
          if(isnan(rho(i)))then

            write(umsg,'(A,i4,",",i4,A,1x,es11.4)') &
              'NaN in SEE solution'//new_line('A')// &
              '(iz,il)=(',iz,i,')'// &
              new_line('A')//'rho00: ',rho(i)

            call abortedS(umsg,urou,tid,.True.,.True.)

          end if

          ! Store the corresponding element
          Atom%crho(iR,iz) = dcmplx(rho(i),0d0)

          ! Store normalized population
          Atom%popu(i,iz) = sqrt(2d0*rJ + 1d0)*rho(i)

          ! And cancel the flag
          Atom%rhonull(iR,iz) = .False.

        end do ! levels
      end do ! terms

      end subroutine rhosolI

!#####################################################################
!#####################################################################
!#####################################################################

      end module seei_mod
