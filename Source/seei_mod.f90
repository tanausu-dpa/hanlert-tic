      !> Populations statistical equilibrium equations
      module seei_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     18/12/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/12/2025:    V4.0.3 - NaN checks now use ieee (TdPA)
!                           - Only notify of negative populations in
!                             SEE if in synthesis mode (TdPA)
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
!  SEEI
!    Manage the solution of the statistical equilibrium equations for
!  the atomic populations
!
!  SEEI_actual
!    Solve the statistical equilibrium equations for the atomic
!  populations
!
!  SEbuildI
!    Build the statistical equilibrium equations
!
!  ALIbuildI
!    Add the Accelerated Lambda Iteration terms to the statistical
!  equilibrium equations
!
!  initrhoI
!    Initialize independent vector and manage additional atomic
!  flags
!
!  densmatrI
!    Solve the statistical equilibrium equations
!
!  rhosolI
!    Transform the solution of the statistical equilibrium equations
!  into atomic populations
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

      !> Manage the solution of the statistical equilibrium equations
      !! for the atomic populations\n
      !!     Atom(Atom_class): Structure with atomic data\n
      !!    Atom0(Rhoc_class): Structure to store the density matrix
      !!                       of the previous iteration\n
      !!      JRad(double(:)): Mean intensity integrated over the
      !!                       absorption profile\n
      !!     JRadS(double(:)): Mean intensity integrated over the
      !!                       emission profile\n
      !!   JPhot(double(:,:)): Intensity integrals in the
      !!                       photoionization rates\n
      !!    LamL(double(:,:)): Lambda operator for bound-bound
      !!                       transitions\n
      !!  LamP(double(:,:,:)): Lambda operator for bound-free
      !!                       transitions\n
      !!          iz(integer): Height index\n
      !!        lALI(logical): If to apply ALI\n
      !!        ALIp(logical): If to apply ALI to bound-free\n
      !!       ALIao(logical): Switch off ALI if negative populations
#ifdef DEBUGSEE
      subroutine SEEI(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP,iz,lALI, &
                      ALIp,ALIao,INPUT)
#else
      subroutine SEEI(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP,iz, &
                      lALI,ALIao,ALIp)
#endif

      ! I/O

#ifdef DEBUGSEE
      type(Input_class), intent(in):: INPUT
#endif
      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(in):: Atom0
      logical, intent(in):: lALI,ALIp,ALIao
      integer, intent(in):: iz
      double precision, dimension(:), intent(in):: JRad,JRadS
      double precision, dimension(:,:), intent(in):: LamL
      double precision, dimension(:,:), intent(in):: Jphot
      double precision, dimension(:,:,:), intent(in):: LamP

      ! Local

      logical:: try_no_ALI


      ! If fixing populations, do not do anything
      if (Atom%fixp) return

      !
      ! Call real SEEI with current ALI configuration
      !
#ifdef DEBUGSEE
      call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                       iz,lALI,ALIp,ALIao,try_no_ALI,INPUT)
#else
      call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                       iz,lALI,ALIp,ALIao,try_no_ALI)
#endif

      !
      ! If we need to try without ALI now because there was a
      ! negative population, do it now
      !
      if (try_no_ALI) then

#ifdef DEBUGSEE
        call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                         iz,.False.,.False.,.False.,try_no_ALI,INPUT)
#else
        call SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                         iz,.False.,.False.,.False.,try_no_ALI)
#endif

      end if

      end subroutine SEEI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the statistical equilibrium equations for the atomic
      !! populations\n
      !!     Atom(Atom_class): Structure with atomic data\n
      !!    Atom0(Rhoc_class): Structure to store the density matrix
      !!                       of the previous iteration\n
      !!      JRad(double(:)): Mean intensity integrated over the
      !!                       absorption profile\n
      !!     JRadS(double(:)): Mean intensity integrated over the
      !!                       emission profile\n
      !!   JPhot(double(:,:)): Intensity integrals in the
      !!                       photoionization rates\n
      !!    LamL(double(:,:)): Lambda operator for bound-bound
      !!                       transitions\n
      !!  LamP(double(:,:,:)): Lambda operator for bound-free
      !!                       transitions\n
      !!          iz(integer): Height index\n
      !!        lALI(logical): If to apply ALI\n
      !!        ALIp(logical): If to apply ALI to bound-free\n
      !!       ALIao(logical): Switch off ALI if negative
      !!                       populations\n
      !!  try_no_ALI(logical): If we need to try again without ALI
#ifdef DEBUGSEE
      subroutine SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                             iz,lALI,ALIp,ALIao,try_no_ALI,INPUT)
#else
      subroutine SEEI_actual(Atom,Atom0,JRad,JRadS,Jphot,LamL,LamP, &
                             iz,lALI,ALIp,ALIao,try_no_ALI)
#endif

      ! I/O

#ifdef DEBUGSEE
      type(Input_class), intent(in):: INPUT
#endif
      type(Atom_class), intent(inout):: Atom
      type(Rhoc_class), intent(in):: Atom0
      logical, intent(in):: lALI,ALIp,ALIao
      logical, intent(out):: try_no_ALI
      integer, intent(in):: iz
      double precision, dimension(:), intent(in):: JRad,JRadS
      double precision, dimension(:,:), intent(in):: LamL
      double precision, dimension(:,:), intent(in):: Jphot
      double precision, dimension(:,:,:), intent(in):: LamP

      ! Local

      double precision:: rho(Atom%nlevel)
      double precision:: STcoeff(Atom%nlevel,Atom%nlevel)


      ! Build SEE equations for multi-level atom
      call SEbuildI(Atom,JRad,JRadS,Jphot,STcoeff,iz)

#ifdef DEBUGSEE
      call dump_see(Atom,STcoeff,INPUT%folder,iz,.False.)
#endif
      ! Add the contributions due to Lambda operator
      if (lALI) call ALIbuildI(Atom,Atom0,LamL,LamP,STcoeff,iz,ALIp)

#ifdef DEBUGSEE
      if (lALI) &
        call dump_see(Atom,STcoeff,INPUT%folder,iz,.True.)
#endif

      ! Initialize rho and fix populations if requested
      call initrhoI(Atom,Atom0,STcoeff,rho,iz)

      ! Solve the SEE
      call densmatrI(rho,Atom%nlevel,STcoeff)

      ! Rearrange the solution into the rhoKQ matrices
      call rhosolI(Atom,rho,iz,lALI,ALIao,try_no_ALI)

      end subroutine SEEI_actual

!#####################################################################
!#####################################################################
!#####################################################################

      !> Build the statistical equilibrium equations\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!       RadJ(double(:)): Mean intensity integrated over the
      !!                        absorption profile\n
      !!      RadJS(double(:)): Mean intensity integrated over the
      !!                        emission profile\n
      !!       JP(double(:,:)): Intensity integrals in the
      !!                        photoionization rates\n
      !!  STcoeff(double(:,:)): Statistical equilibrium equations\n
      !!           iz(integer): Height index
      subroutine SEbuildI(Atom,RadJ,RadJS,JP,STcoeff,iz)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: iz
      double precision,dimension(:), intent(in):: RadJ,RadJS
      double precision,dimension(:,:), intent(in):: JP
      double precision, dimension(:,:), intent(out):: STcoeff

      ! Local

      logical:: zpermit,zpermitJ,zPpermit,zRpermitJ,zJ,zflag
      logical:: zrelax,zupper,zlower,zrelaxJ,zupperJ,zlowerJ

      integer:: itran,ftran,fftran,iphot,iterm,itterm,iJ,iJJ,i,ii

      double precision:: rL,S,rJ,rLL,SS,rJJ
      double precision:: rSCcoeff,rACcoeff,tSCcoeff,tACcoeff
      double precision:: rEcoeff,tEcoeff,rScoeff,rAcoeff
      double precision:: tScoeff,tAcoeff,rEPcoeff,tEPcoeff
      double precision:: rSPcoeff,rAPcoeff,tSPcoeff,tAPcoeff


      ! Routine name
      urou = 'SEbuildI'

      ! Initialize SEE coefficients
      STcoeff = 0d0

      !
      ! Build SEE
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! Skip last ion for zero_ion
        if (Atom%zero_ion) then

          ! If we are in the row for the last ion, skip
          if (Atom%stage(iterm).eq.Atom%stage(Atom%nMulti)) cycle

        end if ! Zero_ion

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

            ! Skip last ion for zero_ion
            if (Atom%zero_ion) then

              ! If we are in the column for the last ion, skip
              if (Atom%stage(itterm).eq.Atom%stage(Atom%nMulti)) cycle

            end if ! No correction

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

                ! If up-down
                if (zupper) then

                  ! Get index
                  ftran = Atom%fst(itran)%irad(iJJ,iJ)

                ! If down-up
                else

                  ! Get index
                  ftran = Atom%fst(itran)%irad(iJ,iJJ)

                end if ! up->down or inverse

                ! If valid, get rolling index
                if (ftran.gt.0) fftran = Atom%ifst_ij(ftran,itran)

              ! If not valid
              else

                ! Dummy indexes
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
      ! Reset the Indentation
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

              ! Up-down b-b transfer rate
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

            ! Down-up b-b transfer rate
            STcoeff(i,ii) = tAcoeff

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! Absorption b-f transfer rate
            call tAPI(rJ,rJJ,JP(iphot,1),tAPcoeff)

            ! Down-up b-f transfer rate
            STcoeff(i,ii) = STcoeff(i,ii) + tAPcoeff

          end if ! b-f transition
        end if ! row > column
      end if ! There is a transition

      ! Diagonal element
      if (zrelax) then

        ! Diagonal level
        if (zJ) then

          ! Spontaneous relaxation rates
          STcoeff(i,ii) = -rEcoeff - rEPcoeff

          ! If stimulated emission
          if (stm) then

            ! Stimulated b-b relaxation rate
            call rSI(iterm,iJ,Atom%irad(:,iterm),Atom%fst, &
                     Atom%ifst_ij,Atom%nJ,Atom%rJval,rJ, &
                     RadJS,rScoeff)

          ! No stimulated
          else

            ! No rate
            rScoeff = 0d0

          end if ! Stimulated emission

          ! Absorption b-b relaxation rate
          call rAI(iterm,iJ,Atom%irad(:,iterm),Atom%fst, &
                   Atom%ifst_ij,Atom%nJ,Atom%nMulti,RadJ,rAcoeff)

          ! Relaxation b-b
          STcoeff(i,ii) = -rScoeff - rAcoeff + STcoeff(i,ii)

          ! Stimulated b-f relaxation rate
          call rSPI(i,Atom%iphot(:,i),JP(:,2),rSPcoeff)

          ! Absorption b-f relaxation rate
          call rAPI(i,Atom%iphot(:,i),Atom%nlevel, &
                    JP(:,1),rAPcoeff)

          ! Relaxation b-f
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

          ! Up-down collisional
          STcoeff(i,ii) = tSCcoeff + STcoeff(i,ii)

        ! row > column
        else if (zlowerJ) then

          ! Colisional excitation transfer rate
          call tAFCI(rJ,rJJ,Atom%CcoeffJ(ii,i,iz),tACcoeff)

          ! Down-up collisional
          STcoeff(i,ii) = STcoeff(i,ii) + tACcoeff

        ! diagonal element
        else if (zrelaxJ) then

          ! Colisional de-excitation relaxation rate
          call rSFCI(i,Atom%CcoeffJ(i,:,iz),rSCcoeff)

          ! Colisional excitation relaxation rate
          call rAFCI(i,Atom%nlevel,Atom%CcoeffJ(i,:,iz),rACcoeff)

          ! Relaxation collisional
          STcoeff(i,ii) = -rSCcoeff - rACcoeff + STcoeff(i,ii)

        end if ! Relation between levels
      end if ! There are collisions

      ! Check if there is any non-zero coefficients in this row
      if (zflag.and.(abs(STcoeff(i,ii)).gt..0D0)) zflag=.false.

              !
              ! Recover the indentation
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

            ! If not aborted yet, fail now
            write(umsg,*) ' # Element',iterm,real(rJ), &
                          ' isolated'
            call abortedS(umsg,urou,.True.,.True.)
            return

          end if ! Singular matrix

        end do ! iJ
      end do ! iterm

      ! Check the matrix is actually square
      if (i.ne.ii) then

        ! Issue error
        umsg = 'STcoeff is not square'
        call abortedS(umsg,urou,.True.,.True.)
        return

      end if ! Non-square matrix

      return

      end subroutine SEbuildI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add the Accelerated Lambda Iteration terms to the statistical
      !! equilibrium equations\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!     Atom0(Rhoc_class): Structure to store the density matrix
      !!                        of the previous iteration\n
      !!     LamL(double(:,:)): Lambda operator for bound-bound
      !!                        transitions\n
      !!   LamP(double(:,:,:)): Lambda operator for bound-free
      !!                        transitions\n
      !!  STcoeff(double(:,:)): Statistical equilibrium equations\n
      !!           iz(integer): Height index\n
      !!         ALIp(logical): If to apply ALI to bound-free
      subroutine ALIbuildI(Atom,Atom0,LamL,LamP,STcoeff,iz,ALIp)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Rhoc_class), intent(in):: Atom0
      logical, intent(in):: ALIp
      integer, intent(in):: iz
      double precision,dimension(:,:), intent(in):: LamL
      double precision,dimension(:,:,:), intent(in):: LamP
      double precision, dimension(:,:), intent(inout):: STcoeff

      ! Local

      logical:: zpermit,zpermitJ,zPpermit,zRpermitJ
      logical:: zrelax,zupper,zlower,zJ

      integer:: itran,ftran,fftran,iphot,iterm,itterm
      integer:: iJ,iJJ,i,ii,iR,iRR

      double precision:: rL,S,rJ,rLL,SS,rJJ,tAcoeff,tAPcoeff


      !
      ! Build ALI contribution to SEE
      !

      ! For each term (row)
      do iterm=1,Atom%nMulti

        ! Skip last ion for zero_ion
        if (Atom%zero_ion) then

          ! If we are in the row for the last ion, skip
          if (Atom%stage(iterm).eq.Atom%stage(Atom%nMulti)) cycle

        end if ! Zero_ion

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

            ! Skip last ion for zero_ion
            if (Atom%zero_ion) then

              ! If we are in the column for the last ion. skip
              if (Atom%stage(itterm).eq.Atom%stage(Atom%nMulti)) cycle

            end if ! Zero_ion

            ! Skip if not doing photoionization ALI
            if (.not.ALIp) then

              ! Skip if different ions
              if (Atom%stage(itterm).eq.Atom%stage(iterm)) cycle

            end if ! Not doing ALI in bound-free

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

                ! If up-down
                if (zupper) then

                  ! Get index
                  ftran = Atom%fst(itran)%irad(iJJ,iJ)

                ! If down-up
                else

                  ! Get index
                  ftran = Atom%fst(itran)%irad(iJ,iJJ)

                end if ! Up-down or reversed

                ! If valid transition, get rolling index
                if(ftran.gt.0) fftran = Atom%ifst_ij(ftran,itran)

              ! No valid transition
              else

                ! Dummy indexes
                ftran = 0
                fftran = 0

              end if ! Valid transitions

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
      ! Reset the Indentation
      !

      ! If there is a transition between levels
      if (zpermitJ) then

        ! If column > row
        if (zupper) then

          ! If it is a b-b transition
          if (zRpermitJ) then

            ! ALI contribution is like absorption transfer rate
            call tAI(0d0,0d0,Atom%fst(itran)%Blu(iJ,iJJ), &
                     LamL(1,fftran),tAcoeff)

            ! ALI correction
            STcoeff(i,ii) = STcoeff(i,ii) - &
                            tAcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! ALI contribution is like absorption transfer rate
            call tAPI(0d0,0d0,LamP(1,iphot,1),tAPcoeff)

            ! ALI correction
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

            ! ALI correction
            STcoeff(i,ii) = STcoeff(i,ii) - &
                            tAcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-b transition

          ! If it is a b-f transition
          if (zPpermit) then

            ! ALI contribution is like absorption transfer rate
            call tAPI(rJ,rJJ,LamP(1,iphot,1),tAPcoeff)

            ! ALI correction
            STcoeff(i,ii)= STcoeff(i,ii) - &
                           tAPcoeff*dble(Atom0%crho(iR,iz))
            STcoeff(i,i) = STcoeff(i,i) + &
                           tAPcoeff*dble(Atom0%crho(iRR,iz))

          end if ! b-f transition
        end if ! row > column
      end if ! There is a transition

              !
              ! Recover indentation
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

      !> Initialize independent vector and manage additional atomic
      !! flags\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!     Atom0(Rhoc_class): Structure to store the density matrix
      !!                        of the previous iteration\n
      !!  STcoeff(double(:,:)): Statistical equilibrium equations\n
      !!        rho(double(:)): Independent vector of the statistical
      !!                        equilibrium equations\n
      !!           iz(integer): Height index
      subroutine initrhoI(Atom,Atom0,STcoeff,rho,iz)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Rhoc_class), intent(in):: Atom0
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

      end if ! Zero ion

      return

      end subroutine initrhoI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the statistical equilibrium equations\n
      !!        rho(double(:)): Independent vector and solution of the
      !!                        statistical equilibrium equations\n
      !!         ndim(integer): Dimensionality of the statistical
      !!                        equilibrium system\n
      !!  STcoeff(double(:,:)): Statistical equilibrium equations
      subroutine densmatrI(rho,ndim,STcoeff)

      ! I/O

      integer, intent(in):: ndim
      double precision, dimension(:), intent(out):: rho
      double precision, dimension(:,:), intent(in):: STcoeff

      ! Local

      integer:: i
      integer, dimension(ndim):: indx


      ! Solve SEE
      call DGESV(ndim,1,STcoeff,ndim,indx,rho,ndim,i)

      end subroutine densmatrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Transform the solution of the statistical equilibrium
      !! equations into atomic populations\n
      !!     Atom(Atom_class): Structure with atomic data\n
      !!       rho(double(:)): Solution of the statistical equilibrium
      !!                       equations\n
      !!          iz(integer): Height index\n
      !!        lALI(logical): If to apply ALI\n
      !!       ALIao(logical): Switch off ALI if negative
      !!                       populations\n
      !!  try_no_ALI(logical): If we need to try again without ALI
      subroutine rhosolI(Atom,rho,iz,lALI,ALIao,try_no_ALI)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      logical, intent(in):: lALI,ALIao
      logical, intent(out):: try_no_ALI
      integer, intent(in):: iz
      double precision, dimension(:), intent(in):: rho

      ! Local

      integer:: i,it,iJ,iR

      double precision:: rJ


      ! Routine name
      urou = 'rhosolI'

      ! Initialize check flag
      try_no_ALI = .False.

      ! Initialize null flag
      Atom%rhonull(:,iz) = .True.

      ! Check negative populations
      if (minval(rho).lt.0d0) then

        ! If allowing
        if (nphysR) then

          ! Issue warning
          write(umsg,'(A)') 'Negative population in SEEI '// &
                            'solution, but allowed'
          call abortedS(umsg,urou,.False.,.True.)

        ! Not allowing
        else

          ! If doing ALI
          if (lALI) then

            ! If can swith if off
            if (ALIao) then

              ! Issue warning only if synthesis
              if (run_mode.ne.-1) then
                write(umsg,'(A)') 'Negative population in SEEI '// &
                                  'solution, will try without ALI'
                call abortedS(umsg,urou,.False.,.True.)
              end if

              ! Flag and go back
              try_no_ALI = .True.
              return

            ! If cannot swith it off
            else

              ! Find
              do i=1,Atom%nlevel

                ! Check negativity
                if(rho(i).lt.0d0)then

                  ! Issue error
                  write(umsg,'(A,i4,",",i4,A,1x,es11.4)') &
                    'Negative population in SEE solution'// &
                    new_line('A')// &
                    '(iz,il)=(',iz,i,')'// &
                    new_line('A')//'rho00: ',rho(i)
                  call abortedS(umsg,urou,.not.nphysR,.True.)

                end if ! Negative population at this heright

              end do ! Levels

              ! Go back, we are aborting
              return

            end if ! Can swith ALI off

          ! Not doing ALI
          else

            ! Find
            do i=1,Atom%nlevel

              ! Check negativity
              if(rho(i).lt.0d0)then

                ! Issue error
                write(umsg,'(A,i4,",",i4,A,1x,es11.4)') &
                  'Negative population in SEE solution'// &
                  new_line('A')// &
                  '(iz,il)=(',iz,i,')'// &
                  new_line('A')//'rho00: ',rho(i)
                call abortedS(umsg,urou,.not.nphysR,.True.)

              end if ! Negative population at this heright

            end do ! Levels

            ! Go back, we are aborting
            return

          end if ! Doing ALI
        end if ! Allowed negative
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
          if(ieee_is_nan(rho(i)))then

            ! Issue error
            write(umsg,'(A,i4,",",i4,A,1x,es11.4)') &
              'NaN in SEE solution'//new_line('A')// &
              '(iz,il)=(',iz,i,')'// &
              new_line('A')//'rho00: ',rho(i)
            call abortedS(umsg,urou,.True.,.True.)

          end if ! NaN

          ! Store the corresponding element
          Atom%crho(iR,iz) = dcmplx(rho(i),0d0)

          ! Store normalized population
          Atom%popu(i,iz) = sqrt(2d0*rJ + 1d0)*rho(i)

          ! And cancel the flag
          Atom%rhonull(iR,iz) = .False.

        end do ! levels
      end do ! terms

      end subroutine rhosolI

#ifdef DEBUGSEE
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump SEE matrix into a file.\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!  folder(character(500)): Path to the output folder\n
      !!             iz(integer): Height index\n
      !!           lALI(logical): If to apply ALI
      subroutine dump_see(Atom,STcoeff,folder,iz,lALI)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: folder
      logical, intent(in):: lALI
      integer, intent(in):: iz
      double precision, dimension(:,:):: STcoeff

      ! Local

      character(len=500):: filename,formating

      logical:: exists

      integer:: ilevel


      ! Get file name for 1D
      if (run_mode.eq.0) then

        ! Get file name
        write(filename,'(A,I0.5)') &
            trim(folder)//'/debug_see_',gpid

      ! Get file name for rest
      else

        ! Get file name
        write(filename,'(A,I0.7,A,I0.5)') &
            trim(folder)//'/debug_see_',icoords(3),'_',gpid

      end if ! Run mode

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)

      ! If no exist
      if(.not.exists) then

        ! Open new
        open(800,file=trim(filename))

      ! If exists
      else

        ! Open old
        open(800,file=trim(filename),position='append')

      endif ! File exists

      ! Write header
      write(800,*) ''
      write(800,*) ''
      if (lALI) then
        write(800,'("Height node (ALI)",2(1x,i4))') iz,iz-Rz0+1
      else
        write(800,'("Height node",2(1x,i4))') iz,iz-Rz0+1
      end if

      ! Define formating
      if (Atom%nlevel.lt.10) then
        write(formating,'("(i1,1x,",i1,"(1x,es13.5))")') Atom%nlevel
      else if (Atom%nlevel.lt.100) then
        write(formating,'("(i2,1x,",i2,"(1x,es13.5))")') Atom%nlevel
      else if (Atom%nlevel.lt.1000) then
        write(formating,'("(i3,1x,",i3,"(1x,es13.5))")') Atom%nlevel
      else
        write(formating,'("(i4,1x,",i4,"(1x,es13.5))")') Atom%nlevel
      end if

      ! For each row
      do ilevel=1,Atom%nlevel
        write(800,trim(formating)) ilevel,STcoeff(ilevel,:)
      end do

      ! Close
      close(800)

      end subroutine dump_see

#endif
!#####################################################################
!#####################################################################
!#####################################################################

      end module seei_mod
