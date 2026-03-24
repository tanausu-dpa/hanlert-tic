      !> Collisions initialization
      module initcols_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     25/09/2019
!  Last version:
!     24/03/2026 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     24/03/2026:    V4.0.3 - Made the necessary modifications to
!                             account for the changes in colinter,
!                             that can fallback to Cubic Herminte
!                             interpolation before linear (TdPA)
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
!    Initcols
!      Computes collisional rates
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use funnj_mod
      use inter_mod
      use parameters_mod , only : VTINY , PI , c , c2
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the collisional rates for in a given model
      !! atmosphere and a given atomic model\n
      !!        Atom(Atom_class): Structure with atomic data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!  folder(character(500)): Path to the output folder\n
      !!      Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                          and J-symbols\n
      !!        keeplog(logical): If keeping an ASCII log file for the
      !!                          collisional rates\n
      !!         active(logical): If the current atom is active
      subroutine Initcols(Atom,Atmo,folder,Flgsg,keeplog,active)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Fctsg_class), intent(in):: Flgsg
      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: folder
      logical, intent(in):: active, keeplog

      ! Local

      character(len=2):: c2dump

      logical:: lin,lcub,free,aout

      integer:: ios,iz,it,itt,i,ii,iJ,iJJ,icol,iJ1,ilevel,iterm,jj
      integer:: ilevel1,iterm1,up,low,stagl,stagu,K,minK,maxK

      double precision:: rJ,rJ1,rL,rL1,S,S1,El,Eu,gl,gu,Dfreq2,W6
      double precision:: d1,d2,d3
      double precision, dimension(nz):: nu, nl
      double precision, dimension(:), allocatable:: p_pop, CulI, Culin

      ! Pointers

      type(tmp_col_box_class), pointer:: p_col, p_col_p
      type(Tbox_class), pointer:: p_T, p_T_p

      ! Initialize pointers
      nullify(p_col,p_col_p,p_T,p_T_p)

      ! If the collisional data can be fred when finished
      free = run_mode.eq.0

      ! If we can output ASCII file with log of collisions
      aout = run_mode.eq.0.and.keeplog


      !
      ! Initialize atomic level inverse lifetime
      !
      allocate(Atom%damp(Atom%nMulti,nz))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%damp)
      Atom%damp = 0d0

      ! For each upper term
      do itt=2,Atom%nMulti

        ! For each lower term
        do it=1,itt-1

          ! If a registered transition exists
          if (Atom%irad(itt,it).gt.0) then

            ! Here the damping parameter is set to the theoretical
            ! minimum based on the value of the Einstein A-coefficient
            Atom%damp(itt,:) = Atom%damp(itt,:) + &
                               1d-8*Atom%Ecoeff(itt,it)/c/(4d0*PI)

          end if

        end do ! Lower term
      end do ! Upper term


      !
      ! Calculate elastic collisional rates
      !

      ! If active or there are rates
      if (active.or.Atom%ngk.ge.1) then

        ! Allocation
        ! g^(K) in the SEE
        allocate(Atom%gk(Atom%nMulti,Atom%nJmax,Atom%nJmax, &
                         0:Atom%nKmax,nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%gk)
        Atom%gk = 0d0

      end if

      ! If there are inputs
      if (Atom%ngk.ge.1) then

        ! For each entry
        do ii=1,Atom%ngk

          ! Get level, term, and sublevel indexes
          ilevel = Atom%elas(ii)%ilevel
          iterm = Atom%term(ilevel)
          iJ = Atom%sublevel(ilevel)

          ! For each line in this entry
          do jj=1,Atom%elas(ii)%nentry

            ! Read the multipole, the type of input and the
            ! dimensionality
            K = Atom%elas(ii)%datum(jj)%K

            ! If it is a fit input
            if(Atom%elas(ii)%datum(jj)%typo.eq.0)then

              ! Read the coefficients
              d1 = Atom%elas(ii)%datum(jj)%a
              d2 = Atom%elas(ii)%datum(jj)%b
              d3 = Atom%elas(ii)%datum(jj)%c

              ! For each height
              do iz=1,nZ

                ! Calculate the rate following the fit
                Atom%gk(iterm,iJ,iJ,K,iz) =  &
                                           d1*1d-9* &
                                           ((Atmo%T(iz)*2d-4)**d2)* &
                                           (d3**(Atmo%T(iz)*2d-4))* &
                                           sum(Atmo%nh(iz,1:5))*1d-8 &
                                          + Atom%gk(iterm,iJ,iJ,K,iz)
              end do ! Heights

            end if ! Fit input

            ! If the multipole was 0
            if (K.eq.0) then

              ! For each height
              do iz=1,nZ

                ! Add this to the damping, taking a weighted average
                Atom%damp(iterm,iz) = 1d-8* &
                                  Atom%gk(iterm,iJ,iJ,0,iz)* &
                                  (2d0*Atom%rJval(iJ,iterm) + 1d0)/ &
                                  Atom%deg(iterm)/c/(4d0*PI) + &
                                  Atom%damp(iterm,iz)
              end do ! Heights

            end if ! Multipole is 0

          end do ! Input sub entry

          ! If can free memory
          if (free) then

            ! Deallocate
            MRAMc = MRAMc - 1d-6*sizeof(Atom%elas(ii)%datum)
            deallocate(Atom%elas(ii)%datum)
            MRAMc = MRAMc - 1d-6*sizeof(Atom%elas(ii))

          end if ! Can free memory

        end do ! Input entry

        ! If can free memory, deallocate
        if (free) deallocate(Atom%elas)

        ! Only if active atoms
        if (active) then

          !
          ! Heuristic definition of non-diagonal terms
          !

          ! For each term
          do iterm=1,Atom%nMulti

            ! For each sublevel
            do iJ=1,Atom%nJ(iterm)

              ! For each other sublevel
              do iJ1=1,Atom%nJ(iterm)

                ! Skip diagonals
                if(iJ.eq.iJ1)cycle

                ! Calculate limits for multipoles
                minK = nint(abs(Atom%rJval(iJ,iterm) - &
                                Atom%rJval(iJ1,iterm)))
                maxK = nint(Atom%rJval(iJ,iterm) + &
                            Atom%rJval(iJ1,iterm))

                ! For every height
                do iz=1,nZ

                  ! For each available multipole
                  do K=minK,maxK

                    ! If both levels have 0 rate
                    if(Atom%gk(iterm,iJ,iJ,K,iz).lt.1d-100.and. &
                       Atom%gk(iterm,iJ1,iJ1,K,iz).lt.1d-100)then

                      ! Cross-rate is zero as well
                      Atom%gk(iterm,iJ,iJ1,K,iz) = 0d0

                    ! If one of the levels have 0 rate
                    else if (Atom%gk(iterm,iJ,iJ,K,iz).lt. &
                             1d-100.or. &
                             Atom%gk(iterm,iJ1,iJ1,K,iz).lt. &
                             1d-100) then

                      ! If the iJ has zero rate
                      if(Atom%gk(iterm,iJ,iJ,K,iz).lt.1d-100)then

                        ! The rate is the one for the iJ1
                        Atom%gk(iterm,iJ,iJ1,K,iz) = &
                                           Atom%gk(iterm,iJ1,iJ1,K,iz)

                      ! If the iJ1 has zero rate
                      else

                        ! The rate is the one for the iJ
                        Atom%gk(iterm,iJ,iJ1,K,iz) = &
                                           Atom%gk(iterm,iJ,iJ,K,iz)

                      end if ! Which rate is zero

                    ! If both levels have non-zero rate
                    else

                      ! Take the average
                      Atom%gk(iterm,iJ,iJ1,K,iz) = &
                              (Atom%gk(iterm,iJ,iJ,K,iz)* &
                               (2d0*Atom%rJval(iJ,iterm) + 1d0) + &
                               Atom%gk(iterm,iJ1,iJ1,K,iz)* &
                               (2d0*Atom%rJval(iJ1,iterm) +.1d0))/ &
                              (2d0*Atom%rJval(iJ,iterm) + 1d0 + &
                               2d0*Atom%rJval(iJ1,iterm) + 1d0)

                    end if ! Which rates are zero

                  end do ! K
                end do ! heights
              end do ! iJ1
            end do ! iJ
          end do ! Term

        ! If passive atom
        else

          ! Deallocate depolarizing rates
          MRAMc = MRAMc - 1d-6*sizeof(Atom%gk)
          deallocate(Atom%gk)

        end if ! Active atom
      end if ! There are elastic collisions


      !
      ! Calculate the down-up ionizing collisions and
      ! up-down exciting collisions
      !


      ! If active atom or there are collisions
      if (active.or.Atom%ncol.ge.1) then

        ! Allocations

        ! Indexing of collisions between terms
        allocate(Atom%icol(Atom%nMulti,Atom%nMulti))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%icol)
        Atom%icol = 0

        ! Collisional rates for collisions between terms
        allocate(Atom%Ccoeff(Atom%nMulti,Atom%nMulti,nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%Ccoeff)
        Atom%Ccoeff = 0d0

        ! Collisional rates for collisions between levels
        allocate(Atom%CcoeffJ(Atom%nlevel,Atom%nlevel,nZ))
        MRAMc = MRAMc + 1d-6*sizeof(Atom%CcoeffJ)
        Atom%CcoeffJ = 0d0

      end if ! Active atom or there are collisions

      !
      ! Inelastic collisions
      !

      ! If there are inelastic collisions
      if (Atom%ncol.ge.1) then

        ! Allocate auxiliar
        allocate(p_pop(nz))
        allocate(CulI(nz))
        allocate(Culin(nz))
        ! Nullify pointer in Ccoeff_special
        nullify(Atom%Ccoeff_special)

        ! Initialize pointer
        p_T => Atom%Tbox

        ! For every collision
        do icol=1,Atom%ncol

          ! Get temperature box
          if (p_T%ind.ne.Atom%inelas(icol)%ind) then

            ! If behind, reset pointer
            if (p_T%ind.gt.Atom%inelas(icol)%ind) then
              p_T => Atom%Tbox
            end if

            ! Go to next until found
            do while (p_T%ind.ne.Atom%inelas(icol)%ind)
              p_T => p_T%next
            end do

          end if ! Incorrect temperature box

          ! Get level or term indexes
          up = Atom%inelas(icol)%up
          low = Atom%inelas(icol)%low

          ! If between terms
          if (p_T%col_type.eq.0) then

            ! Degeneracy
            gu = Atom%deg(up)
            gl = Atom%deg(low)

            ! Index this transition
            Atom%icol(up,low) = icol
            Atom%icol(low,up) = icol

            ! Stages
            stagl = Atom%stage(low)
            stagu = Atom%stage(up)

          ! If between levels
          else if (p_T%col_type.eq.1) then

            ! Get the S, L, and J of levels
            S = Atom%Sval(Atom%term(low))
            S1 = Atom%Sval(Atom%term(up))
            rL = Atom%rLval(Atom%term(low))
            rL1 = Atom%rLval(Atom%term(up))
            rJ = Atom%rJval(Atom%sublevel(low),Atom%term(low))
            rJ1 = Atom%rJval(Atom%sublevel(up),Atom%term(up))

            ! Get energy
            El = Atom%FSfreq(Atom%sublevel(low),Atom%term(low))
            Eu = Atom%FSfreq(Atom%sublevel(up),Atom%term(up))

            ! Degeneracy
            gu = (2d0*rJ1 + 1d0)
            gl = (2d0*rJ + 1d0)

            ! Flag this collision as forbidden or not given the
            ! input
            Atom%fcflag(up,low) = Atom%inelas(icol)%forbid
            Atom%fcflag(low,up) = Atom%inelas(icol)%forbid

            ! Stages
            stagl = Atom%stage(Atom%term(low))
            stagu = Atom%stage(Atom%term(up))

            ! Check if they are "forbidden" because of being
            ! ionizing
            if (stagu.ne.stagl) then

              ! Ionizing is 2
              Atom%fcflag(up,low) = 2
              Atom%fcflag(low,up) = 2

            ! Check if they are forbidden through dipole selection
            ! rules
            else

              ! Check if same term
              if (Atom%term(low).eq.Atom%term(up)) then

                ! Dipole forbidden is 1
                Atom%fcflag(up,low) = 1
                Atom%fcflag(low,up) = 1

              ! Different terms
              else

                ! Check if intercombination
                if (nint(abs(S1-S)).gt.0) then

                  ! Dipole forbidden (for pol.) is 1
                  Atom%fcflag(up,low) = 1
                  Atom%fcflag(low,up) = 1

                ! If spin allowed
                else

                  ! Check orbital selection rule
                  if (nint(abs(rL - rL1)).gt.1.or. &
                      nint(rL + rL1).eq.0) then

                    ! Dipole forbidden is 1
                    Atom%fcflag(up,low) = 1
                    Atom%fcflag(low,up) = 1

                  ! If allowed by orbital angular momentum
                  else

                    ! Check total angular momentum selection rule
                    if (nint(abs(rJ1 - rJ)).gt.1.or. &
                        nint(rJ1 + rJ).eq.0) then

                      ! Dipole forbidden is 1
                      Atom%fcflag(up,low) = 1
                      Atom%fcflag(low,up) = 1

                    ! If allowed by total angular momentum
                    else

                      ! If not radiatively connected
                      if (Atom%irad(Atom%term(up), &
                                    Atom%term(low)).lt.1) then

                        ! Assume same parity, forbidden
                        Atom%fcflag(up,low) = 1
                        Atom%fcflag(low,up) = 1

                      end if ! Radiatively connected
                    end if ! Angular momentum rules
                  end if ! Orbital selection rules
                end if ! Spin selection rules
              end if ! Same term
            end if ! Ionization check
          end if ! If between levels


          !
          ! b-b symmetric
          !

          ! If it is a symmetric b-b excitation
          if (Atom%inelas(icol)%col_type.eq.0.or. &
              Atom%inelas(icol)%col_type.eq.2.or. &
              Atom%inelas(icol)%col_type.eq.3) then

            ! Collision with electrons
            if (Atom%inelas(icol)%col_type.eq.0) then
              c2dump = 'BE'
              p_pop = Atmo%ne
            ! Collision with neutral hydrogen
            else if (Atom%inelas(icol)%col_type.eq.2) then
              c2dump = 'BP'
              p_pop = Atmo%nh(:,6)
            ! Collision with protons
            else if (Atom%inelas(icol)%col_type.eq.3) then
              c2dump = 'BH'
              p_pop = Atmo%nh(:,1)
            end if

            ! If neutral or ion flag
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              ! Get input and multiply by upper level/term degeneracy
              Culin = Atom%inelas(icol)%Cul*gu

              ! If neutral scale mode
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then

                ! Divide by degeneracy lower level/term and
                ! temperature factor
                Culin = Culin/gl
                Culin = Culin/sqrt(p_T%temp)

              ! If ion scale mode
              else

                ! Multiply by temperature factor
                Culin = Culin*sqrt(p_T%temp)

              end if ! Neutral or ion scale mode

            ! No scale mode
            else

              ! Keep just the rate
              Culin = Atom%inelas(icol)%Cul

            end if ! Type of scaling

            ! Interpolate in the tabulation
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,CulI,nZ,p_T%flin, &
                          lcub,lin)

            ! If the interpolation ended up being linear or CH
            ! and it is master
            if ((lin.or.lcub).and.pid.eq.0) then

              ! Linear
              if (lin) then

                ! If term to term transition
                if (p_T%col_type.eq.0) then

                  ! Write message
                  write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                       ' # Inelastic collisional rate '// &
                       c2dump//' ',up,low,'between terms in atom', &
                       Atom%Element,'was negative with Spline '// &
                       'interpolation, did linear.'

                ! If level to level transition
                else

                  ! Write message
                  write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                      ' # Inelastic collisional rate '// &
                      c2dump//' ',up,low,'between levels in atom', &
                      Atom%Element,'was negative with Spline '// &
                      'interpolation, did linear.'

                end if ! Term-term or level-level

              ! CH
              else

                ! If term to term transition
                if (p_T%col_type.eq.0) then

                  ! Write message
                  write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                       ' # Inelastic collisional rate '// &
                       c2dump//' ',up,low,'between terms in atom', &
                       Atom%Element,'was negative with Spline '// &
                       'interpolation, did Cubic Hermite.'

                ! If level to level transition
                else

                  ! Write message
                  write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                      ' # Inelastic collisional rate '// &
                      c2dump//' ',up,low,'between levels in atom', &
                      Atom%Element,'was negative with Spline '// &
                      'interpolation, did Cubir Hermite.'

                end if ! Term-term or level-level
              end if ! Linear or CH

              ! Verbose warning
              call verbose

            end if ! If had to do CH or linear and is master

            ! If neutral or ion scale mode
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              ! Divide by upper level degeneracy
              CulI = CulI/gu

              ! If neutral scale mode
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then

                ! Revert products
                CulI = CulI*gl
                CulI = CulI*sqrt(Atmo%T)

              ! If ion scale mode
              else

                ! Revert products
                CulI = CulI/sqrt(Atmo%T)

              end if ! Scale mode
            end if ! If any scale mode

            ! Convertion to 10^8 s-1
            CulI = CulI*1d-8

            ! If collision between terms
            if (p_T%col_type.eq.0) then

              ! For each height
              do iz=1,nZ

                ! Compute actual rate
                Atom%Ccoeff(up,low,iz) = Atom%Ccoeff(up,low,iz) + &
                                         CulI(iz)*p_pop(iz)

                ! Add to the damping
                Atom%damp(up,iz) = Atom%damp(up,iz) + &
                                   1d-8*CulI(iz)*p_pop(iz)/c/4d0/PI

              end do ! Heights

            ! If collision between levels
            else if (p_T%col_type.eq.1) then

              ! For each height
              do iz=1,nZ

                ! Compute actual rate
                Atom%CcoeffJ(up,low,iz) = CulI(iz)*p_pop(iz) + &
                                          Atom%CcoeffJ(up,low,iz)

              end do ! Heights

            end if ! Term-term or level-level


          !
          ! b-f symmetric
          !
          else if (Atom%inelas(icol)%col_type.eq.1.or. &
                   Atom%inelas(icol)%col_type.eq.4) then

            ! Collision with electrons
            if (Atom%inelas(icol)%col_type.eq.1) then
              c2dump = 'FE'
              p_pop = Atmo%ne
            ! Collision with neutral hydrogen
            else if (Atom%inelas(icol)%col_type.eq.4) then
              c2dump = 'FH'
              p_pop = Atmo%nh(:,1)
            end if

            ! If any scale mode
            if (p_T%nion.ge.0) then

              ! Calculate difference of energies
              Dfreq2 = abs(Eu-El)*c2*1d4

              ! And scale rate
              Culin = Atom%inelas(icol)%Cul* &
                      exp(Dfreq2/p_T%temp)/ &
                      sqrt(p_T%temp)

            ! No scale
            else

              ! Keep just the rate
              Culin = Atom%inelas(icol)%Cul

            end if ! Scale mode

            ! Interpolate in the tabulation
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,CulI,nZ,p_T%flin, &
                          lcub,lin)

            ! If the interpolation ended up being linear or CH
            ! and it is master
            if ((lin.or.lcub).and.pid.eq.0) then

              ! Linear
              if (lin) then

                ! Write message
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '//c2dump//' ', &
                  low,up,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did linear.'

              ! CH
              else

                ! Write message
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '//c2dump//' ', &
                  low,up,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did Cubic Hermite.'

              end if ! Linear or CH

              ! Verbose
              call verbose

            end if ! If had to do CH or linear and is master

            ! If any scale mode
            if (p_T%nion.ge.0) then

              ! Revert products
              CulI = CulI*sqrt(Atmo%T)/exp(Dfreq2/Atmo%T)

            end if ! Any scale mode

            ! Convertion to 10^8 s-1
            CulI = CulI*1d-8

            ! For every height
            do iz=1,nZ

              ! Compute actual rate
              Atom%CcoeffJ(low,up,iz) = CulI(iz)*p_pop(iz) + &
                                        Atom%CcoeffJ(low,up,iz)

            end do ! Heights

          !
          ! Charge transfer rate
          !
          else if (Atom%inelas(icol)%col_type.eq.5.or. &
                   Atom%inelas(icol)%col_type.eq.6) then

            ! If not the first non-symmetric collision
            if (associated(Atom%Ccoeff_special)) then

              ! Create the next box and point p_col to it
              p_col => Atom%Ccoeff_special
              do while (associated(p_col%next))
                p_col => p_col%next
              end do
              allocate(p_col%next)
              p_col => p_col%next
              nullify(p_col%next)

            ! First non-symmetric collision
            else

              ! Allocate Ccoeff_special and point p_col to it
              allocate(Atom%Ccoeff_special)
              p_col => Atom%Ccoeff_special
              nullify(p_col%next)

            end if ! First non-symmetric collision

            ! If up-low collision
            if (Atom%inelas(icol)%col_type.eq.5) then

              ! Configure
              p_col%ifrom = up
              p_col%ito = low
              c2dump = 'c0'
              p_pop = Atmo%nh(:,1)

            ! If low-up collision
            else

              ! Configure
              p_col%ifrom = low
              p_col%ito = up
              c2dump = 'c+'
              p_pop = Atmo%nh(:,6)

            end if ! up->low or low->up

            ! Allocate coefficient
            allocate(p_col%C(nZ))
            MRAMc = MRAMc + 1d-6*(12 + sizeof(p_col%C))

            ! If neutral or ion scale mode
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              ! Get input and multiply by upper level/term degeneracy
              Culin = Atom%inelas(icol)%Cul*gu

              ! If neutral scale mode
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then

                ! Divide by degeneracy lower level/term and
                ! temperature factor
                Culin = Culin/gl
                Culin = Culin/sqrt(p_T%temp)

              ! If ion scale mode
              else

                ! Multiply by temperature factor
                Culin = Culin*sqrt(p_T%temp)

              end if ! Neutral or ion scale mode

            ! No scale mode
            else

                ! Multiply by temperature factor
              Culin = Atom%inelas(icol)%Cul

            end if ! Type of scaling

            ! Interpolate in the tabulation
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,p_Col%C,nZ,p_T%flin, &
                          lcub,lin)

            ! If the interpolation ended up being linear or CH
            ! and it is master
            if ((lin.or.lcub).and.pid.eq.0) then

              ! Linear
              if (lin) then

                ! Write message
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '// &
                  c2dump//' ',up,low,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did linear.'

              ! CH
              else

                ! Write message
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '// &
                  c2dump//' ',up,low,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did Cubic Hermite.'

              end if ! Linear or CH

              ! Verbose
              call verbose

            end if ! If had to do CH or linear and is master

            ! If neutral or ion scale mode
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              ! Divide by upper level degeneracy
              p_Col%C = p_Col%C/gu

              ! If neutral scale mode
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then

                ! Revert products
                p_Col%C = p_Col%C*gl
                p_Col%C = p_Col%C*sqrt(Atmo%T)

              ! If ion scale mode
              else

                ! Revert products
                p_Col%C = p_Col%C/sqrt(Atmo%T)

              end if ! Scale mode
            end if ! If any scale mode

            ! Convertion to 10^8 s-1
            p_Col%C = p_Col%C*1d-8

            ! Get rate
            p_Col%C = p_Col%C*p_pop

            ! Flag 2 (stage change)
            p_Col%flag = 2

            ! Nullify pointer
            nullify(p_col)

          end if ! Type of collision

          ! If memory can be fred
          if (free) then

            ! Deallocate table data
            MRAMc = MRAMc - 1d-6*sizeof(Atom%inelas(icol)%Cul)
            deallocate(Atom%inelas(icol)%Cul)
            MRAMc = MRAMc - 1d-6*sizeof(Atom%inelas(icol))

          end if ! Can free memory

        end do ! Collisional rates

        ! If can free memory, deallocate inelastic database
        if (free) deallocate(Atom%inelas)

        ! Deallocate Tbox if possible
        if (free) then

          ! For every collision
          do while (associated(Atom%Tbox))

            ! Initialize
            p_T => Atom%Tbox

            ! Go to the last one
            do while (associated(p_T%next))
              p_T_p => p_T
              p_T => p_T%next
            end do ! Forward navigation

            ! Deallocate temp
            MRAMc = MRAMc - 1d-6*sizeof(p_T%temp)
            deallocate(p_T%temp)

            ! Remove the last one
            if (associated(p_T_p)) then
              nullify(p_T)
              deallocate(p_T_p%next)
              nullify(p_T_p%next)
              nullify(p_T_p)
            ! We are in the last one
            else
              nullify(p_T)
              deallocate(Atom%Tbox)
              nullify(Atom%Tbox)
            end if

            ! Memory count
            MRAMc = MRAMc - 20d-6

          end do ! Every collision

        ! Not freeing data
        else

          ! Just nullify pointer
          nullify(p_T)

        end if ! Have to free data

        ! Deallocate auxiliars
        deallocate(p_pop)
        deallocate(CulI)
        deallocate(Culin)


        !
        ! Compute term b-b collisional transitions from the level
        ! wise ones
        !

        ! If multilevel atom, we just need to copy
        if (Atom%ML) then

          ! Get current maximum collision index
          icol = maxval(Atom%icol)

          ! For each lower level
          do ilevel=1,Atom%nMulti-1

            ! For each upper level
            do ilevel1=ilevel+1,Atom%nMulti

              ! If ionizing, skip
              if (Atom%stage(ilevel).ne.Atom%stage(ilevel1)) cycle

              ! If there is already a term-term collision, skip
              if (Atom%icol(ilevel,ilevel1).gt.0) cycle

              ! If this collision is forbidden, don't include it
              if (Atom%fcflag(ilevel1,ilevel).gt.0) cycle

              ! Save the rate in the term-term array
              Atom%Ccoeff(ilevel1,ilevel,:) = &
                                        Atom%CcoeffJ(ilevel1,ilevel,:)

              ! Add to the damping parameter
              Atom%damp(ilevel1,1:nZ) = 1d-8* &
                                 Atom%Ccoeff(ilevel1,ilevel,1:nZ)/c/ &
                                 (4d0*PI) + Atom%damp(ilevel1,1:nZ)

              ! Background atoms are done
              if (.not.active) cycle

              ! If the transition exists
              if (maxval(Atom%Ccoeff(ilevel1,ilevel,:)).gt.0d0) then

                ! Advance index
                icol = icol + 1

                ! And index this collisional transition
                Atom%icol(ilevel1,ilevel) = icol
                Atom%icol(ilevel,ilevel1) = icol

              end if ! There is a transition

            end do ! Upper level
          end do ! Lower level

        ! If multi-term, we need to do the averages
        else

          ! Get current maximum transition index
          icol = maxval(Atom%icol)

          ! Lower term
          do iterm=1,Atom%nMulti-1

            ! Upper term
            do iterm1=iterm+1,Atom%nMulti

              ! If ionizing, skip
              if (Atom%stage(iterm).ne.Atom%stage(iterm1)) cycle

              ! If there is already a term-term collision, skip
              if (Atom%icol(iterm,iterm1).gt.0) cycle

              ! Upper level
              do iJ1=1,Atom%nJ(iterm1)

                ! Lower level
                do iJ=1,Atom%nJ(iterm)

                  ! Get continuous level index
                  ilevel = Atom%irho(iterm)%irho_ij(iJ)
                  ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)

                  ! If this collision is forbidden, don't include it
                  ! on in the average
                  if (Atom%fcflag(ilevel1,ilevel).gt.0) cycle

                  ! Add contribution to the average
                  Atom%Ccoeff(iterm1,iterm,:) = &
                                    Atom%Ccoeff(iterm1,iterm,:) + &
                                  (2d0*Atom%rJval(iJ1,iterm1)+1d0)* &
                                  Atom%CcoeffJ(ilevel1,ilevel,:)/ &
                                  Atom%deg(iterm1)

                end do ! Lower level
              end do ! Upper level

              ! Add to the damping parameter
              Atom%damp(iterm1,1:nZ) = 1d-8* &
                                 Atom%Ccoeff(iterm1,iterm,1:nZ)/c/ &
                                (4d0*PI) + Atom%damp(iterm1,1:nZ)

              ! Background atoms are done here
              if (.not.active) cycle

              ! If there is a transition
              if (maxval(Atom%Ccoeff(iterm1,iterm,:)).gt.0d0) then

                ! Advance index
                icol = icol + 1

                ! Save transition index
                Atom%icol(iterm1,iterm) = icol
                Atom%icol(iterm,iterm1) = icol

              end if ! There is a transition

            end do ! Upper term
          end do ! Lower term

        end if ! Multi-level or multi-term

        ! If background atom
        if (.not.active) then

          ! Remove collisional data
          if (allocated(Atom%Ccoeff)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom%Ccoeff)
            deallocate(Atom%Ccoeff)
          end if
          if (allocated(Atom%CcoeffJ)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom%CcoeffJ)
            deallocate(Atom%CcoeffJ)
          end if
          if (allocated(Atom%icol)) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom%icol)
            deallocate(Atom%icol)
          end if

          ! If we can free memory, do it
          if (allocated(Atom%fcflag).and.free) then
            MRAMc = MRAMc - 1d-6*sizeof(Atom%fcflag)
            deallocate(Atom%fcflag)
          end if

          ! While there is data in non-symmetric rates
          do while (associated(Atom%Ccoeff_special))

            ! Initialize
            p_col => Atom%Ccoeff_special
            nullify(p_col_p)

            ! Go to the last one
            do while (associated(p_col%next))
              p_col_p => p_col
              p_col => p_col%next
            end do ! Forward navigation

            ! Deallocate content
            MRAMc = MRAMc - 1d-6*sizeof(p_col%C)
            deallocate(p_col%C)

            ! Remove the last one
            if (associated(p_col_p)) then
              nullify(p_col_p%next)
              nullify(p_col_p)
            ! We are in the last one
            else
              deallocate(Atom%Ccoeff_special)
              nullify(Atom%Ccoeff_special)
              nullify(p_col)
            end if

            ! Memory count
            MRAMc = MRAMc - 12d-6

          end do ! There is data

        ! If active atom
        else

          ! Output file if active, requested, and master
          if (pid.eq.0.and.aout) then

            ! Open file
            open (200,file=trim(folder)//'/'// &
                           trim(Atom%file_label)//'.COL', &
                    status='unknown', iostat=ios, action='write')

            ! Header
            write(200,'(A)') 'Information about the collisional rates'

            ! For each pair of terms
            do iterm=1,Atom%nMulti
              do iterm1=1,Atom%nMulti

                ! Only report up->down
                if (iterm.ge.iterm1) cycle

                ! Do not consider different stages here
                if (Atom%stage(iterm).ne.Atom%stage(iterm1)) cycle

                ! Check if there is not term-term collisional rate
                if (Atom%icol(iterm1,iterm).lt.1) then

                  ! Report
                  write(200,'(A)') ' '
                  write(200,'("No rate from term",1x,i4,'// &
                            '1x,"to term",1x,i4)') iterm1,iterm

                ! There is a term-term rate
                else

                  ! Report
                  write(200,'(A)') ' '
                  write(200,'("Yes rate from term",1x,i4,'// &
                            '1x,"to term",1x,i4)') iterm1,iterm

                end if ! If there is a rate

                ! For each pair of sublevels
                do iJ1=1,Atom%nJ(iterm1)
                  do iJ=1,Atom%nJ(iterm)

                    ! Get level indexes
                    ilevel = Atom%irho(iterm)%irho_ij(iJ)
                    ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)

                    ! Only report up->down
                    if (ilevel.ge.ilevel1) cycle

                    ! If there is a rate
                    if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)) &
                        .gt.0d0) then

                      ! If allowed
                      if (Atom%fcflag(ilevel1,ilevel).eq.0) then

                        ! Report
                        write(200,'("Allowed rate from level",1x,'// &
                                  'i4,1x,"to level",1x,i4)') &
                                  ilevel1,ilevel

                      ! If forbidden
                      else if (Atom%fcflag(ilevel1,ilevel).eq.1) then

                        ! Report
                        write(200,'("Forbidden rate from level",'// &
                                  '1x,i4,1x,"to level",1x,i4)') &
                                  ilevel1,ilevel

                      end if ! Allower/forbidden

                    ! There is no rate
                    else

                      ! Report
                      write(200,'("No rate from level",1x,i4,1x,'// &
                                '"to level",1x,i4)') ilevel1,ilevel

                    end if ! If there is a rate

                  end do ! Lower level
                end do ! Upper level
              end do ! Lower term
            end do ! Upper term

            ! Close report file
            close(200)

          end if ! Master


          !
          ! Compute level b-b collisional transitions from the term
          ! wise ones
          !

          ! If Multilevel atom
          if (Atom%ML) then

            ! For each term pair
            do ilevel1=2,Atom%nMulti
              do ilevel=1,ilevel1-1

                ! Skip if there is already a rate
                if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)).gt.0d0) &
                  cycle

                ! Get transition index
                icol = Atom%icol(ilevel1,ilevel)

                ! Check there is not a collision, skip
                if (icol.lt.1) cycle

                ! Copy rate
                Atom%CcoeffJ(ilevel1,ilevel,:) = &
                                         Atom%Ccoeff(ilevel1,ilevel,:)

                ! Check orbital momentum selection rule
                if (nint(abs(Atom%rLval(ilevel) - &
                             Atom%rLval(ilevel1))).gt.1.or. &
                    nint(Atom%rLval(ilevel) + &
                         Atom%rLval(ilevel1)).eq.0) then

                  ! Flag forbidden
                  Atom%fcflag(ilevel1,ilevel) = 1
                  Atom%fcflag(ilevel,ilevel1) = 1

                ! If allowed by orbital angular momentum
                else

                  ! Check total angular momentum selection rule
                  if(nint(abs(Atom%rJval(1,ilevel1)- &
                     Atom%rJval(1,ilevel))).gt.1.or. &
                     nint(Atom%rJval(1,ilevel1)+ &
                     Atom%rJval(1,ilevel)).eq.0) then

                    ! Flag forbidden
                    Atom%fcflag(ilevel1,ilevel) = 1
                    Atom%fcflag(ilevel,ilevel1) = 1

                  ! Passed selection rules
                  else

                    ! Not connected radiatively
                    if (Atom%irad(ilevel,ilevel1).lt.1) then

                      ! Assume same parity, forbidden
                      Atom%fcflag(ilevel1,ilevel) = 1
                      Atom%fcflag(ilevel,ilevel1) = 1

                    end if ! Allowed transition
                  end if ! Angular momentum rules
                end if ! Orbital selection rules

              end do ! Every lower level
            end do ! Every upper level

          ! If multi-term atom
          else

            ! For every pair of terms
            do iterm1=2,Atom%nMulti
              do iterm=1,iterm1-1

                ! Get transition index
                icol = Atom%icol(iterm1,iterm)

                ! Skip if there is no rate
                if (icol.lt.1) cycle

                ! Get quantum numbers
                rL1 = Atom%rLval(iterm1)
                rL = Atom%rLval(iterm1)
                S = Atom%Sval(iterm1)

                ! For every upper level (within term)
                do iJ1=1,Atom%nJ(iterm1)

                  ! Get level index and angular momentum
                  ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)
                  rJ1 = Atom%rJval(iJ1,iterm1)

                  ! For every lower level (within term)
                  do iJ=1,Atom%nJ(iterm)

                    ! Get level index and angular momentum
                    ilevel = Atom%irho(iterm)%irho_ij(iJ)
                    rJ = Atom%rJval(iJ,iterm)

                    ! If there is a rate already, skip
                    if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)).gt. &
                        0d0) cycle

                    ! Compute 6J symbol
                    W6 = fun6j(rL1,rL,1d0,rJ,rJ1,S,Flgsg)

                    ! If the symbol is zero, skip
                    if (abs(W6).le.0d0) cycle

                    ! Scale factors
                    W6 = (2d0*rL1+1d0)*(2d0*rJ+1d0)*W6*W6

                    ! Compute level to level rate
                    Atom%CcoeffJ(ilevel1,ilevel,:) = &
                                        Atom%Ccoeff(iterm1,iterm,:)*W6

                  end do ! Lower level
                end do ! Upper level
              end do ! Lower term
            end do ! Upper term

          end if ! Multi-level or multi-term


          !
          ! Calculate the up-down ionizing collisions and
          ! low-up exciting collisions
          !

          ! For each lower term
          do it=1,Atom%nMulti

            ! For each upper term
            do itt=it,Atom%nMulti

              ! If they are in different ions
              if (Atom%stage(it).ne.Atom%stage(itt)) then

                ! For each level in lower term
                do iJ=1,Atom%nJ(it)

                  ! For each level in upper term
                  do iJJ=1,Atom%nJ(itt)

                    ! Determine the level indexes
                    i = Atom%irho(it)%irho_ij(iJ)
                    ii = Atom%irho(itt)%irho_ij(iJJ)

                    ! If there is not a collisional rate associated
                    ! with the collisional ionization, skip
                    if (maxval(Atom%CcoeffJ(i,ii,:)).le.0d0) cycle

                    ! For each height
                    do iz=1,nz

                      ! Compute the recombination rate
                      Atom%CcoeffJ(ii,i,iz) = Atom%CcoeffJ(i,ii,iz)* &
                                              Atom%populte(i,iz)/ &
                                              (Atom%populte(ii,iz) + &
                                               VTINY)
                    end do ! Heights
                  end do ! J for upper level
                end do ! J for lower level

              ! If they are in the same ion
              else

                ! For each level in lower term
                do iJ=1,Atom%nJ(it)

                  ! For each level in upper term
                  do iJJ=1,Atom%nJ(itt)

                    ! If it is the same level, skip
                    if(it.eq.itt.and.iJ.ge.iJJ) cycle

                    ! Determine level indexes
                    i = Atom%irho(it)%irho_ij(iJ)
                    ii = Atom%irho(itt)%irho_ij(iJJ)

                    ! If there is not a collisional rate associated
                    ! to the de-excitation, skip
                    if (maxval(Atom%CcoeffJ(ii,i,:)).le.0d0) cycle

                    ! For each height
                    do iz=1,nz

                      ! Compute the collisional excitation rate
                      Atom%CcoeffJ(i,ii,iz) = Atom%CcoeffJ(ii,i,iz)* &
                                              Atom%populte(ii,iz)/ &
                                              (Atom%populte(i,iz) + &
                                               VTINY)

                    end do ! Heights
                  end do ! J for upper level
                end do ! J for lower level

                ! Check if there is a term-term collisional transition
                ! indexed
                icol = Atom%icol(it,itt)

                ! If there is not, skip
                if (icol.lt.1) cycle

                ! Reset population vectors
                nu = 0d0
                nl = 0d0

                ! For each level in upper term
                do iJJ=1,Atom%nJ(itt)

                  ! Get level index
                  i = Atom%irho(itt)%irho_ij(iJJ)

                  ! Accumulate the LTE population of the upper term
                  nu = nu + Atom%populte(i,:)

                end do ! Levels in upper term

                ! For each level in lower term
                do iJ=1,Atom%nJ(it)

                  ! Get level index
                  i = Atom%irho(it)%irho_ij(iJ)

                  ! Accumulate the LTE population of the lower term
                  nl = nl + Atom%populte(i,:)

                end do ! Levels in lower term

                ! Compute the LTE population ratio
                nu = nu/nl

                ! For each height
                do iz=1,nz

                  ! Determine the excitation rate
                  Atom%Ccoeff(it,itt,iz) = Atom%Ccoeff(itt,it,iz)* &
                                           nu(iz)
                end do ! Heights

              end if ! Ionization stage comparison

            end do ! Upper term
          end do ! Lower term

          ! If there are non-symmetric collisions
          if (associated(Atom%Ccoeff_special)) then

            ! Do while there is data
            do while (associated(Atom%Ccoeff_special))

              ! Initialize
              p_col => Atom%Ccoeff_special
              nullify(p_col_p)

              ! Go to the last one
              do while (associated(p_col%next))
                p_col_p => p_col
                p_col => p_col%next
              end do ! Forward navigation

              ! Add rate
              Atom%CcoeffJ(p_col%ifrom,p_col%ito,:) = p_col%C + &
                            Atom%CcoeffJ(p_col%ifrom,p_col%ito,2)

              ! Flag ionization
              Atom%fcflag(p_col%ifrom,p_col%ito) = 2
              Atom%fcflag(p_col%ito,p_col%ifrom) = 2

              ! Deallocate data
              MRAMc = MRAMc - 1d-6*sizeof(p_col%C)
              deallocate(p_col%C)

              ! Remove the last one
              if (associated(p_col_p)) then
                nullify(p_col_p%next)
                nullify(p_col_p)
              ! We are in the last one
              else
                nullify(Atom%Ccoeff_special)
                nullify(p_col)
              end if

              ! Memory count
              MRAMc = MRAMc - 12d-6

            end do ! There is data

          end if ! non-symmetric collisional data
        end if ! Active atom
      end if ! There are inelastic collisions

      ! Check if everything is fine
      call control

      return

      end subroutine Initcols

!#####################################################################
!#####################################################################
!#####################################################################

      end module initcols_mod
