      !> Collisions initialization
      module initcols_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     09/25/2019
!  Last version:
!     09/29/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.5 - The COL log file is only created by
!                             request (TdPA)
!
!     09/25/2023:    V3.0.4 - Changed name for COL file (TdPA)
!
!     10/26/2022:    V3.0.3 - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.2 - Initializing every pointer when
!                             initcols starts (TdPA)
!                           - Nullifying p_col after each use (TdPA)
!                           - Bugfix: Correctly freeing memory in
!                             the Tbox data structure (TdPA)
!                           - Bugfix: Correctly freeing memory in
!                             the Ccoeff_special data structure (TdPA)
!                           - Bugfix: Atom%Ccoeff must be freed in
!                             every case, and not only in 1Ds (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: The input data was being
!                             directly modified, but this is not
!                             suitable for non-1D runs (TdPA)
!                           - Bugfix: Can only free Ccoeff_special
!                             if in the pure 1D case (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o The deallocation of the input data
!                                and the writing of a file is
!                                conditioned by the type of run with
!                                the free and aout variables.
!                              o The damping of the terms is now
!                                allocated and initialized here.
!                              o The option to input explicit elastic
!                                collisions has been removed.
!                             (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Removed domain decomposition (TdPA)
!
!     09/28/2020:    V1.0.5 - Charge transfer collisions were using
!                             electron density instead of hydrogen
!                             number density (TdPA)
!
!     09/11/2020:    V1.0.4 - Correctly deallocate collisional data
!                             before nullifying the pointer (TdPA)
!
!     03/05/2020:    V1.0.3 - When writing the file with collisional
!                             information, now the space in elements
!                             with singular letters is avoided (TdPA)
!
!     11/19/2019:    V1.0.2 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     10/02/2019:    V1.0.1 - Bugfix: Missed one col_type in Atom
!                             that should reference Atom%inelas (TdPA)
!
!     09/25/2019:    V1.0.0 - First version (TdPA)
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
!    This subroutine initializes the populations
!
!    Initcols:
!      Computes collisional rates
!
!    setNSCcoeff:
!      Adds the collisional rates of the special collisions (charge
!      transfer)
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

      !> Initializes the atomic populations\n
      !!         Atom(Atom_class): Structure with the atomic data\n
      !!        Atom0(Atom_class): A copy of Atom\n
      !!   filename(character(:)): Name of the file to read, if any\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Flgsg(Fctsg_class): Structure with factorials and
      !!                           signs\n
      !!         keeplog(logical): Bool to specify if keeping the
      !!                           collisional log file\n
      !!          active(logical): Bool to specify if this atom is
      !!                           active or not
      subroutine Initcols(Atom,Atmo,folder,Flgsg,keeplog,active)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Fctsg_class), intent(in):: Flgsg
      type(Atom_class), intent(inout):: Atom
      character(len=500), intent(in):: folder
      logical, intent(in):: active, keeplog

      ! Local

      character(len=2):: c2dump
      logical:: lin,free,aout
      integer:: ios,iz,it,itt,i,ii,iJ,iJJ,icol,iJ1,ilevel,iterm,jj
      integer:: ilevel1,iterm1,up,low,stagl,stagu
      integer:: K,minK,maxK
      double precision:: rJ,rJ1,rL,rL1,S,S1,El,Eu,gl,gu,Dfreq2,W6
      double precision:: d1,d2,d3
      double precision, dimension(nz):: nu, nl
      double precision, dimension(:), allocatable:: p_pop, CulI, Culin

      ! Pointer
      type(tmp_col_box_class), pointer:: p_col, p_col_p
      type(Tbox_class), pointer:: p_T, p_T_p

      ! Initialize pointers
      nullify(p_col,p_col_p,p_T,p_T_p)

      !
      ! Condition to free input data and output ascii collisions
      !
      free = run_mode.eq.0
      aout = run_mode.eq.0.and.keeplog


      !
      ! Initialize atomic level damping
      !
      allocate(Atom%damp(Atom%nMulti,nz))
      Atom%damp = 0d0

      ! For each upper term
      do itt=2,Atom%nMulti

        ! For each lower term
        do it=1,itt-1

          ! If transition exists
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
        Atom%gk = 0d0

      end if

      ! If there are inputs
      if (Atom%ngk.ge.1) then

        ! For each entry
        do ii=1,Atom%ngk

          ! Read to which level corresponds and find term and sublevel
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

              ! Reproduce the fit
              do iz=1,nZ

                Atom%gk(iterm,iJ,iJ,K,iz) =  &
                                           d1*1d-9* &
                                           ((Atmo%T(iz)*2d-4)**d2)* &
                                           (d3**(Atmo%T(iz)*2d-4))* &
                                           sum(Atmo%nh(iz,1:5))*1d-8 &
                                          + Atom%gk(iterm,iJ,iJ,K,iz)
              end do

            end if

            ! If the multipole was 0, add this to the damping, making
            ! a weighted average
            if (K.eq.0) then

              do iz=1,nZ

                Atom%damp(iterm,iz) = 1d-8* &
                                  Atom%gk(iterm,iJ,iJ,0,iz)* &
                                  (2d0*Atom%rJval(iJ,iterm) + 1d0)/ &
                                  Atom%deg(iterm)/c/(4d0*PI) + &
                                  Atom%damp(iterm,iz)
              end do

            end if

          end do ! Input sub entry

          ! Deallocate
          if (free) deallocate(Atom%elas(ii)%datum)

        end do ! Input entry

        ! Deallocate
        if (free) deallocate(Atom%elas)

        ! Only if active atoms
        if (active) then

          ! Heuristic definition of non-diagonal terms
          do iterm=1,Atom%nMulti
            do iJ=1,Atom%nJ(iterm)
              do iJ1=1,Atom%nJ(iterm)

                if(iJ.eq.iJ1)cycle

                minK = nint(abs(Atom%rJval(iJ,iterm) - &
                                Atom%rJval(iJ1,iterm)))
                maxK = nint(Atom%rJval(iJ,iterm) + &
                            Atom%rJval(iJ1,iterm))

                do iz=1,nZ
                  do K=minK,maxK

                    ! If both levels have 0 rate, the cross rate is 0
                    if(Atom%gk(iterm,iJ,iJ,K,iz).lt.1d-100.and. &
                       Atom%gk(iterm,iJ1,iJ1,K,iz).lt.1d-100)then

                      Atom%gk(iterm,iJ,iJ1,K,iz) = 0d0

                    ! If one of the levels have 0 rate, the cross rate
                    ! is the non-zero one
                    else if(Atom%gk(iterm,iJ,iJ,K,iz).lt.1d-100.or. &
                            Atom%gk(iterm,iJ1,iJ1,K,iz).lt.1d-100)then

                      if(Atom%gk(iterm,iJ,iJ,K,iz).lt.1d-100)then

                        Atom%gk(iterm,iJ,iJ1,K,iz) = &
                                           Atom%gk(iterm,iJ1,iJ1,K,iz)

                      else

                        Atom%gk(iterm,iJ,iJ1,K,iz) = &
                                           Atom%gk(iterm,iJ,iJ,K,iz)

                      end if

                    ! If both levels have non-zero rate, average them
                    else

                      Atom%gk(iterm,iJ,iJ1,K,iz) = &
                              (Atom%gk(iterm,iJ,iJ,K,iz)* &
                               (2d0*Atom%rJval(iJ,iterm) + 1d0) + &
                               Atom%gk(iterm,iJ1,iJ1,K,iz)* &
                               (2d0*Atom%rJval(iJ1,iterm) +.1d0))/ &
                              (2d0*Atom%rJval(iJ,iterm) + 1d0 + &
                               2d0*Atom%rJval(iJ1,iterm) + 1d0)

                    end if

                  end do ! K
                end do ! heights

              end do ! iJ1
            end do ! iJ
          end do ! Term

        ! If passive atom, deallocate
        else

          deallocate(Atom%gk)

        end if ! Active atom
      end if ! There are elastic collisions


      !
      ! Calculate the down-up ionizing collisions and
      ! up-down exciting collisions
      !


      ! If active atom or there are collisions, allocate
      if (active.or.Atom%ncol.ge.1) then
        ! Allocations
        ! Indexing of collisions between terms
        allocate(Atom%icol(Atom%nMulti,Atom%nMulti))
        Atom%icol = 0
        ! Collisional rates for collisions between terms
        allocate(Atom%Ccoeff(Atom%nMulti,Atom%nMulti,nZ))
        Atom%Ccoeff = 0d0
        ! Collisional rates for collisions between levels
        allocate(Atom%CcoeffJ(Atom%nlevel,Atom%nlevel,nZ))
        Atom%CcoeffJ = 0d0
      end if

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

          end if

          ! Translate
          up = Atom%inelas(icol)%up
          low = Atom%inelas(icol)%low

          ! If between terms, get degeneracy
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

          ! If between levels, get degeneracy and flags
          else if (p_T%col_type.eq.1) then

            ! Get the S, L, and J of levels
            S = Atom%Sval(Atom%term(low))
            S1 = Atom%Sval(Atom%term(up))
            rL = Atom%rLval(Atom%term(low))
            rL1 = Atom%rLval(Atom%term(up))
            rJ = Atom%rJval(Atom%sublevel(low),Atom%term(low))
            rJ1 = Atom%rJval(Atom%sublevel(up),Atom%term(up))
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

              Atom%fcflag(up,low) = 2
              Atom%fcflag(low,up) = 2

            ! Check if they are forbidden through dipole selection
            ! rules
            else

              ! Check if same term
              if (Atom%term(low).eq.Atom%term(up)) then

                Atom%fcflag(up,low) = 1
                Atom%fcflag(low,up) = 1

              else

                if (nint(abs(S1-S)).gt.0) then

                  Atom%fcflag(up,low) = 1
                  Atom%fcflag(low,up) = 1

                ! If spin allowed, check orbital angular momentum
                else

                  if (nint(abs(rL - rL1)).gt.1.or. &
                      nint(rL + rL1).eq.0) then
                    Atom%fcflag(up,low) = 1
                    Atom%fcflag(low,up) = 1

                  ! If allowed by orbital angular momentum, check
                  ! angular momentum
                  else

                    if (nint(abs(rJ1 - rJ)).gt.1.or. &
                        nint(rJ1 + rJ).eq.0) then

                      Atom%fcflag(up,low) = 1
                      Atom%fcflag(low,up) = 1

                    else

                      ! Not connected radiatively, assume same
                      ! parity
                      if (Atom%irad(Atom%term(up), &
                                    Atom%term(low)).lt.1) then

                        Atom%fcflag(up,low) = 1
                        Atom%fcflag(low,up) = 1

                      end if ! Allowed transition
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

            if (Atom%inelas(icol)%col_type.eq.0) then
              c2dump = 'BE'
              p_pop = Atmo%ne
            else if (Atom%inelas(icol)%col_type.eq.2) then
              c2dump = 'BP'
              p_pop = Atmo%nh(:,6)
            else if (Atom%inelas(icol)%col_type.eq.3) then
              c2dump = 'BH'
              p_pop = Atmo%nh(:,1)
            end if

            ! If neutral or ion flag, factorize weights
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              Culin = Atom%inelas(icol)%Cul*gu

              ! If neutral
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then
                Culin = Culin/gl
                Culin = Culin/sqrt(p_T%temp)
              ! If ion
              else
                Culin = Culin*sqrt(p_T%temp)
              end if
            ! No flag
            else
              Culin = Atom%inelas(icol)%Cul
            end if

            ! Interpolate the table
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,CulI,nZ,p_T%flin,lin)
            if (lin.and.pid.eq.0) then
              if (p_T%col_type.eq.0) then
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                       ' # Inelastic collisional rate '// &
                       c2dump//' ',up,low,'between terms in atom', &
                       Atom%Element,'was negative with Spline '// &
                       'interpolation, did linear.'
              else
                write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                      ' # Inelastic collisional rate '// &
                      c2dump//' ',up,low,'between levels in atom', &
                      Atom%Element,'was negative with Spline '// &
                      'interpolation, did linear.'
              end if
              call verbose
            end if

            ! If neutral or ion flag, factorize weights
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              CulI = CulI/gu

              ! If neutral
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then
                CulI = CulI*gl
                CulI = CulI*sqrt(Atmo%T)
              ! If ion
              else
                CulI = CulI/sqrt(Atmo%T)
              end if
            end if

            ! Convertion to 10^8 s-1
            CulI = CulI*1d-8

            ! Terms
            if (p_T%col_type.eq.0) then

              ! Compute actual rate and add to the damping
              do iz=1,nZ

                Atom%Ccoeff(up,low,iz) = Atom%Ccoeff(up,low,iz) + &
                                         CulI(iz)*p_pop(iz)

                Atom%damp(up,iz) = Atom%damp(up,iz) + &
                                   1d-8*CulI(iz)*p_pop(iz)/c/4d0/PI
              end do

            ! Levels
            else if (p_T%col_type.eq.1) then

              ! Compute actual
              do iz=1,nZ

                Atom%CcoeffJ(up,low,iz) = CulI(iz)*p_pop(iz) + &
                                          Atom%CcoeffJ(up,low,iz)
              end do

            end if


          !
          ! b-f
          !
          else if (Atom%inelas(icol)%col_type.eq.1.or. &
                   Atom%inelas(icol)%col_type.eq.4) then

            if (Atom%inelas(icol)%col_type.eq.1) then
              c2dump = 'FE'
              p_pop = Atmo%ne
            else if (Atom%inelas(icol)%col_type.eq.4) then
              c2dump = 'FH'
              p_pop = Atmo%nh(:,1)
            end if

            ! If any flag for this table
            if (p_T%nion.ge.0) then

              ! Calculate difference of energies
              Dfreq2 = abs(Eu-El)*c2*1d4
              Culin = Atom%inelas(icol)%Cul* &
                      exp(Dfreq2/p_T%temp)/ &
                      sqrt(p_T%temp)
            else
              Culin = Atom%inelas(icol)%Cul
            end if

            ! Interpolate the table
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,CulI,nZ,p_T%flin,lin)
            if (lin.and.pid.eq.0) then
              write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '//c2dump//' ', &
                  low,up,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did linear.'
              call verbose
            end if

            ! If any flag for this table
            if (p_T%nion.ge.0) then
              CulI = CulI*sqrt(Atmo%T)/exp(Dfreq2/Atmo%T)
            end if

            ! Convertion to 10^8 s-1
            CulI = CulI*1d-8

            ! Compute actual rate
            do iz=1,nZ

              Atom%CcoeffJ(low,up,iz) = CulI(iz)*p_pop(iz) + &
                                        Atom%CcoeffJ(low,up,iz)
            end do

          !
          ! Charge transfer
          !
          else if (Atom%inelas(icol)%col_type.eq.5.or. &
                   Atom%inelas(icol)%col_type.eq.6) then

            ! Create the next box and point p_col to it
            if (associated(Atom%Ccoeff_special)) then
              p_col => Atom%Ccoeff_special
              do while (associated(p_col%next))
                p_col => p_col%next
              end do
              allocate(p_col%next)
              p_col => p_col%next
              nullify(p_col%next)
            else
              allocate(Atom%Ccoeff_special)
              p_col => Atom%Ccoeff_special
              nullify(p_col%next)
            end if

            ! Direction
            if (Atom%inelas(icol)%col_type.eq.5) then
              p_col%ifrom = up
              p_col%ito = low
              c2dump = 'c0'
              p_pop = Atmo%nh(:,1)
            else
              p_col%ifrom = low
              p_col%ito = up
              c2dump = 'c+'
              p_pop = Atmo%nh(:,6)
            end if

            ! Coefficient
            allocate(p_col%C(nZ))

            ! If neutral or ion flag, factorize weights
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              Culin = Atom%inelas(icol)%Cul*gu

              ! If neutral
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then
                Culin = Culin/gl
                Culin = Culin/sqrt(p_T%temp)
              ! If ion
              else
                Culin = Culin*sqrt(p_T%temp)
              end if
            else
              Culin = Atom%inelas(icol)%Cul
            end if

            ! Interpolate the table
            call colinter(p_T%temp,Culin,p_T%nTmp, &
                          Atmo%T,p_Col%C,nZ,p_T%flin,lin)
            if (lin.and.pid.eq.0) then
              write(umsg,'(A,1x,i4,1x,"-->",i4,1x,A,1x,A,1x,A)') &
                  ' # Inelastic collisional rate '// &
                  c2dump//' ',up,low,'between levels in atom', &
                  Atom%Element,'was negative with Spline '// &
                  'interpolation, did linear.'
              call verbose
            end if

            ! If neutral or ion flag, factorize weights
            if (p_T%nion.ge.0.and.p_T%nion.le.3) then

              p_Col%C = p_Col%C/gu

              ! If neutral
              if (p_T%nion.eq.0.or.p_T%nion.eq.2) then
                p_Col%C = p_Col%C*gl
                p_Col%C = p_Col%C*sqrt(Atmo%T)
              ! If ion
              else
                p_Col%C = p_Col%C/sqrt(Atmo%T)
              end if
            end if

            ! Convertion to 10^8 s-1
            p_Col%C = p_Col%C*1d-8

            ! Get rate
            p_Col%C = p_Col%C*p_pop
            p_Col%flag = 2

            ! Nullify pointer
            nullify(p_col)

          end if ! Type of collision

          ! Deallocate table data
          if (free) deallocate(Atom%inelas(icol)%Cul)

        end do ! Collisional rates

        ! Deallocate inelastic database
        if (free) deallocate(Atom%inelas)

        ! Deallocate Tbox
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
            deallocate(p_T%temp)

            ! Remove the last one
            if (associated(p_T_p)) then
              nullify(p_T)
              deallocate(p_T_p%next)
              nullify(p_T_p%next)
              nullify(p_T_p)
            else
              nullify(p_T)
              deallocate(Atom%Tbox)
              nullify(Atom%Tbox)
            end if

          end do ! Every collision

        ! Not freeing data, just nullify pointer
        else

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

        ! If multilevel, just copy
        if (Atom%ML) then

          icol = maxval(Atom%icol)
          do ilevel=1,Atom%nMulti-1
            do ilevel1=ilevel+1,Atom%nMulti

              ! If ionizing, skip
              if (Atom%stage(ilevel).ne.Atom%stage(ilevel1)) cycle

              ! If there is already a term-term collision, skip
              if (Atom%icol(ilevel,ilevel1).gt.0) cycle

              ! If this collision is forbidden, don't include it
              if (Atom%fcflag(ilevel1,ilevel).gt.0) cycle

              Atom%Ccoeff(ilevel1,ilevel,:) = &
                                        Atom%CcoeffJ(ilevel1,ilevel,:)

              ! Add to the damping parameter
              Atom%damp(ilevel1,1:nZ) = 1d-8* &
                                 Atom%Ccoeff(ilevel1,ilevel,1:nZ)/c/ &
                                 (4d0*PI) + Atom%damp(ilevel1,1:nZ)

              if (.not.active) cycle

              ! And index the new transition if it exists
              if (maxval(Atom%Ccoeff(ilevel1,ilevel,:)).gt.0d0) then

                icol = icol + 1

                Atom%icol(ilevel1,ilevel) = icol
                Atom%icol(ilevel,ilevel1) = icol

              end if

            end do
          end do

        ! If multi-term, do it proper
        else

          icol = maxval(Atom%icol)
          do iterm=1,Atom%nMulti-1
            do iterm1=iterm+1,Atom%nMulti

              ! If ionizing, skip
              if (Atom%stage(iterm).ne.Atom%stage(iterm1)) cycle

              ! If there is already a term-term collision, skip
              if (Atom%icol(iterm,iterm1).gt.0) cycle

              do iJ1=1,Atom%nJ(iterm1)
                do iJ=1,Atom%nJ(iterm)

                  ilevel = Atom%irho(iterm)%irho_ij(iJ)
                  ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)

                  ! If this collision is forbidden, don't include it
                  ! on in the average
                  if (Atom%fcflag(ilevel1,ilevel).gt.0) cycle

                  Atom%Ccoeff(iterm1,iterm,:) = &
                                    Atom%Ccoeff(iterm1,iterm,:) + &
                                  (2d0*Atom%rJval(iJ1,iterm1)+1d0)* &
                                  Atom%CcoeffJ(ilevel1,ilevel,:)/ &
                                  Atom%deg(iterm1)

                end do
              end do

              ! Add to the damping parameter
              Atom%damp(iterm1,1:nZ) = 1d-8* &
                                 Atom%Ccoeff(iterm1,iterm,1:nZ)/c/ &
                                (4d0*PI) + Atom%damp(iterm1,1:nZ)

              if (.not.active) cycle

              ! And index the new transition if it exists
              if (maxval(Atom%Ccoeff(iterm1,iterm,:)).gt.0d0) then

                icol = icol + 1

                Atom%icol(iterm1,iterm) = icol
                Atom%icol(iterm,iterm1) = icol

              end if

            end do
          end do

        end if ! MT or ML


        ! If passive, remove all collisional data
        if (.not.active) then

          deallocate(Atom%Ccoeff)
          deallocate(Atom%CcoeffJ)
          if (allocated(Atom%fcflag).and. &
              free) deallocate(Atom%fcflag)
          deallocate(Atom%icol)

          ! Do while there is data
          do while (associated(Atom%Ccoeff_special))

            ! Initialize
            p_col => Atom%Ccoeff_special

            ! Go to the last one
            do while (associated(p_col%next))
              p_col_p => p_col
              p_col => p_col%next
            end do ! Forward navigation

            ! Deallocate content
            deallocate(p_col%C)

            ! Remove the last one
            if (associated(p_col_p)) then
              nullify(p_col_p%next)
              nullify(p_col_p)
            else
              deallocate(Atom%Ccoeff_special)
              nullify(Atom%Ccoeff_special)
              nullify(p_col)
            end if

          end do ! There is data

        ! If active
        else

          ! Output file if active and master
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

                if (iterm.ge.iterm1) cycle
                if (Atom%stage(iterm).ne.Atom%stage(iterm1)) cycle

                ! Check if there is term collisional rate
                if (Atom%icol(iterm1,iterm).lt.1) then
                  write(200,'(A)') ' '
                  write(200,'("No rate from term",1x,i4,'// &
                            '1x,"to term",1x,i4)') iterm1,iterm
                else
                  write(200,'(A)') ' '
                  write(200,'("Yes rate from term",1x,i4,'// &
                            '1x,"to term",1x,i4)') iterm1,iterm
                end if

                do iJ1=1,Atom%nJ(iterm1)
                  do iJ=1,Atom%nJ(iterm)

                    ilevel = Atom%irho(iterm)%irho_ij(iJ)
                    ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)

                    if (ilevel.ge.ilevel1) cycle

                    if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)) &
                        .gt.0d0) then
                      if (Atom%fcflag(ilevel1,ilevel).eq.0) then
                        write(200,'("Allowed rate from level",1x,'// &
                                  'i4,1x,"to level",1x,i4)') &
                                  ilevel1,ilevel
                      else if (Atom%fcflag(ilevel1,ilevel).eq.1) then
                        write(200,'("Forbidden rate from level",'// &
                                  '1x,i4,1x,"to level",1x,i4)') &
                                  ilevel1,ilevel
                      end if
                    else
                      write(200,'("No rate from level",1x,i4,1x,'// &
                                '"to level",1x,i4)') ilevel1,ilevel
                    end if
                  end do
                end do
              end do
            end do

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

                if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)).gt.0d0) &
                  cycle

                icol = Atom%icol(ilevel1,ilevel)

                ! Check there is a collision
                if (icol.lt.1) cycle

                Atom%CcoeffJ(ilevel1,ilevel,:) = &
                                         Atom%Ccoeff(ilevel1,ilevel,:)

                if (nint(abs(Atom%rLval(ilevel) - &
                             Atom%rLval(ilevel1))).gt.1.or. &
                    nint(Atom%rLval(ilevel) + &
                         Atom%rLval(ilevel1)).eq.0) then
                  Atom%fcflag(ilevel1,ilevel) = 1
                  Atom%fcflag(ilevel,ilevel1) = 1

                ! If allowed by orbital angular momentum, check
                ! angular momentum
                else

                  if(nint(abs(Atom%rJval(1,ilevel1)- &
                     Atom%rJval(1,ilevel))).gt.1.or. &
                     nint(Atom%rJval(1,ilevel1)+ &
                     Atom%rJval(1,ilevel)).eq.0) then

                    Atom%fcflag(ilevel1,ilevel) = 1
                    Atom%fcflag(ilevel,ilevel1) = 1

                  else

                    ! Not connected radiatively, assume same
                    ! parity
                    if (Atom%irad(ilevel,ilevel1).lt.1) then

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

                icol = Atom%icol(iterm1,iterm)

                ! Check there is a collision
                if (icol.lt.1) cycle

                rL1 = Atom%rLval(iterm1)
                rL = Atom%rLval(iterm1)
                S = Atom%Sval(iterm1)

                ! For every pair of levels
                do iJ1=1,Atom%nJ(iterm1)
                  ilevel1 = Atom%irho(iterm1)%irho_ij(iJ1)
                  rJ1 = Atom%rJval(iJ1,iterm1)

                  do iJ=1,Atom%nJ(iterm)
                    ilevel = Atom%irho(iterm)%irho_ij(iJ)
                    rJ = Atom%rJval(iJ,iterm)

                    if (maxval(Atom%CcoeffJ(ilevel1,ilevel,:)).gt. &
                        0d0) cycle

                    W6 = fun6j(rL1,rL,1d0,rJ,rJ1,S,Flgsg)

                    if (abs(W6).le.0d0) cycle

                    W6 = (2d0*rL1+1d0)*(2d0*rJ+1d0)*W6*W6

                    Atom%CcoeffJ(ilevel1,ilevel,:) = &
                                        Atom%Ccoeff(iterm1,iterm,:)*W6

                  end do
                end do ! Every pair of J levels
              end do
            end do ! Every pair of terms

          end if ! ML or MT


          !
          ! Calculate the up-down ionizing collisions and
          ! low-up exciting collisions
          !

          ! For each pair of terms
          do it=1,Atom%nMulti-1
            do itt=it,Atom%nMulti

              ! If they are in different ions
              if (Atom%stage(it).ne.Atom%stage(itt)) then

                ! For each pair of FS levels
                do iJ=1,Atom%nJ(it)
                  do iJJ=1,Atom%nJ(itt)

                    ! Determine the level indexes
                    i = Atom%irho(it)%irho_ij(iJ)
                    ii = Atom%irho(itt)%irho_ij(iJJ)

                    ! If there is a collisional rate associated
                    ! with the collisional ionization
                    if (maxval(Atom%CcoeffJ(i,ii,:)).le.0d0) cycle

                    ! Compute for each height the collisional
                    ! recombination rate
                    do iz=1,nz
                      Atom%CcoeffJ(ii,i,iz) = Atom%CcoeffJ(i,ii,iz)* &
                                              Atom%populte(i,iz)/ &
                                              (Atom%populte(ii,iz) + &
                                               VTINY)
                    end do ! Heights
                  end do ! J for upper level
                end do ! J for lower level

              ! If they are in the same ion
              else

                ! For each pair of FS levels
                do iJ=1,Atom%nJ(it)
                  do iJJ=1,Atom%nJ(itt)

                    ! If it is the same level, skip
                    if(it.eq.itt.and.iJ.ge.iJJ) cycle

                    ! Determine level indexes
                    i = Atom%irho(it)%irho_ij(iJ)
                    ii = Atom%irho(itt)%irho_ij(iJJ)

                    ! If there is a collisional rate associated
                    ! to the de-excitation
                    if (maxval(Atom%CcoeffJ(ii,i,:)).le.0d0) cycle

                    ! Compute for each height the collisional
                    ! excitation rate
                    do iz=1,nz
                      Atom%CcoeffJ(i,ii,iz) = Atom%CcoeffJ(ii,i,iz)* &
                                              Atom%populte(ii,iz)/ &
                                              (Atom%populte(i,iz) + &
                                               VTINY)
                    end do ! Heights
                  end do ! J for upper level
                end do ! J for lower level

                ! Check if there is a term collisional transition
                ! indexed
                icol = Atom%icol(it,itt)
                if (icol.lt.1) cycle

                ! Reset population vectors
                nu = 0d0
                nl = 0d0

                ! Get the LTE population of the upper term
                do iJJ=1,Atom%nJ(itt)

                  i = Atom%irho(itt)%irho_ij(iJJ)
                  nu = nu + Atom%populte(i,:)

                end do

                ! Get the LTE population of the upper term
                do iJ=1,Atom%nJ(it)

                  i = Atom%irho(it)%irho_ij(iJ)
                  nl = nl + Atom%populte(i,:)

                end do

                ! Compute the quotient
                nu = nu/nl

                ! For each height determine the excitation rate
                do iz=1,nz

                  Atom%Ccoeff(it,itt,iz) = Atom%Ccoeff(itt,it,iz)* &
                                           nu(iz)
                end do

              end if ! Ionization stage comparison

            end do ! Upper term
          end do ! Lower term

          ! If there are non-symmetric collisions
          if (associated(Atom%Ccoeff_special)) then

            ! Do while there is data
            do while (associated(Atom%Ccoeff_special))

              ! Initialize
              p_col => Atom%Ccoeff_special

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
              deallocate(p_col%C)

              ! Remove the last one
              if (associated(p_col_p)) then
                nullify(p_col_p%next)
                nullify(p_col_p)
              else
                nullify(Atom%Ccoeff_special)
                nullify(p_col)
              end if

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
