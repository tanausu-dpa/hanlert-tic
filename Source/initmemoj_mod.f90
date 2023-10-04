      !> Flow control of the code
      module initmemoj_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     03/11/2021
!  Last version:
!     09/29/2022 V3.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.4 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!
!     08/30/2023:    V3.0.3 - Added LTE lines to the initialization
!                             of 3J symbols (TdPA)
!
!     10/26/2022:    V3.0.2 - Changed the storage structure of the
!                             rdip variable (TdPA)
!                           - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.1 - Added limitation to the height do loop,
!                             although I am not sure this is the
!                             right decision (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     03/11/2021:    V1.0.0 - First version (TdPA)
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
!  initmemoJ:
!    Initializes the J symbols so later we can call everything with
!  OpenMP without racing conditions
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use funnj_mod
      use parameters_mod , only: TINYB , TINYJS , TINYEV , k2f , &
                                 TINYCO
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Secondary main with that controls the execution flow.\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !! LTElines(LTEline_class): Structure with the LTE line data\n
      !!      Flgsg(Fctsg_class): Structure with factorials and
      !!                          signs\n
      !!    Bstrength(dfloat(:)): Magnetic field strength\n
      !!             lp(logical): If doing formal solution in this run
      subroutine initmemoJ(Atom,LTElines,Flgsg,Bstrength,lp)

      ! I/O
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), allocatable:: LTElines
      type(Fctsg_class), intent(inout):: Flgsg
      logical, intent(in):: lp
      double precision, dimension(:), intent(in):: Bstrength


      ! Local

      logical:: BfieldY,BfieldN,zflag,zconj,zRpermit,zCpermit,zTS,zTA
      logical:: zrelax,zupper,zlower,zion,zpermit,zJ,zJ_JJ,zdiag,zMK
      logical:: zJ1,zJ1_JJ1,zRpermitJ,zCpermitJ,zMK_J,zK,zK1KK,zK2KK
      logical:: zR,zQ,zconj1,zupperJ,zlowerJ,zpermitJ,zrelaxJ
      logical, dimension(:,:),allocatable:: lfill

      integer:: ia,iterm,iJ,ir,iJ1,ir1,K,iQ,i,i1,itterm,itran,icol
      integer:: iJJ,iir,iJJ1,iir1,KK,iQQ,ii,ii1,Kr,iQr,iterml,itermu
      integer:: jtran,Ktilde,iL,iL1,iU,iU1,Kmin,Kmax,ios
      integer:: mF,K1,iQ1,iQl,Kl,nMu,nMl,iMl,iMu,iMl1,iMu1,iJl,iJu
      integer:: Ku,kL1,iJl1,iJu1,kU1,kUb,iJub
      integer:: iQdif,Kaux,iraux,iz,kLb,iJlb,iMf,ip,ipp
      integer:: ip1,iJlb1,itermf,klb1,nMf

      double precision, parameter:: dSLJ = .25d0

      double precision:: rL,S,rJ,dFS,rJsgn,rK,Q,tJsgn,rLL,SS,rJJ,rJJ1
      double precision:: rJJsgn,rKK,QQ,tJJsgn,rKr,rLu,rKtilde,Clb1
      double precision:: rJl,rJu,rJu1,rJf,rK1,Q1,Ql,rKl,rJumax,rJlmax
      double precision:: rMl,rMu,rMl1,rMu1,CC,CC1,Cu,Cl,Cu1,Cl1,rJl1
      double precision:: Cub,rJub,Qdif,rKaux,rJlb,Clb,rJfmax,rMf,p,p1
      double precision:: pp,Bmax,Qr,rJ1,rJlb1,rLf
      double precision:: JSym


      ! Determine if magnetic, non-magnetic or both
      if (maxval(Bstrength).gt.TINYB) then
        BfieldY = .True.
      else
        BfieldY = .False.
      end if
      if (minval(Bstrength).le.TINYB) then
        BfieldN = .True.
      else
        BfieldN = .False.
      end if

      ! For each atom
      do ia=1,nA

        ! If formal solution
        if (lp) then

          ! SEE

          allocate(lfill(Atom(ia)%ndim,Atom(ia)%ndim))
          lfill = .True.

          ! For each term (row)
          do iterm=1,Atom(ia)%nMulti

            ! Get term quantities
            rL = Atom(ia)%rLval(iterm)
            S = Atom(ia)%Sval(iterm)

            ! For each level (row_a)
            do iJ=1,Atom(ia)%nJ(iterm)

              ! Get level momentum
              rJ = Atom(ia)%rJval(iJ,iterm)

              ! Get level index
              ir = Atom(ia)%irho(iterm)%irho_ij(iJ)

              ! We use the conjugation properties of the SE matrix
              ! (see zconj), so we can restrict the range of iJ1
              ! For each level (row_b)
              do iJ1=iJ,Atom(ia)%nJ(iterm)  ! use with zconj
             !do iJ1=1,Atom(ia)%nJ(iterm)
             !do iJ1=iJ,iJ              ! no J,J'

                ! Get level momentum
                rJ1 = Atom(ia)%rJval(iJ1,iterm)

                ! Get level index
                ir1 = Atom(ia)%irho(iterm)%irho_ij(iJ1)

                ! Get energy diference in the appropriate unnits
                dFS = k2f*(Atom(ia)%FSfreq(iJ,iterm) - &
                           Atom(ia)%FSfreq(iJ1,iterm))

                ! Get sign from the difference of J,J'
                rJsgn = Flgsg%sg(nint(rJ-rJ1))

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(iterm))

                  ! Get real number
                  rK = dble(K)

                  ! For each Q
                  do iQ=-K,K

                    ! Get real number
                    Q = dble(iQ)

                    ! Get sign from the difference of J,J' and Q
                    tJsgn = rJsgn*Flgsg%sg(iQ)

                    ! Get the indexes of the elements J,J' and J',J
                    i = Atom(ia)%irho(iterm)%Jrho(iJ1,iJ)%kq(iQ,K)
                    i1 = Atom(ia)%irho(iterm)%Jrho(iJ,iJ1)%kq(-iQ,K)

                    ! Flag of filled line
                    zflag = .True.

                    ! If J!=J' or Q!=0, use conjugation properties
                    zconj = (iJ.ne.iJ1).or.(iQ.ne.0)

                    ! For each term (column)
                    do itterm=1,Atom(ia)%nMulti

                      ! Get the term quantities
                      rLL = Atom(ia)%rLval(itterm)
                      SS = Atom(ia)%Sval(itterm)

                      ! Check if there is a permitted transition
                      ! between these terms
                      itran = Atom(ia)%irad(itterm,iterm)
                      zRpermit = Atom(ia)%irad(itterm,iterm).ne.0 &
                                 .and. &
                                 nint(abs(rL-rLL)).le.1.and. &
                                 nint(rL+rLL).gt.0.and. &
                                 nint(abs(S-SS)).lt.1

                      ! Check if there is a permitted collision
                      ! between these terms
                      icol = Atom(ia)%icol(itterm,iterm)
                      zCpermit = Atom(ia)%icol(itterm,iterm).ne.0 &
                                 .and. &
                                 nint(abs(rL-rLL)).le.1.and. &
                                 nint(rL+rLL).gt.0.and. &
                                 nint(abs(S-SS)).lt.1

                      ! Check how the column relates to the row
                      zrelax = itterm.eq.iterm
                      zupper = itterm.gt.iterm
                      zlower = itterm.lt.iterm

                      ! Check ionization
                      zion = abs(Atom(ia)%stage(iterm) - &
                                 Atom(ia)%stage(itterm)).gt.0

                      ! Check if there is any permitted transition
                      ! between these terms
                      zpermit = zRpermit.or.zCpermit

                      ! For each level (column_a)
                      do iJJ=1,Atom(ia)%nJ(itterm)

                        ! Get level momentum
                        rJJ = Atom(ia)%rJval(iJJ,itterm)

                        ! Get level index
                        iir = Atom(ia)%irho(itterm)%irho_ij(iJJ)

                        ! Check if J = J''
                        zJ = abs(rJ-rJJ).lt.dSLJ

                        ! Check if J and J'' differ in less than 1
                        zJ_JJ = nint(abs(rJ-rJJ)).le.1

                        ! For each level (column_b)
                        do iJJ1=iJJ,Atom(ia)%nJ(itterm)

                          ! Get level momentun
                          rJJ1 = Atom(ia)%rJval(iJJ1,itterm)

                          ! Get level index
                          iir1 = Atom(ia)%irho(itterm)%irho_ij(iJJ1)

                          ! Get sign given J'' and J'''
                          rJJsgn = Flgsg%sg(nint(rJJ-rJJ1))

                          ! Check if J' = J'''
                          zJ1 = abs(rJ1-rJJ1).lt.dSLJ

                          ! Check if J' and J''' differ in less than 1
                          zJ1_JJ1 = nint(abs(rJ1-rJJ1)).le.1

                          ! Determine if there is a permitted FS tran.
                          zRpermitJ = zRpermit.and.zJ_JJ.and.zJ1_JJ1

                          ! Determine if there is a permitted FS coll.
                          zCpermitJ = zCpermit.and.zJ_JJ.and.zJ1_JJ1

                          ! Determine if there is any kind of
                          ! permitted transition
                          zpermitJ = zRpermitJ.or.zCpermitJ

                          ! Check if one of the pairs of J differs in
                          ! one, but the other is the same
                          zMK_J = (zJ_JJ.and.zJ1).or.(zJ1_JJ1.and.zJ)

                          ! Spontaneous emission transfer rate
                          if (zupper.and.zRpermitJ) then
                            JSym = fun6j(rJ,rJ1,rK,rJJ1,rJJ,1d0,Flgsg)
                            JSym = fun6j(rLL,rL,1d0,rJ,rJJ,S,Flgsg)
                            JSym = fun6j(rLL,rL,1d0,rJ1,rJJ1,S,Flgsg)
                          end if

                          ! For each K'
                          do KK=nint(abs(rJJ-rJJ1)), &
                                min(nint(rJJ+rJJ1), &
                                    Atom(ia)%Kcut(itterm))

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

                            ! Check if stimulated transfer rate is to
                            ! be calculated from the flag and K values
                            zTS = stm.and.zK2KK

                            ! Check if absorption transfer rate is to
                            ! be calculated from the term order and K
                            ! values
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

                              ! Get the indexes of the elements J'',
                              ! J''' and J''',J''
                              ii = Atom(ia)%irho(itterm)% &
                                            Jrho(iJJ1,iJJ)%kq(iQQ,KK)
                              ii1 = Atom(ia)%irho(itterm)% &
                                            Jrho(iJJ,iJJ1)%kq(-iQQ,KK)

                              ! If J!=J' or Q!=0, use conjugation
                              ! properties
                              zconj1 = (iJJ.ne.iJJ1).or.(iQQ.ne.0)

        !
        ! Reset the Identation
        !

        !
        ! If we have not touched this element yet, proceed
        !
        if (lfill(i,ii)) then

          ! Skip first row
          if (i.ne.1) then

            ! Transition rates are nonvanishing only for radiatively
            ! permitted transitions between the pairs (iJ,iJ1) and
            ! (iJJ,iJJ1)
            if (zpermitJ) then

              ! Contribute the transition rate for emission
              ! (espontaneous and stimulated) from the upper terms
              if (zupper) then

                ! If there is a permitted radiative transition
                if(zRpermitJ)then

                  ! If stimulated emission transfer rate to calculate
                  if (zTS) then

                    !
                    ! Stimulated emission transfer rate

                    ! For each K radiation
                    do Kr=0,Atom(ia)%Krad(itran)

                      rKr = dble(Kr)

                      ! For each Q radiation
                      do iQr=-Kr,Kr

                        Qr = dble(iQr)

                        JSym = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)

                      end do

                      JSym = fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0, &
                                   rK,rKK,rKr,Flgsg)
                    end do

                    JSym = fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)
                    JSym = fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)

                  end if ! stimulated emission
                end if ! There is a permitted radiative transition

                ! If K and K' differ in 2 as maximum and there is a
                ! permitted collisional transition
                if (zK2KK.and.zCpermitJ) then

                  ! De-excitation collisional transfer rate
                  JSym = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)
                  JSym = fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0, &
                               rK,rKK,0d0,Flgsg)

                  JSym = fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)
                  JSym = fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)

                end if ! If K-K' <= 2 and there is a collisional tran.
              end if ! If column > row

              ! If absorption transfer rate to be calculated
              if (zTA) then

                ! If there is a radiative permitted transition
                if (zRpermitJ) then

                  !
                  ! Absorption transfer rate

                  ! For each K radiation
                  do Kr=0,Atom(ia)%Krad(itran)

                    rKr = dble(Kr)

                    if(nint(rKr+rKK).gt. &
                       max(Atom(ia)%Kcut(iterm), &
                           Atom(ia)%Kcut(itterm))) cycle

                    ! For each Q radiation
                    do iQr=-Kr,Kr

                      Qr = dble(iQr)

                      JSym = fun3j(rK,rKK,rKr,-Q,QQ,-Qr,Flgsg)

                    end do

                    JSym = fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0, &
                                 rK,rKK,rKr,Flgsg)
                  end do

                  JSym = fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)
                  JSym = fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)

                end if ! permitted radiative transition

                ! If there is a permitted collisional transition
                if (zCpermitJ) then

                  ! Excitation collisional transfer rate
                  JSym = fun3j(rK,rKK,0d0,-Q,QQ,0d0,Flgsg)
                  JSym = fun9j(rJ,rJJ,1d0,rJ1,rJJ1,1d0, &
                         rK,rKK,0d0,Flgsg)

                  JSym = fun6j(rL,rLL,1d0,rJJ,rJ,S,Flgsg)
                  JSym = fun6j(rL,rLL,1d0,rJJ1,rJ1,S,Flgsg)

                end if ! If there is a permitted collisional trans
              end if ! If absorption transfer rate to be calculated
            end if ! If there is any permitted transition

            ! If it is a diagonal element (in terms)
            if (zrelax) then

              ! If Q=Q'
              if (zQ) then

                ! If magnetic kernel to be calculated
                if (zMK) then

                  ! Multiterm
                  if (.not.Atom(ia)%ML) then

                    ! If zJ
                    if (zJ) then
                      JSym = fun6j(rK,rKK,1d0,rJJ1,rJ1,rJ,Flgsg)
                      JSym = fun6j(rJJ1,rJ1,1d0,S,S,rL,Flgsg)
                    end if

                    ! If zJ1
                    if (zJ1) then
                      JSym = fun6j(rK,rKK,1d0,rJJ,rJ,rJ1,Flgsg)
                      JSym = fun6j(rJ,rJJ,1d0,S,S,rL,Flgsg)
                    end if

                    ! If any
                    if (zJ.or.zJ1) &
                      JSym = fun3j(rK,rKK,1d0,-Q,Q,0d0,Flgsg)

                  end if ! Multilevel or multiterm
                end if ! magnetic kernel to be calculated
              end if ! Q=Q'

              ! If diagonal in one of the J pairs and K-K'<=2
              if (zR) then

                ! If stimulated emission
                if (stm) then

                  !
                  ! Stimulated emission relaxation rate

                  ! Each lower term
                  do iterml=1,iterm-1

                    ! Transition
                    jtran = Atom(ia)%irad(iterml,iterm)

                    ! Cycle
                    if (jtran.eq.0) cycle

                    ! Lower level L
                    rLl = Atom(ia)%rLval(iterml)

                    ! Radiation field K
                    do Kr=0,Atom(ia)%Krad(jtran)

                      rKr = dble(Kr)

                      ! Radiation field Q
                      do iQr=-Kr,Kr

                        Qr=dble(iQr)

                        JSym = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)

                      end do

                      if (zJ) then
                        JSym = fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)
                        JSym = fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)
                      end if

                      if (zJ1) then
                        JSym = fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)
                        JSym = fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg)
                      end if

                      JSym = fun6j(rL,rL,rKr,1d0,1d0,rLl,Flgsg)

                    end do
                  end do ! Lower term

                end if ! stimulated emission

                !
                ! Absorption relaxation rate

                ! Each upper term
                do itermu=iterm+1,Atom(ia)%nMulti

                  ! Transition
                  jtran = Atom(ia)%irad(itermu,iterm)

                  ! Cycle
                  if (jtran.eq.0) cycle

                  ! Upper level angular momentum
                  rLu = Atom(ia)%rLval(itermu)

                  ! K radiation
                  do Kr=0,Atom(ia)%Krad(jtran)

                    rKr = dble(Kr)

                    if(nint(rKr+rKK).gt. &
                       max(Atom(ia)%Kcut(iterm), &
                           Atom(ia)%Kcut(itermu))) cycle

                    ! Q radiation
                    do iQr=-Kr,Kr

                      Qr = dble(iQr)

                      JSym = fun3j(rK,rKK,rKr,Q,-QQ,Qr,Flgsg)

                    end do

                    if (zJ) then
                      JSym = fun6j(rL,rL,rKr,rJJ1,rJ1,S,Flgsg)
                      JSym = fun6j(rK,rKK,rKr,rJJ1,rJ1,rJ,Flgsg)
                    end if

                    if (zJ1) then
                      JSym = fun6j(rL,rL,rKr,rJJ,rJ,S,Flgsg)
                      JSym = fun6j(rK,rKK,rKr,rJJ,rJ,rJ1,Flgsg)
                    end if

                    JSym = fun6j(rL,rL,rKr,1d0,1d0,rLu,Flgsg)

                  end do
                end do ! Upper level

                ! Only if diagonal
                if (abs(rK-rKK).lt..4d0) then

                  !
                  ! De-excitation collisional transfer rate

                  ! Each lower term
                  do iterml=1,iterm-1

                    ! Cycle
                    if (Atom(ia)%icol(iterml,iterm).eq.0) cycle

                    ! Lower level angular momentum
                    rLl = Atom(ia)%rLval(iterml)

                    if (zJ) then
                      JSym = fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)
                      JSym = fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)
                    end if

                    if (zJ1) then
                      JSym = fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)
                      JSym = fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg)
                    end if

                    JSym = fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)
                    JSym = fun6j(rL,rL,0d0,1d0,1d0,rLl,Flgsg)

                  end do ! Lower terms

                  !
                  ! Excitation collisional transfer rate

                  ! For each upper term
                  do itermu=iterm+1,Atom(ia)%nMulti

                    ! Cycle
                    if (Atom(ia)%icol(itermu,iterm).eq.0) cycle

                    ! Upper level angular momentum
                    rLu = Atom(ia)%rLval(itermu)

                    if (zJ) then
                      JSym = fun6j(rL,rL,0d0,rJJ1,rJ1,S,Flgsg)
                      JSym = fun6j(rK,rKK,0d0,rJJ1,rJ1,rJ,Flgsg)
                    end if

                    if (zJ1) then
                      JSym = fun6j(rL,rL,0d0,rJJ,rJ,S,Flgsg)
                      JSym = fun6j(rK,rKK,0d0,rJJ,rJ,rJ1,Flgsg)
                    end if

                    JSym = fun3j(rK,rKK,0d0,Q,-QQ,0d0,Flgsg)
                    JSym = fun6j(rL,rL,0d0,1d0,1d0,rLu,Flgsg)

                  end do ! Upper terms

                end if ! Diagonal

              end if ! if diagonal in one pair of J and K-K'<=2
            end if ! If it is a diagonal element (in terms)

            ! If diagonal in K, Q, and in rhoKQ coherences
            if (zK.and.zQ.and.(iJ.eq.iJ1).and.(iJJ.eq.iJJ1)) then

              ! Check relative positions of levels
              zrelaxJ = iir.eq.ir

              !
              ! Multi-level contributions
              !

              ! If not relaxation rate
              if (.not.zrelaxJ) then

                ! Check relative positions of levels
                zupperJ = iir.gt.ir
                zlowerJ = iir.lt.ir

                ! Check there are collisions
                if (maxval(Atom(ia)%CcoeffJ(iir,ir,:)).gt.0d0) then

                  ! Check there is a forbidden collision
                  if (Atom(ia)%fcflag(iir,ir).gt.0.and. &
                      (K.eq.0.or..not.zion)) then

                    ! If column > row
                    if (zupperJ) then

                      !
                      ! De-excitation collisional transfer rate
                      if (K.gt.0.and.fcol_transfer) then

                        ! Get K tilde
                        Ktilde = nint(abs(rJ - rJJ))

                        ! Ktilde >= 1
                        if (Ktilde.lt.1) Ktilde = 1
                        rKtilde = dble(Ktilde)

                        JSym = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

                        ! If invalid 6J
                        if (abs(JSym).ge.TINYJS) &
                          Jsym = fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)
                      end if

                    ! If column < row
                    else if (zlowerJ) then

                      !
                      ! Excitation collisional transfer rate
                      if (K.gt.0.and.fcol_transfer) then

                        ! Get K tilde
                        Ktilde = nint(abs(rJ - rJJ))

                        ! Ktilde >= 1
                        if (Ktilde.lt.1) Ktilde = 1
                        rKtilde = dble(Ktilde)

                        JSym = fun6j(rJ,rJ,0d0,rJJ,rJJ,rKtilde,Flgsg)

                        ! If invalid 6J
                        if (abs(JSym).ge.TINYJS) &
                          JSym = fun6j(rJ,rJ,rK,rJJ,rJJ,rKtilde,Flgsg)
                      end if
                    end if ! Direction of transition
                  end if ! Forbidden transfer rate
                end if ! Collisions
              end if ! Relaxation rate
            end if ! Level diagonal element
          end if ! first row?
        end if ! Filled element

                              !
                              ! Recover the identation
                              !

                            end do ! Q'
                          end do ! K'
                        end do ! J''' (column_b)
                      end do ! J'' (column_a)
                    end do ! terms (column)
                  end do ! Q
                end do ! K
              end do ! J' (row_b)
            end do ! J (row_a)
          end do ! Terms (row)

          ! Deallocate filler
          deallocate(lfill)

        end if ! Formal solution

        !
        ! RT Coefficients
        !

        ! Magnetic field
        if (BfieldY) then

          ! Initialize search
          Bmax = 0
          iz = -1

          ! For each height
          do ios=Rz0,Rz1

            ! If no field, skip
            if (Bstrength(ios).lt.TINYB) cycle

            ! If larger field, get representative
            if (Bstrength(ios).gt.Bmax) then
              Bmax = Bstrength(ios)
              iz = ios
            end if

          end do ! Heights

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Angular momentum
            S = Atom(ia)%Sval(itermu)
            rLu = Atom(ia)%rLval(itermu)
            rLl = Atom(ia)%rLval(iterml)

            !
            ! emiss
            !

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu+S
      nMu = nint(2d0*rJumax+1d0)
      rJlmax = rLl + S
      nMl = nint(2d0*rJlmax+1d0)

      !
      ! Compute emissivity
      !

      ! For each Ml
      do iMl=1,nMl

        ! Value of Ml
        rMl = -rJlmax + dble(iMl-1)

        ! For each mu_l
        do iL=1,Atom(ia)%nblk(iMl,iterml)

          ! For each Mu
          do iMu=1,nMu

            ! Value of Mu
            rMu = -rJumax + dble(iMu-1)

            ! If not pi nor sigma, skip
            if (nint(abs(rMu-rMl)).gt.1) cycle

            ! Get difference between M momentums in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_u
            do iU=1,Atom(ia)%nblk(iMu,itermu) ! sum over mu_u

              ! For each Mu'
              do iMu1=1,nMu

                ! Value of Mu'
                rMu1 = -rJumax + dble(iMu1-1)

                ! If not pi nor sigma, skip
                if (nint(abs(rMu1-rMl)).gt.1) cycle

                ! Get the difference between M momentums
                q1 = rMu1-rMl
                QQ = q1-q

                ! Make the difference integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! For each K
                do K=abs(iQQ),2

                  ! Get the real number
                  rK = dble(K)

                  ! Racah algebra
                  JSym = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

                  ! If not allowed (3j-sym=0) skip
                  if (abs(JSym).lt.TINYJS) cycle

                  ! For each Jl
                  do kL=1,Atom(ia)%nblk(iMl,iterml)

                    ! Get eigenvector lower level
                    Cl = Atom(ia)%evec(kL,iL,iMl,iterml,iz)

                    ! If coefficient too small, skip
                    if (abs(Cl).lt.TINYEV) cycle

                    ! Get J level index
                    iJl = Atom(ia)%iJval(kL,iMl,iterml)

                    ! Get angular momentum
                    rJl = Atom(ia)%rJval(iJl,iterml)

                    ! For each Ju
                    do kU=1,Atom(ia)%nblk(iMu,itermu)

                      ! Get eigenvector lower level
                      Cu = Atom(ia)%evec(kU,iU,iMu,itermu,iz)

                      ! If coefficient too small, skip
                      if (abs(Cu).lt.TINYEV) cycle

                      ! Get J level index
                      iJu = Atom(ia)%iJval(kU,iMu,itermu)

                      ! Get angular momentum
                      rJu = Atom(ia)%rJval(iJu,itermu)

                      ! Coefficients times the dipolar matrix
                      CC = Cu*Cl*Atom(ia)%rdip(jtran)% &
                                 rdip(iq,iMu,iMl,iJu,iJl)

                      ! If coefficient small, skip
                      if (abs(CC).lt.TINYCO) cycle

                      ! For each Jl'
                      do kL1=1,Atom(ia)%nblk(iMl,iterml)

                        ! Get eigenvector upper level'
                        Cl1 = Atom(ia)%evec(kL1,iL,iMl,iterml,iz)

                        ! If coefficient too small, skip
                        if (abs(Cl1).lt.TINYEV) cycle

                        ! Get J level index
                        iJl1 = Atom(ia)%iJval(kL1,iMl,iterml)

                        ! Get angular momentum
                        rJl1 = Atom(ia)%rJval(iJl1,iterml)

                        ! For each Ju'
                        do kU1=1,Atom(ia)%nblk(iMu1,itermu)

                          ! Get J level index
                          iJu1 = Atom(ia)%iJval(kU1,iMu1,itermu)

                          ! Get angular momentum
                          rJu1 = Atom(ia)%rJval(iJu1,itermu)

                          ! Coefficient times dipolar matrix
                          CC1 = Cl1* &
                                Atom(ia)%rdip(jtran)% &
                                         rdip(iq1,iMu1,iMl,iJu1,iJl1)

                          ! If coefficient small, skip
                          if (abs(CC1).lt.TINYCO) cycle

                          ! For each Jub
                          do kUb=1,Atom(ia)%nblk(iMu,itermu) ! sum Jub

                            ! Get eigenvector for upper level b
                            Cub = Atom(ia)%evec(kUb,iU,iMu,itermu,iz)

                            ! If coefficient too small, skip
                            if (abs(Cub).lt.TINYEV) cycle

                            ! Get J level index
                            iJub = Atom(ia)%iJval(kUb,iMu,itermu)

                            ! Get angular momentum
                            rJub = Atom(ia)%rJval(iJub,itermu)

                            !
                            ! Sum over (Kl,Ql)

                            ! Difference between magnetic momentums
                            Qdif = rMu1-rMu

                            ! Convert to integer
                            iQdif = nint(Qdif)

                            ! Determine the limits in K
                            Kmin = max(abs(iQdif), &
                                       nint(abs(rJu1-rJub)))
                            Kmax = min(nint(rJu1+rJub), &
                                       Atom(ia)%Kcut(itermu))

                            ! For each K
                            do Kaux=Kmin,Kmax

                              ! Get the real number
                              rKaux = dble(Kaux)

                              ! Get the SEE index
                              iRaux = Atom(ia)%irho(itermu)% &
                                               Jrho(iJub,iJu1)% &
                                               kq(iQdif,Kaux)

                              ! If flagged as small, skip
                              if (iRaux.le.0) cycle

                              ! Racah algebra
                              JSym = fun3j(rJu1,rJub,rKaux, &
                                           rMu1,-rMu,-Qdif,Flgsg)
                            end do ! K
                          end do ! kUb
                        end do ! kU1
                      end do ! kL1
                    end do ! kU
                  end do ! kL
                end do ! K
              end do ! Mu1
            end do ! iU
          end do ! Mu
        end do ! iL
      end do ! Ml

            !
            ! absorb
            !

      ! For each Mu
      do iMu=1,nMu

        ! Value of Mu
        rMu = -rJumax + dble(iMu-1)

        ! For each mu_u
        do iU=1,Atom(ia)%nblk(iMu,itermu)

          ! For each Ml
          do iMl=1,nMl

            ! Value of Ml
            rMl = -rJlmax + dble(iMl-1)

            ! If not pi nor sigma, skip
            if (nint(abs(rMu-rMl)).gt.1) cycle

            ! Get difference between M momentums in integer
            q = rMu-rMl
            iq = nint(q)

            ! For each mu_l
            do iL=1,Atom(ia)%nblk(iMl,iterml)

              ! For each Ml'
              do iMl1=1,nMl

                ! Value of Ml'
                rMl1 = -rJlmax + dble(iMl1-1)

                ! If not pi nor sigma, skip
                if (nint(abs(rMu-rMl1)).gt.1) cycle

                ! Get the difference between M momentums
                q1 = rMu - rMl1
                QQ = q1-q

                ! Make the difference integers
                iq1 = nint(q1)
                iQQ = nint(QQ)

                ! For each K
                do K=abs(iQQ),2

                  ! Get the real number
                  rK = dble(K)

                  ! Racah algebra
                  JSym = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

                  ! If not allowed (3j-sym=0) skip
                  if (abs(JSym).lt.TINYJS) cycle

                  ! For each Ju
                  do kU=1,Atom(ia)%nblk(iMu,itermu)

                    ! Get eigenvector upper level
                    Cu = Atom(ia)%evec(kU,iU,iMu,itermu,iz)

                    ! If coefficient too small, skip
                    if (abs(Cu).lt.TINYEV) cycle

                    ! Get J level index
                    iJu = Atom(ia)%iJval(kU,iMu,itermu)

                    ! Get angular momentum
                    rJu = Atom(ia)%rJval(iJu,itermu)

                    ! For each Jl
                    do kL=1,Atom(ia)%nblk(iMl,iterml)

                      ! Get eigenvector lower level
                      Cl = Atom(ia)%evec(kL,iL,iMl,iterml,iz)

                      ! If coefficient too small, skip
                      if (abs(Cl).lt.TINYEV) cycle

                      ! Get J level index
                      iJl = Atom(ia)%iJval(kL,iMl,iterml)

                      ! Get angular momentum
                      rJl = Atom(ia)%rJval(iJl,iterml)

                      ! Coefficients times the dipolar matrix
                      CC = Cl*Cu*Atom(ia)%rdip(jtran)% &
                                 rdip(iq,iMu,iMl,iJu,iJl)

                      ! If coefficient small, skip
                      if (abs(CC).lt.TINYEV) cycle

                      ! For each Ju'
                      do kU1=1,Atom(ia)%nblk(iMu,itermu)

                        ! Get eigenvector upper level'
                        Cu1 = Atom(ia)%evec(kU1,iU,iMu,itermu,iz)

                        ! If coefficient too small, skip
                        if (abs(Cu1).lt.TINYEV) cycle

                        ! Get J level index
                        iJu1 = Atom(ia)%iJval(kU1,iMu,itermu)

                        ! Get angular momentum
                        rJu1 = Atom(ia)%rJval(iJu1,itermu)

                        ! For each Jl'
                        do kL1=1,Atom(ia)%nblk(iMl1,iterml)

                          ! Get J level index
                          iJl1 = Atom(ia)%iJval(kL1,iMl1,iterml)

                          ! Get angular momentum
                          rJl1 = Atom(ia)%rJval(iJl1,iterml)

                          ! Coefficient times dipolar matrix
                          CC1 = Cu1*Atom(ia)%rdip(jtran)% &
                                    rdip(iq1,iMu,iMl1, &
                                         iJu1,iJl1)

                          ! If coefficient small, skip
                          if (abs(CC1).lt.TINYCO) cycle

                          ! For each Jlb
                          do kLb=1,Atom(ia)%nblk(iMl,iterml)

                            ! Get eigenvector for lower level b
                            Clb = Atom(ia)%evec(kLb,iL,iMl,iterml,iz)

                            ! If coefficient too small, skip
                            if (abs(Clb).lt.TINYEV) cycle

                            ! Get J level index
                            iJlb = Atom(ia)%iJval(kLb,iMl,iterml)

                            ! Get angular momentum
                            rJlb = Atom(ia)%rJval(iJlb,iterml)

                            ! Coefficient and sign
                            Clb = Clb*Flgsg%sg(nint(rJlb-rMl))

                            ! Difference between magnetic momentums
                            Qdif = rMl-rMl1

                            ! Convert to integer
                            iQdif = nint(Qdif)

                            ! Determine the limits in K
                            Kmin = max(abs(iQdif), &
                                       nint(abs(rJlb-rJl1)))
                            Kmax = min(nint(rJlb+rJl1), &
                                       Atom(ia)%Kcut(iterml))

                            ! For each K
                            do Kaux=Kmin,Kmax

                              ! Check for Kcut
                              if (KcutAB.and.abs(Kaux).gt.Kcut) cycle

                              ! Get the real number
                              rKaux = dble(Kaux)

                              ! Get the SEE index
                              iRaux = Atom(ia)%irho(iterml)% &
                                            Jrho(iJl1,iJlb)% &
                                            kq(iQdif,Kaux)

                              ! If flagged as small, skip
                              if (iRaux.le.0) cycle

                              ! Racah algebra
                              JSym = fun3j(rJlb,rJl1,rKaux, &
                                           rMl,-rMl1,-Qdif,Flgsg)
                            end do ! K
                          end do ! kLb
                        end do ! kL1
                      end do ! kU1
                    end do ! kL
                  end do ! kU
                end do ! K
              end do ! Ml1
            end do ! iL
          end do ! Ml
        end do ! iU
      end do ! Mu


            !
            ! emiss2ord
            !

      ! Final term
      itermf = iterml
      rLf = Atom(ia)%rLval(itermf)

      ! Determine the maximum angular momentum and the number
      ! of magnetic sublevels for that maximum momentum
      rJumax = rLu + S
      nMu = nint(2d0*rJumax+1d0)
      rJfmax = rLf + S
      nMf = nint(2d0*rJfmax+1d0)

      ! For all the possible lower terms
      do i=1,Atom(ia)%nMulti-1

        ! If there is no transition or this term is larger
        ! than the upper term of the output transition, skip
        if(i.ge.itermu.or.Atom(ia)%irad(i,itermu).eq.0) cycle

        ! Store the input lower term index
        iterml = i

        ! Get index of input transition
        itran = Atom(ia)%irad(iterml,itermu)

        ! Angular momentum input lower level
        rLl = Atom(ia)%rLval(iterml)

        ! Determine maximum value of J and number of magnetic
        ! sublevels for this maximum J
        rJlmax = rLl + S
        nMl = nint(2d0*rJlmax+1d0)

        ! For each Mf
        do iMf=1,nMf

          ! Value of Mf
          rMf = -rJfmax + dble(iMf-1)

          ! For each mu_f
          do mF=1,Atom(ia)%nblk(iMf,itermf)

            ! For each Mu
            do iMu=1,nMu

              ! Value of Mu
              rMu = -rJumax + dble(iMu-1)

              ! Difference between M momentums, done integer
              q = rMu - rMf
              iq = nint(q)

              ! If not pi nor sigma, skip
              if(abs(iq).gt.1) cycle

              ! For each mu_u
              do iU=1,Atom(ia)%nblk(iMu,itermu)

                ! For each Mu'
                do iMu1=1,nMu

                  ! Value of Mu'
                  rMu1 = -rJumax + dble(iMu1-1)

                  ! Difference between M momentums
                  q1 = rMu1-rMf
                  QQ = q1-q

                  ! Convert to integers
                  iq1 = nint(q1)
                  iQQ = nint(QQ)

                  ! If not pi or sigma, skip
                  if(abs(iq1).gt.1) cycle

                  ! For each mu_u'
                  do iU1=1,Atom(ia)%nblk(iMu1,itermu)

                    ! For each Ml
                    do iMl=1,nMl

                      ! Value of Ml
                      rMl = -rJlmax + dble(iMl-1)

                      ! Difference between M momentums, in integer
                      p = rMu-rMl
                      ip = nint(p)

                      ! If not pi nor sigma, skip
                      if(abs(ip).gt.1) cycle

                      ! For each mu_l
                      do iL=1,Atom(ia)%nblk(iMl,iterml)

                        ! For each Ml'
                        do iMl1=1,nMl

                          ! Value of Ml'
                          rMl1 = -rJlmax + dble(iMl1-1)

                          ! Difference between M momentums
                          p1 = rMu1-rMl1
                          PP = p1-p

                          ! Convert to integer
                          ip1 = nint(p1)
                          iPP = nint(PP)

                          ! If not pi nor sigma, skip
                          if(abs(ip1).gt.1) cycle

                          ! For each mu_l'
                          do iL1 = 1,Atom(ia)%nblk(iMl1,iterml)

        !
        ! Reset indexing
        !

        ! For each K'
        do K1=abs(iPP),2

          ! Get real value
          rK1 = dble(K1)

          ! For each K
          do K=abs(iQQ),2

            ! Get real value
            rK = dble(K)

            ! Racah algebra
            JSym = fun3j(1d0,1d0,rK,-q,q1,-QQ,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(JSym).lt.TINYJS) cycle

            ! Racah algebra
            JSym = fun3j(1d0,1d0,rK1,-p,p1,-PP,Flgsg)

            ! If forbidden (3j-sym=0), skip
            if(abs(JSym).lt.TINYJS) cycle

            ! For each Jlb
            do kLb=1,Atom(ia)%nblk(iMl,iterml)

              ! Get eigenvector for lower b level
              Clb = Atom(ia)%evec(kLb,iL,iMl,iterml,iz)

              ! If coefficient too small, skip
              if(abs(Clb).lt.TINYEV) cycle

              ! Get J level index
              iJlb = Atom(ia)%iJval(kLb,iMl,iterml)

              ! Get angular momentum
              rJlb = Atom(ia)%rJval(iJlb,iterml)

              ! For each Jlb'
              do kLb1=1,Atom(ia)%nblk(iMl1,iterml)

                ! Get eigenvector for lower b level
                Clb1 = Atom(ia)%evec(kLb1,iL1,iMl1,iterml,iz)

                ! If coefficient too small, skip
                if(abs(Clb1).lt.TINYEV) cycle

                ! Get J level index
                iJlb1 = Atom(ia)%iJval(kLb1,iMl1,iterml)

                ! Get angular momentum
                rJlb1 = Atom(ia)%rJval(iJlb1,iterml)

                ! sum over (Kl,Ql)

                ! Difference between magnetic momentums
                Qdif = rMl-rMl1

                ! Convert to integer
                iQdif = nint(Qdif)

                ! Determine the limits in K
                Kmin = max(abs(iQdif),nint(abs(rJlb-rJlb1)))
                Kmax = min(nint(rJlb+rJlb1),Kcut)

                ! For each K
                do Kaux=Kmin,Kmax

                  ! Check for Kcut
                  if (KcutAB.and.abs(Kaux-K1).gt.Kcut) cycle

                  ! Get the real number
                  rKaux = dble(Kaux)

                  ! Get the SEE index
                  iRaux = Atom(ia)%irho(iterml)% &
                                   Jrho(iJlb1,iJlb)% &
                                   kq(iQdif,Kaux)

                  ! If flagged as small, skip
                  if (iRaux.le.0) cycle

                  ! Racah algebra
                  JSym = fun3j(rJlb,rJlb1,rKaux, &
                               rMl,-rMl1,-Qdif,Flgsg)

                end do ! K
              end do ! kLb1
            end do ! kLb
          end do ! K1
        end do ! K
                          end do ! iL1
                        end do ! Ml1
                      end do ! iL
                    end do ! Ml
                  end do ! iU1
                end do ! Mu1
              end do ! iU
            end do ! Mu
          end do ! mF
        end do ! Mf
      end do ! Terms

          end do ! Transitions

        end if ! Yes magnetic field

        ! No magnetic field
        if (BfieldN) then

          ! Initialize search
          Bmax = 0
          iz = -1

          ! For each height
          do ios=1,nz

            ! If Yes field, skip
            if (Bstrength(ios).gt.TINYB) cycle

            iz = ios
            exit

          end do ! Heights

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Quantum numbers
            S = Atom(ia)%Sval(itermu)
            rLu = Atom(ia)%rLval(itermu)
            rLl = Atom(ia)%rLval(iterml)

            !
            ! emissNB
            !

            ! For each Jl
            do iL=1,Atom(ia)%nJ(iterml)

              ! Jl
              rJl = Atom(ia)%rJval(iL,iterml)

              ! For each Ju
              do iU=1,Atom(ia)%nJ(itermu)

                ! Ju
                rJu = Atom(ia)%rJval(iU,itermu)

                ! 6-j
                JSym = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

                ! Check size
                if (abs(JSym).lt.TINYJS) cycle

                ! For each Ju'
                do iU1=1,Atom(ia)%nJ(itermu)

                  ! Ju'
                  rJu1 = Atom(ia)%rJval(iU1,itermu)

                  ! 6-j
                  JSym = fun6j(rLu,rLl,1d0,rJl,rJu1,S,Flgsg)

                  ! Check size
                  if (abs(JSym).lt.TINYJS) cycle

                  ! Determine the limits in K
                  Kmin = nint(abs(rJu-rJu1))
                  Kmax = min(nint(rJu+rJu1),Atom(ia)%Kcut(itermu),2)

                  ! For each K
                  do K=Kmin,Kmax

                    ! Get the real number
                    rK = dble(K)

                    ! For each Q
                    do iQ=-K,K

                      ! Get the SEE index
                      iR = Atom(ia)%irho(itermu)% &
                                    Jrho(iU,iU1)%kq(iQ,K)

                      ! If flagged as small, skip
                      if (iR.le.0) cycle

                      ! 6-j
                      JSym = fun6j(1d0,1d0,rK,rJu,rJu1,rJl,Flgsg)

                    end do ! Q
                  end do ! K
                end do ! iU1
              end do ! iU
            end do ! iL

            !
            ! absorbNB
            !

            ! For each Ju
            do iU=1,Atom(ia)%nJ(itermu)

              ! Get Ju
              rJu = Atom(ia)%rJval(iU,itermu)

              ! For each mu_l
              do iL=1,Atom(ia)%nJ(iterml)

                ! Get Jl
                rJl = Atom(ia)%rJval(iL,iterml)

                ! 6-j
                JSym = fun6j(rLu,rLl,1d0,rJl,rJu,S,Flgsg)

                ! Check size
                if (abs(JSym).lt.TINYJS) cycle

                ! For each mu_l
                do iL1=1,Atom(ia)%nJ(iterml)

                  ! Get Jl
                  rJl1 = Atom(ia)%rJval(iL1,iterml)

                  ! 6-j
                  JSym = fun6j(rLu,rLl,1d0,rJl1,rJu,S,Flgsg)

                  if (abs(JSym).lt.TINYJS) cycle

                  ! Determine the limits in K
                  Kmin = nint(abs(rJl-rJl1))
                  Kmax = min(nint(rJl+rJl1),Atom(ia)%Kcut(iterml),2)

                  ! For each K
                  do K=Kmin,Kmax

                    ! Get the real number
                    rK = dble(K)

                    ! For each Q
                    do iQ=-K,K

                      ! Get the SEE index
                      iR = Atom(ia)%irho(iterml)%Jrho(iL1,iL)% &
                                    kq(iQ,K)

                      ! If flagged as small, skip
                      if (iR.le.0) cycle

                      ! Racah algebra
                      JSym = fun6j(1d0,1d0,rK,rJl1,rJl,rJu,Flgsg)

                    end do ! Q
                  end do ! K
                end do ! iL1
              end do ! iL
            end do ! iU

            !
            ! emiss2orbNB
            !

            ! Final term
            itermf = iterml
            rLf = Atom(ia)%rLval(itermf)

            ! For all the possible lower terms
            do i=1,Atom(ia)%nMulti-1

              ! If there is no transition or this term is larger
              ! than the upper term of the output transition, skip
              if(i.ge.itermu.or.Atom(ia)%irad(i,itermu).eq.0) cycle

              ! Store the input lower term index
              iterml = i

              ! Get index of input transition
              itran = Atom(ia)%irad(iterml,itermu)

              ! Angular momentum input lower level
              rLl = Atom(ia)%rLval(iterml)

              ! For each Jf
              do mF=1,Atom(ia)%nJ(itermf)

                ! Get Jf
                rJf = Atom(ia)%rJval(mF,itermf)

                ! For each Ju
                do iU=1,Atom(ia)%nJ(itermu)

                  ! Get Ju
                  rJu = Atom(ia)%rJval(iU,itermu)

                  JSym = fun6j(rJu,rJf,1d0,rLf,rLu,S,Flgsg)

                  ! Check size
                  if (abs(JSym).lt.TINYJS) cycle

                  ! For each Ju'
                  do iU1=1,Atom(ia)%nJ(itermu)

                    ! Get Ju'
                    rJu1 = Atom(ia)%rJval(iU1,itermu)

                    JSym = fun6j(rJu1,rJf,1d0,rLf,rLu,S,Flgsg)

                    ! Check size
                    if (abs(JSym).lt.TINYJS) cycle

                    ! For each Jl
                    do iL=1,Atom(ia)%nJ(iterml)

                      ! Get Jl
                      rJl = Atom(ia)%rJval(iL,iterml)

                      JSym = fun6j(rJu,rJl,1d0,rLl,rLu,S,Flgsg)

                      ! Check size
                      if (abs(JSym).lt.TINYJS) cycle

                      ! For each Jl'
                      do iL1=1,Atom(ia)%nJ(iterml)

                        ! Get Jl1
                        rJl1 = Atom(ia)%rJval(iL1,iterml)

                        JSym = fun6j(rJu1,rJl1,1d0,rLl,rLu,S,Flgsg)

                        ! Check size
                        if (abs(JSym).lt.TINYJS) cycle

       !
       ! Reset indexing
       !

       ! For each K
       do K=0,Krad

         ! Get real value
         rK = dble(K)

         ! For each Q
         do iQ=-K,K

           Q = dble(iQ)

           ! For each K'
           do K1=0,2

             ! Get real value
             rK1 = dble(K1)

             ! Racah algebra
             JSym = fun6j(rK1,rJu,rJu1,rJf,1d0,1d0,Flgsg)

             ! If forbidden (6j-sym=0), skip
             if(abs(JSym).lt.TINYJS) cycle

             ! For each Q'
             do iQ1=-K1,K1

               Q1 = dble(iQ1)

               iQl = iQ + iQ1
               Ql = dble(iQl)

               ! Determine the limits in K
               Kmin = max(abs(iQl),nint(abs(rK-rK1)), &
                          nint(abs(rJl-rJl1)))
               Kmax = min(nint(rJl+rJl1),nint(rK+rK1),Kcut)

               ! For each Kl
               do Kl=Kmin,Kmax

                 ! Check sum K + Kl
                 if((Kl+K).gt.Kcut) cycle

                 ! Get the real number
                 rKl = dble(Kl)

                 ! Get the SEE index
                 iR = Atom(ia)%irho(iterml)% &
                               Jrho(iL1,iL)%kq(iQl,Kl)

                 ! If flagged as small, skip
                 if (iR.le.0) cycle

                 ! Racah algebra
                 JSym = fun3j(rK,rK1,rKl,Q,Q1,-Ql,Flgsg)
                 JSym = fun9j(rK,rK1,rKl,1d0,rJu1,rJl1, &
                              1d0,rJu,rJl,Flgsg)

               end do ! Kl
             end do !Q'
           end do !K'
         end do ! Q
       end do ! K
                      end do ! iL1
                    end do ! iL
                  end do ! iU1
                end do ! iU
              end do ! mF
            end do ! Terms
          end do ! b-b transitions

        end if ! No magnetic field

      end do ! Each atom


      ! For each LTE line
      do ia=1,nLTEl

        ! For each Mu
        do iMu=1,LTElines(ia)%nMu

          ! Mu
          rMu = -LTElines(ia)%Ju + dble(iMu-1)

          ! For each Ml
          do iMl=1,LTElines(ia)%nMl

            ! Ml
            rMl = -LTElines(ia)%Jl + dble(iMl-1)

            ! q value
            q = rMl - rMu
            iq = nint(q)

            ! Selection rule
            if (abs(q).gt.1) cycle

            ! 3J
            Jsym = fun3j(LTElines(ia)%Ju,LTElines(ia)%Jl,1d0, &
                         -rMu,rMl,-q,Flgsg)

          end do ! Ml
        end do ! Mu

      end do ! NLTE lines

#ifdef _OPENMP
      ! Make memoization read_only
      Flgsg%can_write = .False.
#endif

      return

      end subroutine initmemoJ

!#####################################################################
!#####################################################################
!#####################################################################

      end module initmemoj_mod
