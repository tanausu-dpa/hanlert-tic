      !> Maximum relative change computation
      module mrc_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     20/04/2017
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  MRC_sb
!    Calculate the Maximum Relative Change for rhoKQ
!
!  MRCI_sb
!    Calculate the Maximum Relative Change for rho00
!
!  MRCJ_sb
!    Calculate the Maximum Relative Change for the frequency dependent
!  mean intensity
!
!  MRCJKQ_sb
!    Calculate the Maximum Relative Change for the frequency dependent
!  radiation field tensors
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use parameters_mod , only : TINYMRC0 , TINYMRCR, TINYMRCJ
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the Maximum Relative Change for rhoKQ\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atom0(Rhoc_class(:)): Structure to store the density matrix
      !!                        of the previous iteration\n
      !!        MRC(MRC_class): Structure with the MRC data
      subroutine MRC_sb(Atom,Atom0,MRC)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Rhoc_class), dimension(:), intent(in):: Atom0
      type(MRC_class), intent(out):: MRC

      ! Local

      integer:: iz,ia,it,iJ,iJ1,K,iQ,iR,iI

      double precision:: MRCL
      double precision:: rho01,rhoo01,rhoo02,rho02,rho0,rhoo0,rJ,rJ1

      complex(kind=8):: MRHO,MRHO_old


      ! Initialize quantities
      MRC%values(2,1) = -1d99
      MRC%values(2,2) = -1d99
      MRC%indexes = 0
      MRC%indexes(2,:) = 1


      !
      ! Calculate MRC
      !

      ! For each atom
      do ia=1,nA

        ! For each height
        do iz=Rz0,Rz1

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J
            do iJ=1,Atom(ia)%nJ(it)

              ! Get angular momentum
              rJ = Atom(ia)%rJval(iJ,it)

              ! For each level J'
              do iJ1=1,Atom(ia)%nJ(it)

                ! Get angular momentum
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! If different level
                if(iJ.ne.iJ1)then

                  ! Get the indexes of the two levels
                  iR = Atom(ia)%irho(it)% &
                                Jrho(iJ,iJ)%kq(0,0)
                  iI = Atom(ia)%irho(it)% &
                                Jrho(iJ1,iJ1)%kq(0,0)

                  ! Get the rho00 of both levels
                  rho01 = dble(Atom(ia)%crho(iR,iz))
                  rho02 = dble(Atom(ia)%crho(iI,iz))
                  rhoo01 = dble(Atom0(ia)%crho(iR,iz))
                  rhoo02 = dble(Atom0(ia)%crho(iI,iz))

                  ! If both rho00 are small, do not bother
                  if (rho01.lt.TINYMRC0.and.rho02.lt.TINYMRC0.and. &
                      rhoo01.lt.TINYMRC0.and.rhoo02.lt.TINYMRC0) cycle

                  ! The rho00 are taken as a gemetric average
                  rho0 = sqrt(rho01*rho02)
                  rhoo0 = sqrt(rhoo01*rhoo02)

                  ! Get inverses if non-zero
                  if (rho0.gt.0d0) rho0 = 1d0/rho0
                  if (rhoo0.gt.0d0) rhoo0 = 1d0/rhoo0

                ! If same levels
                else

                  ! Get the index of the level
                  iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(0,0)

                  ! Get the rho00
                  rho0 = dble(Atom(ia)%crho(iR,iz))
                  rhoo0 = dble(Atom0(ia)%crho(iR,iz))

                  ! If both rho00 are small, don't bother
                  if (rho0.lt.TINYMRC0.and.rhoo0.lt.TINYMRC0) cycle

                  ! Get inverses if non-zero
                  if (rho0.gt.0d0) rho0 = 1d0/rho0
                  if (rhoo0.gt.0d0) rhoo0 = 1d0/rhoo0

                endif ! Different levels

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(it))

                  ! For each Q
                  do iQ=-K,K

                    ! Get the index of the rhoKQ
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! If very small (relative to rho00), ignore
                    if (abs(Atom(ia)%crho(iR,iz))*rho0.le. &
                        TINYMRCR) cycle

                    ! Simplify the rhoKQ involved in easier to handle
                    ! variables
                    MRHO = Atom(ia)%crho(iR,iz)
                    MRHO_old = Atom0(ia)%crho(iR,iz)

                    ! Determine the MRC avoiding dividing by 0
                    if(abs(MRHO_old*rhoo0).GT.TINYMRCR)then
                      MRCL = abs(MRHO - MRHO_old)/abs(MRHO_old)
                    else if(abs(MRHO*rho0).gt.TINYMRCR)then
                      MRCL = abs(MRHO - MRHO_old)/abs(MRHO)
                    else
                      MRCL = abs(MRHO - MRHO_old)
                    end if

                    ! If K is 0, population MRC
                    if (K.eq.0) then

                      ! If the new difference is larger than the last
                      ! stored one, update data
                      if(MRCL.gt.MRC%values(2,1))then
                        MRC%values(2,1) = MRCL
                        MRC%indexes(1,1) = ia
                        MRC%indexes(2,1) = iz
                        MRC%indexes(3,1) = it
                        MRC%indexes(4,1) = nint(2d0*rJ)
                        MRC%indexes(5,1) = nint(2d0*rJ1)
                      endif

                    ! If K is not 0, polarization MRC
                    else

                      ! If the new difference is larger than the last
                      ! stored one, update data
                      if(MRCL.gt.MRC%values(2,2))then

                        MRC%values(2,2) = MRCL
                        MRC%indexes(1,2) = ia
                        MRC%indexes(2,2) = iz
                        MRC%indexes(3,2) = it
                        MRC%indexes(4,2) = nint(2d0*rJ)
                        MRC%indexes(5,2) = nint(2d0*rJ1)
                        MRC%indexes(6,1) = K
                        MRC%indexes(6,2) = iQ

                      endif ! Larger RC
                    endif ! K = 0

                  end do ! Q
                end do ! K
              end do ! J1
            end do ! J
          end do ! term
        end do ! iz
      end do ! atom

      ! Avoid the -1d99 if you did not find anything not
      ! null
      if (MRC%indexes(3,1).eq.0) MRC%values(2,1) = 0d0
      if (MRC%indexes(3,2).eq.0) MRC%values(2,2) = 0d0

      end subroutine MRC_sb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the Maximum Relative Change for rho00\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atom0(Rhoc_class(:)): Structure to store the density matrix
      !!                        of the previous iteration\n
      !!        MRC(MRC_class): Structure with the MRC data
      subroutine MRCI_sb(Atom,Atom0,MRC)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Rhoc_class), dimension(:), intent(in):: Atom0
      type(MRC_class), intent(out):: MRC

      ! Local

      integer:: iz,ia,i,it,iJ,iR

      double precision:: MRHO,MRHO_old,MRCL,rJ


      ! Initialize quantities
      MRC%values(2,1) = -1d99
      MRC%indexes = 0
      MRC%indexes(2,1) = 1


      !
      ! Calculate MRC
      !

      ! For each atom
      do ia=1,nA

        ! For each height
        do iz=Rz0,Rz1

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J
            do iJ=1,Atom(ia)%nJ(it)

              ! Get angular momentum
              rJ = Atom(ia)%rJval(iJ,it)

              ! Get level index
              i = Atom(ia)%irho(it)%irho_ij(iJ)

              ! Get the index of the rho00
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! Simplify the rhoKQ involved in easier to handle
              ! variables
              MRHO = dble(Atom(ia)%crho(iR,iz))
              MRHO_old = dble(Atom0(ia)%crho(iR,iz))

              ! If both rho00 are small, do not bother
              if (MRHO.lt.TINYMRC0.and.MRHO_old.lt.TINYMRC0) cycle

              ! Determine the MRC avoiding dividing by 0
              if(MRHO_old.GT.TINYMRCR)then
                MRCL = abs(MRHO - MRHO_old)/MRHO_old
              else if(MRHO.gt.TINYMRCR)then
                MRCL = abs(MRHO - MRHO_old)/MRHO
              else
                MRCL = abs(MRHO - MRHO_old)
              end if

              ! If the new difference is larger than the last
              ! stored one
              if(MRCL.gt.MRC%values(2,1))then

                ! Update data
                MRC%values(2,1) = MRCL
                MRC%indexes(1,1) = ia
                MRC%indexes(2,1) = iz
                MRC%indexes(3,1) = i

              end if ! Larger RC

            end do !J
          end do !term
        end do !iz
      end do !atom

      ! Avoid the -1d99 if you did not find anything not
      ! null
      if (MRC%indexes(3,1).eq.0) MRC%values(2,1) = 0d0

      end subroutine MRCI_sb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the Maximum Relative Change for the frequency
      !! dependent mean intensity\n
      !!     JC(double(:,:)): Frequency dependent mean intensity\n
      !!  JCold(double(:,:)): Frequency dependent mean intensity in
      !!                      the previous iteration\n
      !!      MRC(MRC_class): Structure with the MRC data
      subroutine MRCJ_sb(JC,JCold,MRC)

      ! I/O

      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: JC
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: JCold
      type(MRC_class), intent(out):: MRC

      ! Local

      integer:: iz,ifreq

      double precision:: new,old,MRCL


      !
      ! Initialize quantities
      !
      MRC%values(2,1) = -1d99
      MRC%indexes = 0
      MRC%indexes(2,1) = 1


      !
      ! Calculate MRC
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! Simplify the J involved in easier to handle variables
          new = JC(ifreq,iz)
          old = JCold(ifreq,iz)

          ! Determine the MRC avoiding dividing by 0
          if(old.GT.TINYMRCJ)then
            MRCL = abs(new - old)/old
          else if(old.gt.TINYMRCJ)then
            MRCL = abs(new - old)/new
          else
            MRCL = abs(new - old)
          end if

          ! If the new difference is larger than the last
          ! stored one
          if(MRCL.gt.MRC%values(2,1))then

            ! Update data
            MRC%values(2,1) = MRCL
            MRC%indexes(1,1) = ifreq
            MRC%indexes(2,1) = iz

          end if ! Larger RC

        end do !frequency
      end do !iz

      end subroutine MRCJ_sb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the Maximum Relative Change for the frequency
      !! dependent radiation field tensors\n
      !!     JC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!  JCold(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                            frequency dependence in the
      !!                            previous iteration\n
      !!            MRC(MRC_class): Structure with the MRC data
      subroutine MRCJKQ_sb(JC,JCold,MRC)

      ! I/O

      complex(kind=8), dimension(-2:2,0:2,nfreq,Rz0:Rz1), &
                       intent(in):: JC
      complex(kind=8), dimension(0:2,0:2,nfreq,Rz0:Rz1), &
                        intent(in):: JCold
      type(MRC_class), intent(out):: MRC

      ! Local

      integer:: iz,ifreq,K,iQ

      double precision:: JC0,JCold0,MRCL

      complex(kind=8):: JCK,JColdK


      !
      ! Initialize quantities
      !
      MRC%values(2,1) = -1d99
      MRC%values(2,2) = -1d99
      MRC%indexes = 0
      MRC%indexes(2,:) = 1


      !
      ! Calculate MRC
      !

      ! For each height
      do iz=Rz0,Rz1

        ! For each frequency
        do ifreq=1,nfreq

          ! K = 0
          JC0 = dble(JC(0,0,ifreq,iz))
          JCold0 = dble(JCold(0,0,ifreq,iz))

          ! Skip if too small
          if (JC0.lt.TINYMRCJ.and.JCold0.lt.TINYMRCJ) cycle

          ! Get inverse
          JC0 = 1d0/JC0
          JCold0 = 1d0/JCold0

          ! For each K
          do K=0,Krad

            ! For each Q
            do iQ=0,K

              ! Skip if axial
              if (axial.and.iQ.gt.0) cycle

              ! Get values
              JCK = JC(iQ,K,ifreq,iz)
              JColdK = JCold(iQ,K,ifreq,iz)

              ! If very small, ignore
              if (abs(JCK*JC0).le.TINYMRCR.and. &
                  abs(JColdK*JCold0).le.TINYMRCR) cycle

              ! Determine the MRC avoiding dividing by 0
              if (abs(JColdK*JCold0).gt.TINYMRCR) then
                MRCL = abs(JCK - JColdK)/abs(JColdK)
              else if (abs(JCK*JC0).gt.TINYMRCR) then
                MRCL = abs(JCK - JColdK)/abs(JCK)
              else
                MRCL = abs(JCK - JColdK)
              end if

              ! Mean intensity
              if (K.eq.0) then

                ! The new difference is larger
                if (MRCL.gt.MRC%values(2,1)) then

                  ! Update data
                  MRC%values(2,1) = MRCL
                  MRC%indexes(1,1) = ifreq
                  MRC%indexes(2,1) = iz

                end if ! Larger RC

              ! Not mean intensity
              else

                ! The new difference is larger
                if (MRCL.gt.MRC%values(2,2)) then

                  ! Update data
                  MRC%values(2,2) = MRCL
                  MRC%indexes(1,2) = ifreq
                  MRC%indexes(2,2) = iz
                  MRC%indexes(3,2) = K
                  MRC%indexes(4,2) = iQ

                end if ! Larger RC
              end if ! Mean intensity
            end do ! Q
          end do ! K
        end do !frequency
      end do !iz

      end subroutine MRCJKQ_sb

!#####################################################################
!#####################################################################
!#####################################################################

      end module mrc_mod
