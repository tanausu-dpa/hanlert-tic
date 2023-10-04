      !> Maximum relative change computation
      module mrc_mod
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
!     09/29/2023 V3.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.3 - Updated to term-wise K cut limits (TdPA)
!
!     10/26/2022:    V3.0.2 - Changed the indexing of atomic levels
!                             in Atom (TdPA)
!
!     10/25/2022:    V3.0.1 - Implemented restriction of the height
!                             axis (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     02/20/2019:    V1.2.0 - Now uses specific TINY variables (TdPA)
!
!     11/01/2017:    V1.1.1 - MRC for intensity gives directly the
!                             level index (TdPA)
!
!     09/22/2017:    V1.1.0 - Possibility to limit K (TdPA)
!
!     09/08/2017:    V1.0.1 - Ignore rho00 below double precision
!                             and avoid dividing by 0 (TdPA)
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
!  MRC_sb
!    Calculate the Maximum Relative Change for rhoKQ
!
!  MRCI_sb
!    Calculate the Maximum Relative Change for rho00
!
!  MRCJ_sb
!    Calculate the Maximum Relative Change for J00C
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

      !> Calculates Maximum Relative Change for rhoKQ\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atom0(Atom_class): Structure with the atomic data for
      !!                       previous iteration\n
      !!       MRC(MRC_class): Structure with the MRC data
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
                    if(K.eq.0)then

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

                  end do ! K
                end do ! Q
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

      !> Calculates Maximum Relative Change for rho00\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!    Atom0(Atom_class): Structure with the atomic data for
      !!                       previous iteration\n
      !!       MRC(MRC_class): Structure with the MRC data
      subroutine MRCI_sb(Atom,Atom0,MRC)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Rhoc_class), dimension(:), intent(in):: Atom0
      type(MRC_class), intent(out):: MRC

      ! Local

      integer:: iz,ia,i,it,iJ,iR

      double precision:: MRHO,MRHO_old,MRCL,rJ


      !
      ! Initialize quantities
      !
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
              ! stored one, update data
              if(MRCL.gt.MRC%values(2,1))then

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

      !> Calculates Maximum Relative Change for J00 with frequency
      !! dependence\n
      !!     JC(dfloat(:,:)): Mean intensity\n
      !!  JCold(dfloat(:,:)): Mean intensity in previous iteration\n
      !!      MRC(MRC_class): Structure with the MRC data
      subroutine MRCJ_sb(JC,JCold,MRC)

      ! I/O

      double precision, dimension(nfreq,Rz0:Rz1):: JC, JCold
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
          ! stored one, update data
          if(MRCL.gt.MRC%values(2,1))then

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

      end module mrc_mod
