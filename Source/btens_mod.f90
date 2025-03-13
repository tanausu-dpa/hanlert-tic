      !> Geometrical tensors in magnetic field reference frame
      module btens_mod
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
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Revised headers (TdPA)
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
!  Btens
!    Rotate geometrical TKQ tensors from the vertical to the magnetic
!  reference frames
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use rdmat_mod
      use parameters_mod , only : cZero, cOne, cImag
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rotate geometrical TKQ tensors from the vertical to the
      !! magnetic reference frames\n
      !   TS(dcomplx(:,:,:)): Geometrical tensors in the vertical
      !!                      reference frame\n
      !!  TB(dcomplx(:,:,:)): Geometrical tensors in the magnetic
      !!                      field reference frame\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!      thetaB(double): Polar angle of the magnetic field in the
      !!                      vertical reference frame\n
      !!        phiB(double): Azimuthal angle of the magnetic field in
      !!                      the vertical reference frame
      subroutine Btens(TS,TB,Flgsg,thetaB,phiB)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: thetaB, phiB
      complex(kind=8), intent(in):: TS(0:3,-2:2,0:2)
      complex(kind=8), intent(out):: TB(0:3,-2:2,0:2)

      ! Local

      integer:: i,K,iQ,iQ1

      double precision:: rK,Q,Q1

      complex(kind=8), dimension(-2:2):: cexpPh
      complex(kind=8), dimension(-2:2,-2:2,0:2):: D


      !
      ! Initializations
      !

      ! Initialize exponentials in the rotation matrix
      cexpPh(0) = cOne

      ! chiB
      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      ! 2chiB
      cexpPh(2) = cexpPh(1)*cexpPh(1)
      cexpPh(-2) = conjg(cexpPh(2))

      ! Initialize TB(i,Q,K)
      TB = cZero

      ! Initialize rotation matrix D(Q,Q1,K) for K=1,2
      do K=1,2

        ! Double precision K
        rK = dble(K)

        ! For each Q (only non-negative values are used below)
        do iQ=0,K

          ! Double precision Q
          Q = dble(iQ)

          ! For each Q'
          do iQ1=-K,K

            ! Double precision Q'
            Q1 = dble(iQ1)

            ! Rotation matrix
            D(iQ1,iQ,K) = cexpPh(-iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

          end do ! Qp
        end do ! Q
      end do ! K


      !
      ! Rotate
      !

      ! K=0 does not rotate
      TB(0,0,0) = TS(0,0,0)

      ! K=1,2
      do K=1,2

        !
        ! Q=0 component
        !

        ! Stokes parameters
        do i=0,3

          ! Rotate
          TB(i,0,K) = sum(D(-K:K,0,K)*TS(i,-K:K,K))

        end do ! Stokes components

        ! Rest of Q
        do iQ=1,K

          ! Stokes parameters
          do i=0,3

            ! Rotate Q>0
            TB(i,iQ,K) = sum(D(-K:K,iQ,K)*TS(i,-K:K,K))

            ! Use dependence relations for Q<0
            TB(i,-iQ,K) = Flgsg%sg(iQ)*conjg(TB(i,iQ,K))

          end do ! Stokes parameters
        end do ! Q > 0
      end do ! K

      end subroutine Btens

!#####################################################################
!#####################################################################
!#####################################################################

      end module btens_mod
