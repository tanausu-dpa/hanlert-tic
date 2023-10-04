      !> Geometrical tensors in magnetic field reference frame
      module btens_mod
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
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     04/20/2017:    V1.0.0 - First working version (TdPA)
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
!     This subroutine rotates the TKQ in the vertical reference frame
!   (TS) into the magnetic field reference frame (TB)
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

      !> Rotates geometrical tensors from the vertical to the
      !! magnetic field reference frame\n
      !!    TS(dcmplx(:,:,:)): Geometrical tensors in the vertical
      !!                       reference frame\n
      !!    TB(dcmplx(:,:,:)): Geometrical tensors in the magnetic
      !!                       field reference frame\n
      !!   Flgsg(Fctsg_class): Structure with factorials and
      !!                       signs\n
      !!       thetaB(dfloat): Polar angle to rotate\n
      !!         phiB(dfloat): Azimuth to rotate
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

      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      cexpPh(2) = cexpPh(1)*cexpPh(1)
      cexpPh(-2) = conjg(cexpPh(2))

      ! Initialize TB(i,Q,K)
      TB = cZero

      ! Initialize rotation matrix D(Q,Q1,K) for K=1,2
      do K=1,2

        rK = dble(K)

        do iQ=-K,K

          Q = dble(iQ)

          do iQ1=-K,K

            Q1 = dble(iQ1)
            D(iQ1,iQ,K) = cexpPh(-iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

          end do ! Qp

        end do ! Q

      end do ! K


      !
      ! Rotate
      !

      ! K=0
      TB(0,0,0) = TS(0,0,0)

      ! K=1,2
      do K=1,2

        do i=0,3
            TB(i,0,K) = sum(D(-K:K,0,K)*TS(i,-K:K,K))
        end do ! Stokes components

        do iQ=1,K
          do i=0,3

            TB(i,iQ,K) = sum(D(-K:K,iQ,K)*TS(i,-K:K,K))
            TB(i,-iQ,K) = Flgsg%sg(iQ)*conjg(TB(i,iQ,K))

          end do ! Stokes components
        end do ! Q

      end do ! K

      end subroutine Btens

!#####################################################################
!#####################################################################
!#####################################################################

      end module btens_mod
