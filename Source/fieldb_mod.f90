      !> Rotation of radiation field tensors and density matrix
      module fieldb_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Contributors:
!     John Dennis (NCAR)
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
!     06/05/2020:    V1.0.3 - Changed the index order in JRad in
!                             fieldB to accomodate the optimization
!                             in emiss2ord (JD)
!
!     01/23/2019:    V1.0.2 - Bugfix: The rotation of rhoKQ was wrong
!                             for quantum interference, it had the
!                             wrong relation for rhoK-Q (TdPA)
!
!     04/24/2018:    V1.0.1 - Added check of arcos argument in
!                             atom2lab (TdPA)
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
!  fieldB:
!    Rotates JKQ in the vertical frame into the magnetic frame and the
!    opposite
!
!  rhoB:
!    Rotates rhoKQ in the vertical frame into the magnetic frame and
!    the opposite
!
!  atom2lab:
!    calculate the Theta angle between two directions
!    (th1,ph1) and (th2,ph2)
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use rdmat_mod
      use parameters_mod , only : cZero , cOne , cImag , pi
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rotates radiation field tensors.\n
      !!    JRad(dcmplx(:,:,:)): Radiation field tensor to rotate\n
      !!            nn(integer): Secondary dimensions in JRad besides
      !!                         K and Q\n
      !!     Flgsg(Fctsg_class): Structure with factorials and
      !!                         signs\n
      !!         thetaB(dfloat): Polar angle to rotate\n
      !!           phiB(dfloat): Azimuth to rotate\n
      !!           dir(integer): Direction of rotation (rotating
      !!                         forth or back)
      subroutine fieldB(JRad,nn,Flgsg,thetaB,phiB,dir)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: nn,dir
      double precision, intent(in):: thetaB, phiB
      complex(kind=8), dimension(nn,-2:2,0:2), intent(inout):: JRad

      ! Local

      integer:: K,iQ,iQ1,ii

      double precision:: rK,Q,Q1

      complex(kind=8), dimension(-2:2):: cexpPh
      complex(kind=8), dimension(-2:2,1:2):: Jaux
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

      ! Initialize
      Jaux = cZero

      ! Initialize rotation matrix D(Q,Q1,K) for K=1,2

      ! If going from vertical to magnetic field
      if (dir.gt.0) then

        do K=1,2
          rK = dble(K)
          do iQ = -K,K
            Q = dble(iQ)
            do iQ1=-K,K

              Q1 = dble(iQ1)
              D(iQ1,iQ,K) = cexpPh(-iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do
          end do
        end do

      ! If going from magnetic field to vertical
      else

        do K=1,2
          rK = dble(K)
          do iQ = -K,K
            Q = dble(iQ)
            do iQ1=-K,K

              Q1 = dble(iQ1)
              D(iQ1,iQ,K) = cexpPh(-iQ)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do
          end do
        end do

      end if ! Direction of rotation


      !
      ! Rotate
      !

      ! For each input element
      do ii=1,nn

        do K=1,2

          Jaux(0,K) = sum(D(-K:K,0,K)*JRad(ii,-K:K,K))

          do iQ=1,K

            Jaux(iQ,K) = sum(D(-K:K,iQ,K)*JRad(ii,-K:K,K))
            Jaux(-iQ,K) = Flgsg%sg(iQ)*conjg(Jaux(iQ,K))

          end do
        end do

        ! Notice that K=0 is not changed
        JRad(ii,-2:2,1:2) = Jaux

      end do

      end subroutine fieldB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rotates density matrix tensors.\n
      !!       rho(dcmplx(:,:)): Density matrix tensor to rotate\n
      !!            nn(integer): Secondary dimensions in rho besides
      !!                         Q\n
      !!             K(integer): K value of the incoming rho\n
      !!     Flgsg(Fctsg_class): Structure with factorials and
      !!                         signs\n
      !!         thetaB(dfloat): Polar angle to rotate\n
      !!           phiB(dfloat): Azimuth to rotate\n
      !!           dir(integer): Direction of rotation (rotating
      !!                         forth or back)
      subroutine rhoB(rho,nn,K,Flgsg,thetaB,phiB,dir)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: nn,K,dir
      double precision, intent(in):: thetaB, phiB
      complex(kind=8), dimension(-K:K,nn), intent(inout):: rho

      ! Local

      integer:: iQ,iQ1,ii

      double precision:: rK,Q,Q1

      complex(kind=8), dimension(-K:K):: cexpPh
      complex(kind=8), dimension(-K:K):: rhoaux
      complex(kind=8), dimension(-K:K,-K:K):: D


      ! The multipole 0 does not rotate
      if(K.eq.0)return


      !
      ! Initializations
      !

      ! Initialize exponentials in the rotation matrix
      cexpPh(0) = cOne

      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      do iQ=2,K
        cexpPh(iQ) = cexpPh(1)*cexpPh(iQ-1)
        cexpPh(-iQ) = conjg(cexpPh(iQ))
      end do

      ! Initialize
      rhoaux = cZero

      ! Convert K to real
      rK = dble(K)

      ! If going from vertical to magnetic field
      if (dir.gt.0) then

        do iQ = -K,K
          Q = dble(iQ)
          do iQ1=-K,K

            Q1 = dble(iQ1)
            D(iQ1,iQ) = cexpPh(iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

          end do
        end do

      ! If going from magnetic field to vertical
      else

        do iQ = -K,K
          Q = dble(iQ)
          do iQ1=-K,K

            Q1 = dble(iQ1)
            D(iQ1,iQ) = cexpPh(iQ)*rdmat(rK,Q1,Q,Flgsg,thetaB)

          end do
        end do

      end if ! Direction of rotation


      !
      ! Rotate
      !

      ! For each input element
      do ii=1,nn

        rhoaux(0) = sum(D(-K:K,0)*rho(-K:K,ii))

        do iQ=1,K

          rhoaux(iQ) = sum(D(-K:K,iQ)*rho(-K:K,ii))
          rhoaux(-iQ) = sum(D(-K:K,-iQ)*rho(-K:K,ii))

        end do

        rho(-K:K,ii) = rhoaux

      end do

      end subroutine rhoB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes angle between two directions\n
      !!   th1(dfloat): Polar angle direction 1\n
      !!   ph1(dfloat): Azimuth direction 1\n
      !!   th2(dfloat): Polar angle direction 2\n
      !!   ph2(dfloat): Azimuth direction 2
      double precision function atom2lab(th1,ph1,th2,ph2)

      ! I/O
      double precision,intent(in):: th1,ph1,th2,ph2

      ! Local
      double precision:: CTheta

      CTheta = cos(th1)*cos(th2) + &
               sin(th1)*sin(th2)*cos(ph2-ph1)

      ! If over 1 by numerical noise
      if (CTheta.ge.1d0) then

        atom2lab = 0d0

      ! If below -1 by numerical noise
      elseif (CTheta.le.-1d0) then

        atom2lab = pi

      ! Normal case
      else

        atom2lab = acos(CTheta)

      endif

      end function atom2lab

!#####################################################################
!#####################################################################
!#####################################################################

      end module fieldb_mod
