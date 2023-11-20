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
!     11/14/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/14/2023:    V3.0.2 - Added fieldB_alt (TdPA)
!
!     10/31/2023:    V3.0.1 - Added get_scattering,
!                             get_scattering_los, and
!                             scattering_manage subroutines (TdPA)
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
!  fieldB_alt:
!    Rotates JKQ in the vertical frame into the magnetic frame and the
!    opposite, but expects a different index ordering
!
!  rhoB:
!    Rotates rhoKQ in the vertical frame into the magnetic frame and
!    the opposite
!
!  get_scattering:
!    Compute all possible (unique) scattering angles for a given
!    quadrature, sorted by their angle value.
!
!  get_scattering_los:
!    Compute all possible (unique) scattering angles for a given
!    quadrature as input directions and a given output direction,
!    sorted by their angle value.
!
!  scattering_manage:
!    Decide which scattering angles have to be avoided for a given
!    quadrature (output) direction.
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
      use parameters_mod , only : cZero , cOne , cImag , pi , TINYA
      use qsort_mod
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

      !> Rotates radiation field tensors for an alternative index
      !! order.\n
      !!    JRad(dcmplx(:,:,:)): Radiation field tensor to rotate\n
      !!            nn(integer): Secondary dimensions in JRad besides
      !!                         K and Q\n
      !!     Flgsg(Fctsg_class): Structure with factorials and
      !!                         signs\n
      !!         thetaB(dfloat): Polar angle to rotate\n
      !!           phiB(dfloat): Azimuth to rotate\n
      !!           dir(integer): Direction of rotation (rotating
      !!                         forth or back)
      subroutine fieldB_alt(JRad,nn,Flgsg,thetaB,phiB,dir)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      integer, intent(in):: nn,dir
      double precision, intent(in):: thetaB, phiB
      complex(kind=8), dimension(-2:2,0:2,nn), intent(inout):: JRad

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

          Jaux(0,K) = sum(D(-K:K,0,K)*JRad(-K:K,K,ii))

          do iQ=1,K

            Jaux(iQ,K) = sum(D(-K:K,iQ,K)*JRad(-K:K,K,ii))
            Jaux(-iQ,K) = Flgsg%sg(iQ)*conjg(Jaux(iQ,K))

          end do
        end do

        ! Notice that K=0 is not changed
        JRad(-2:2,1:2,ii) = Jaux

      end do

      end subroutine fieldB_alt

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

      !> Generate an array of unique scattering angles for a given
      !! quadrature and index it.\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              los(logical): Indicates if we are normalizing
      !!                            LOS directions
      subroutine get_scattering(Geom,los)

      ! I/O
      type(Geometry_class):: Geom
      logical, intent(in):: los

      ! Local

      type(box_class), target:: box
      type(box_class), pointer:: boxr

      logical:: nofound

      integer:: ii,jj,kk,ll,mm,nn

      double precision:: Co,Coi,So,Soi,Ctheta,theta

      ! Initialize
      nullify(boxr,box%next)
      box%val = -1

      ! Clean
      if (allocated(Geom%i_scatt)) deallocate(Geom%i_scatt)
      if (allocated(Geom%j_scatt)) deallocate(Geom%j_scatt)
      if (allocated(Geom%k_scatt)) deallocate(Geom%k_scatt)
      if (allocated(Geom%skip_jsc)) deallocate(Geom%skip_jsc)
      if (allocated(Geom%skip_ksc)) deallocate(Geom%skip_ksc)
      if (allocated(Geom%V_CScatt)) deallocate(Geom%V_CScatt)
      if (allocated(Geom%V_SScatt)) deallocate(Geom%V_SScatt)
      Geom%nScatt = 0

      ! If LOS, you should not be computing this here
      if (los) return

      ! For each quadrature angle (output)
      do ii=1,Geom%nTh
        Co = Geom%V_mu(ii)
        So = sin(Geom%V_theta(ii))
        do jj=1,Geom%nPh
          ! For each quadrature angle (input)
          do kk=1,Geom%nTh
            Coi = Co*Geom%V_mu(kk)
            Soi = So*sin(Geom%V_theta(kk))
            do ll=1,Geom%nPh2

              ! Cosine scattering angle
              Ctheta = Coi + Soi*cos(Geom%V_phi(jj) - &
                                     Geom%V_phi(ll))

              ! Get angle
              if (Ctheta.ge.1.0) then
                theta = 0d0
              else if (Ctheta.le.-1.0) then
                theta = PI
              else
                theta = acos(Ctheta)
              end if

              ! If not started
              if (box%val.lt.0d0) then

                ! Add to first
                box%val = theta
                Geom%nScatt = 1

              ! Started already
              else

                ! Start finder
                nofound = .True.
                boxr => box

                ! Check till done
                do while (.True.)

                  ! Check current
                  if (abs(theta-boxr%val).lt.TINYA) then
                    nofound = .False.
                    exit
                  end if

                  ! Check next
                  if (.not.associated(boxr%next)) exit

                  ! Point to next
                  boxr => boxr%next

                end do

                ! If not found, add
                if (nofound) then

                  ! Add
                  Geom%nScatt = Geom%nScatt + 1

                  ! Allocate next
                  allocate(boxr%next)
                  boxr => boxr%next
                  nullify(boxr%next)
                  boxr%val = theta

                end if ! Add new scattering angle
              end if ! Need to check previous angles

            end do
          end do
        end do
      end do

      ! Allocate scattering angles
      allocate(Geom%V_CScatt(Geom%nScatt))
      allocate(Geom%V_SScatt(Geom%nScatt))

      ! Initialize runner
      boxr => box
      ii = 0

      ! Copy angles
      do while (.True.)

        ! Store
        ii = ii + 1
        Geom%V_CScatt(ii) = boxr%val

        ! If no more
        if (.not.associated(boxr%next)) exit

        ! Shift
        boxr => boxr%next

      end do

      ! Clean boxes
      do while (associated(box%next))

        ! Initialize
        boxr => box%next

        ! If more
        if (associated(boxr%next)) then

          ! Point to second to last
          do while (.True.)

            ! Second-to-last
            if (.not.associated(boxr%next%next)) exit

            ! Shift one
            boxr => boxr%next

          end do

          ! Remove tail
          deallocate(boxr%next)
          nullify(boxr%next)

        ! No more
        else

          ! Deallocate
          deallocate(boxr)
          nullify(box%next)
          nullify(boxr)

        end if

      end do ! Cleaning box

      ! Order scattering angles
      call QsortC(Geom%V_CScatt)

      !
      ! Index the directions
      !

      ! Allocate indexing
      allocate(Geom%i_scatt(Geom%nPh2,Geom%nTh, &
                            Geom%nPh2*Geom%nTh))

      ! Initialize running output direction
      mm = 0

      ! For each quadrature angle (output)
      do ii=1,Geom%nTh
        Co = Geom%V_mu(ii)
        So = sin(Geom%V_theta(ii))
        do jj=1,Geom%nPh
          mm = mm + 1
          ! For each quadrature angle (input)
          do kk=1,Geom%nTh
            Coi = Co*Geom%V_mu(kk)
            Soi = So*sin(Geom%V_theta(kk))
            do ll=1,Geom%nPh2

              ! Cosine scattering angle
              Ctheta = Coi + Soi*cos(Geom%V_phi(jj) - &
                                     Geom%V_phi(ll))

              ! Get angle
              if (Ctheta.ge.1.0) then
                theta = 0d0
              else if (Ctheta.le.-1.0) then
                theta = PI
              else
                theta = acos(Ctheta)
              end if

              ! Look for position
              do nn=1,Geom%nScatt
                if (abs(theta-Geom%V_CScatt(nn)).le.TINYA) then
                  Geom%i_scatt(ll,kk,mm) = nn
                  exit
                end if
              end do
            end do
          end do
        end do
      end do

      ! Compute cosines and sines
      Geom%V_SScatt = sin(Geom%V_CScatt)
      Geom%V_CScatt = cos(Geom%V_CScatt)

      ! Allocate
      allocate(Geom%skip_jsc(Geom%nScatt))
      allocate(Geom%j_scatt(Geom%nScatt))

      ! Initialize
      Geom%skip_jsc = .False.

      ! Order
      do ii=1,Geom%nScatt

        ! Trivial indexing
        Geom%j_scatt(ii) = ii

      end do

      ! No skip
      Geom%nskip = 0

      ! And done
      return

      end subroutine get_scattering

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate an array of unique scattering angles for a given
      !! quadrature and index it, for a given LOS.\n
      !!      Geom(Geometry_class): Structure with geometry data
      !!              ith(integer): Index of current polar LOS\n
      !!              iph(integer): Index of current azimuth LOS
      subroutine get_scattering_los(Geom,ith,iph)

      ! I/O
      type(Geometry_class):: Geom
      integer, intent(in):: ith,iph

      ! Local

      type(box_class), target:: box
      type(box_class), pointer:: boxr

      logical:: nofound

      integer:: ii,kk,ll,mm,nn

      double precision:: Co,Coi,So,Soi,Ctheta,theta

      ! Initialize
      nullify(boxr,box%next)
      box%val = -1

      ! Clean
      if (allocated(Geom%i_scatt)) deallocate(Geom%i_scatt)
      if (allocated(Geom%j_scatt)) deallocate(Geom%j_scatt)
      if (allocated(Geom%k_scatt)) deallocate(Geom%k_scatt)
      if (allocated(Geom%skip_jsc)) deallocate(Geom%skip_jsc)
      if (allocated(Geom%skip_ksc)) deallocate(Geom%skip_ksc)
      if (allocated(Geom%V_CScatt)) deallocate(Geom%V_CScatt)
      if (allocated(Geom%V_SScatt)) deallocate(Geom%V_SScatt)
      Geom%nScatt = 0

      ! For each quadrature angle (output)
      Co = Geom%L_mu(ith)
      So = sin(Geom%L_theta(ith))
      ! For each quadrature angle (input)
      do kk=1,Geom%nTh
        Coi = Co*Geom%V_mu(kk)
        Soi = So*sin(Geom%V_theta(kk))
        do ll=1,Geom%nPh2

          ! Cosine scattering angle
          Ctheta = Coi + Soi*cos(Geom%L_phi(iph) - &
                                 Geom%V_phi(ll))

          ! Get angle
          if (Ctheta.ge.1.0) then
            theta = 0d0
          else if (Ctheta.le.-1.0) then
            theta = PI
          else
            theta = acos(Ctheta)
          end if

          ! If not started
          if (box%val.lt.0d0) then

            ! Add to first
            box%val = theta
            Geom%nScatt = 1

          ! Started already
          else

            ! Start finder
            nofound = .True.
            boxr => box

            ! Check till done
            do while (.True.)

              ! Check current
              if (abs(theta-boxr%val).lt.TINYA) then
                nofound = .False.
                exit
              end if

              ! Check next
              if (.not.associated(boxr%next)) exit

              ! Point to next
              boxr => boxr%next

            end do

            ! If not found, add
            if (nofound) then

              ! Add
              Geom%nScatt = Geom%nScatt + 1

              ! Allocate next
              allocate(boxr%next)
              boxr => boxr%next
              nullify(boxr%next)
              boxr%val = theta

            end if ! Add new scattering angle
          end if ! Need to check previous angles

        end do
      end do

      ! Allocate scattering angles
      allocate(Geom%V_CScatt(Geom%nScatt))
      allocate(Geom%V_SScatt(Geom%nScatt))

      ! Initialize runner
      boxr => box
      ii = 0

      ! Copy angles
      do while (.True.)

        ! Store
        ii = ii + 1
        Geom%V_CScatt(ii) = boxr%val

        ! If no more
        if (.not.associated(boxr%next)) exit

        ! Shift
        boxr => boxr%next

      end do

      ! Clean boxes
      do while (associated(box%next))

        ! Initialize
        boxr => box%next

        ! If more
        if (associated(boxr%next)) then

          ! Point to second to last
          do while (.True.)

            ! Second-to-last
            if (.not.associated(boxr%next%next)) exit

            ! Shift one
            boxr => boxr%next

          end do

          ! Remove tail
          deallocate(boxr%next)
          nullify(boxr%next)

        ! No more
        else

          ! Deallocate
          deallocate(boxr)
          nullify(box%next)

        end if

      end do ! Cleaning box

      ! Order scattering angles
      call QsortC(Geom%V_CScatt)

      !
      ! Index the directions
      !

      ! Allocate indexing
      allocate(Geom%i_scatt(Geom%nPh2,Geom%nTh,1))

      ! Initialize running output direction
      mm = 0

      ! For each quadrature angle (output)
      Co = Geom%L_mu(ith)
      So = sin(Geom%L_theta(ith))
      ! For each quadrature angle (input)
      do kk=1,Geom%nTh
        Coi = Co*Geom%V_mu(kk)
        Soi = So*sin(Geom%V_theta(kk))
        do ll=1,Geom%nPh2

          ! Cosine scattering angle
          Ctheta = Coi + Soi*cos(Geom%L_phi(iph) - &
                                 Geom%V_phi(ll))

          ! Get angle
          if (Ctheta.ge.1.0) then
            theta = 0d0
          else if (Ctheta.le.-1.0) then
            theta = PI
          else
            theta = acos(Ctheta)
          end if

          ! Look for position
          do nn=1,Geom%nScatt
            if (abs(theta-Geom%V_CScatt(nn)).le.TINYA) then
              Geom%i_scatt(ll,kk,1) = nn
              exit
            end if
          end do
        end do
      end do

      ! Compute cosines and sines
      Geom%V_SScatt = sin(Geom%V_CScatt)
      Geom%V_CScatt = cos(Geom%V_CScatt)

      ! Allocate
      allocate(Geom%skip_jsc(Geom%nScatt))
      allocate(Geom%j_scatt(Geom%nScatt))

      ! Initialize
      Geom%skip_jsc = .False.

      ! Order
      do ii=1,Geom%nScatt

        ! Trivial indexing
        Geom%j_scatt(ii) = ii

      end do

      ! No skip
      Geom%nskip = 0

      ! And done
      return

      end subroutine get_scattering_los

!#####################################################################
!#####################################################################
!#####################################################################

      !> For a given output direction, manage what scattering angles
      !! are to be avoided.\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              ith(integer): Polar index of quadrature\n
      !!              iph(integer): Azimuthal index of quadrature
      subroutine scattering_manage(Geom,ith,iph)

      ! I/O
      type(Geometry_class):: Geom
      integer, intent(in):: ith,iph

      ! Local
      integer:: jdir,ith1,iph1,ish1


      ! Allocate
      if (allocated(Geom%k_scatt)) deallocate(Geom%k_scatt)
      if (allocated(Geom%skip_ksc)) deallocate(Geom%skip_ksc)
      allocate(Geom%skip_ksc(Geom%nScatt))
      allocate(Geom%k_scatt(Geom%nScatt))

      ! Initialize
      Geom%skip_ksc = .True.
      Geom%nskip = Geom%nScatt

      ! Output direction
      jdir = Geom%i_geom(iph,ith)

      ! Check scattering angles for this output direction
      do ith1=1,Geom%nTh
        do iph1=1,Geom%nPh2

          ! Index
          ish1 = Geom%i_scatt(iph1,ith1,jdir)

          ! If skipping
          if (Geom%skip_ksc(ish1)) then

            ! Flag no skip
            Geom%skip_ksc(ish1) = .False.
            Geom%nskip = Geom%nskip - 1

          end if

        end do
      end do

      ! Reindex scattering angles
      jdir = 0

      ! Initialize
      Geom%k_scatt = -1

      ! Run over scattering angles
      do ish1=1,Geom%nScatt

        ! Skip?
        if (Geom%skip_ksc(ish1)) cycle

        ! Get index
        jdir = jdir + 1
        Geom%k_scatt(ish1) = jdir

      end do

      end subroutine scattering_manage

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
