      !> Rotation of radiation field tensors and density matrix
      module fieldb_mod
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
!     03/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     03/12/2024:    V4.0.0 - The determination of the scattering
!                             does not need to consider particular
!                             cases for the different propagation
!                             directions (TdPA)
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
!  fieldB
!    Rotate the radiation field tensors to and from the magnetic
!  field reference frame
!
!  fieldB_alt:
!    Rotate the radiation field tensors to and from the magnetic
!  field reference frame expecting a different indexing
!
!  rhoB:
!    Rotate the density matrix tensors for a given multipole K to and
!  from the magnetic field reference frame
!
!  get_scattering:
!    Generate an array of unique scattering angles for the quadrature
!
!  get_scattering_los:
!    Generate an array of unique scattering angles for the a given
!  line of sight
!
!  atom2lab:
!    Compute the angle between two directions
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

      !> Rotate the radiation field tensors to and from the magnetic
      !! field reference frame\n
      !!  JRad(dcomplx(:,:,:)): Radiation field tensor to rotate\n
      !!           nn(integer): Secondary dimensions in JRad besides
      !!                        K and Q (first dimension)\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!       thetaB(double): Polar angle to rotate\n
      !!         phiB(double): Azimuth to rotate\n
      !!         dir(integer): Direction of rotation (rotating
      !!                       forth [1] or back [-1] the magnetic
      !!                       field reference frame)
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
      complex(kind=8), dimension(-2:2,0:2,0:2):: D


      !
      ! Initializations
      !

      ! Initialize exponentials in the rotation matrix
      cexpPh(0) = cOne

      ! First power
      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      ! Second power
      cexpPh(2) = cexpPh(1)*cexpPh(1)
      cexpPh(-2) = conjg(cexpPh(2))

      ! Initialize auxiliar tensor variable
      Jaux = cZero

      ! Initialize rotation matrix D(Q,Q1,K) for K=1,2

      ! If going from vertical to magnetic field
      if (dir.gt.0) then

        ! For each K
        do K=1,2

          ! Get double precision value
          rK = dble(K)

          ! For each Q (only using >=0)
          do iQ=0,K

            ! Get double precision value
            Q = dble(iQ)

            ! For each Q'
            do iQ1=-K,K

              ! Get double precision value
              Q1 = dble(iQ1)

              ! Get rotation matrix
              D(iQ1,iQ,K) = cexpPh(-iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do ! Q'
          end do ! Q
        end do ! K

      ! If going from magnetic field to vertical
      else

        ! For each K
        do K=1,2

          ! Get double precision value
          rK = dble(K)

          ! For each Q (only using >=0)
          do iQ=0,K

            ! Get double precision value
            Q = dble(iQ)

            ! For each Q'
            do iQ1=-K,K

              ! Get double precision value
              Q1 = dble(iQ1)

              ! Get rotation matrix
              D(iQ1,iQ,K) = cexpPh(-iQ)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do ! Q'
          end do ! Q
        end do ! K

      end if ! Direction of rotation


      !
      ! Rotate
      !

      ! For each input element
      do ii=1,nn

        ! For each K multipole
        do K=1,2

          ! Get Q=0 component
          Jaux(0,K) = sum(D(-K:K,0,K)*JRad(ii,-K:K,K))

          ! For the positive Q values
          do iQ=1,K

            ! Get proper rotation
            Jaux(iQ,K) = sum(D(-K:K,iQ,K)*JRad(ii,-K:K,K))

            ! Use tensor relations
            Jaux(-iQ,K) = Flgsg%sg(iQ)*conjg(Jaux(iQ,K))

          end do ! Q
        end do ! K

        ! Notice that K=0 is not changed
        JRad(ii,-2:2,1:2) = Jaux

      end do ! Array elements

      end subroutine fieldB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rotate the radiation field tensors to and from the magnetic
      !! field reference frame expecting a different indexing\n
      !!  JRad(dcomplx(:,:,:)): Radiation field tensor to rotate\n
      !!           nn(integer): Secondary dimensions in JRad besides
      !!                        K and Q (third dimension)\n
      !!    Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                        J-symbols\n
      !!        thetaB(double): Polar angle to rotate\n
      !!          phiB(double): Azimuth to rotate\n
      !!          dir(integer): Direction of rotation (rotating
      !!                        forth [1] or back [-1] the magnetic
      !!                        field reference frame)
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
      complex(kind=8), dimension(-2:2,0:2,0:2):: D


      !
      ! Initializations
      !

      ! Initialize exponentials in the rotation matrix
      cexpPh(0) = cOne

      ! First power
      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      ! Second power
      cexpPh(2) = cexpPh(1)*cexpPh(1)
      cexpPh(-2) = conjg(cexpPh(2))

      ! Initialize
      Jaux = cZero

      ! Initialize rotation matrix D(Q,Q1,K) for K=1,2

      ! If going from vertical to magnetic field
      if (dir.gt.0) then

        ! For each K
        do K=1,2

          ! Get double precision value
          rK = dble(K)

          ! For each Q (only using >=0)
          do iQ=0,K

            ! Get double precision value
            Q = dble(iQ)

            ! For each Q'
            do iQ1=-K,K

              ! Get double precision value
              Q1 = dble(iQ1)

              ! Get rotation matrix
              D(iQ1,iQ,K) = cexpPh(-iQ1)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do ! Q'
          end do ! Q
        end do ! K

      ! If going from magnetic field to vertical
      else

        ! For each K
        do K=1,2

          ! Get double precision value
          rK = dble(K)

          ! For each Q (only using >=0)
          do iQ=0,K

            ! Get double precision value
            Q = dble(iQ)

            ! For each Q'
            do iQ1=-K,K

              ! Get double precision value
              Q1 = dble(iQ1)

              ! Get rotation matrix
              D(iQ1,iQ,K) = cexpPh(-iQ)*rdmat(rK,Q1,Q,Flgsg,thetaB)

            end do ! Q'
          end do ! Q
        end do ! K

      end if ! Direction of rotation


      !
      ! Rotate
      !

      ! For each input element
      do ii=1,nn

        ! For each K multipole
        do K=1,2

          ! Get Q=0 component
          Jaux(0,K) = sum(D(-K:K,0,K)*JRad(-K:K,K,ii))

          ! For the positive Q values
          do iQ=1,K

            ! Get proper rotation
            Jaux(iQ,K) = sum(D(-K:K,iQ,K)*JRad(-K:K,K,ii))

            ! Use tensor relations
            Jaux(-iQ,K) = Flgsg%sg(iQ)*conjg(Jaux(iQ,K))

          end do ! Q
        end do ! K

        ! Notice that K=0 is not changed
        JRad(-2:2,1:2,ii) = Jaux

      end do ! Array elements

      end subroutine fieldB_alt

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rotate the density matrix tensors for a given multipole K to
      !! and from the magnetic field reference frame\n
      !!    rho(dcomplx(:,:)): Density matrix tensor to rotate\n
      !!          nn(integer): Secondary dimensions in JRad besides
      !!                       K and Q (second dimension)\n
      !!   Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                       J-symbols\n
      !!       thetaB(double): Polar angle to rotate\n
      !!         phiB(double): Azimuth to rotate\n
      !!         dir(integer): Direction of rotation (rotating
      !!                       forth [1] or back [-1] the magnetic
      !!                       field reference frame)
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

      ! First power
      cexpPh(1) = exp(cImag*phiB)
      cexpPh(-1) = conjg(cexpPh(1))

      ! Second and beyond powers
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

        ! Get Q=0 component
        rhoaux(0) = sum(D(-K:K,0)*rho(-K:K,ii))

        ! For rest of Q values
        do iQ=1,K

          rhoaux(iQ) = sum(D(-K:K,iQ)*rho(-K:K,ii))
          rhoaux(-iQ) = sum(D(-K:K,-iQ)*rho(-K:K,ii))

        end do

        ! Save in actual variable
        rho(-K:K,ii) = rhoaux

      end do ! Input elements

      end subroutine rhoB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate an array of unique scattering angles for the
      !! quadrature\n
      !!  Geom(Geometry_class): Structure with geometric data
      subroutine get_scattering(Geom)

      ! I/O

      type(Geometry_class), intent(inout):: Geom

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
      if (allocated(Geom%i_scatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%i_scatt)
        deallocate(Geom%i_scatt)
      end if
      if (allocated(Geom%V_CScatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%V_CScatt)
        deallocate(Geom%V_CScatt)
      end if
      if (allocated(Geom%V_SScatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%V_SScatt)
        deallocate(Geom%V_SScatt)
      end if
      Geom%nScatt = 0

      ! For each quadrature polar angle (output)
      do ii=1,Geom%nTh

        ! Get trigonometry
        Co = Geom%V_mu(ii)
        So = sin(Geom%V_theta(ii))

        ! For each quadrature azimuth (output)
        do jj=1,Geom%nPh

          ! For each quadrature polar angle (input)
          do kk=1,Geom%nTh

            ! Get trigonometry
            Coi = Co*Geom%V_mu(kk)
            Soi = So*sin(Geom%V_theta(kk))

            ! For each quadrature azimuth (input)
            do ll=1,Geom%nPh2

              ! Cosine scattering angle
              Ctheta = Coi + Soi*cos(Geom%V_phi(jj) - &
                                     Geom%V_phi(ll))

              ! If cosine exact 1 or overflown
              if (Ctheta.ge.1.0) then

                ! Angle is zero
                theta = 0d0

              ! If cosine exact -1 or overflown
              else if (Ctheta.le.-1.0) then

                ! Angle is 180
                theta = PI

              ! Intermediate values
              else

                ! Arccos
                theta = acos(Ctheta)

              end if ! Control invalid cosines

              ! If list not started
              if (box%val.lt.0d0) then

                ! Add to first
                box%val = theta
                Geom%nScatt = 1

              ! List already started
              else

                ! Initialize value finder
                nofound = .True.
                boxr => box

                ! Check till done
                do while (.True.)

                  ! If the new angle is close enough to the one
                  ! in this box
                  if (abs(theta-boxr%val).lt.TINYA) then

                    ! Signal found and leave
                    nofound = .False.
                    exit

                  end if ! Close angle

                  ! If last one, exit loop
                  if (.not.associated(boxr%next)) exit

                  ! Point to next
                  boxr => boxr%next

                end do ! Till done with the list

                ! If not found, add
                if (nofound) then

                  ! Add one new angle
                  Geom%nScatt = Geom%nScatt + 1

                  ! Allocate next
                  allocate(boxr%next)
                  boxr => boxr%next
                  nullify(boxr%next)

                  ! Save value
                  boxr%val = theta

                end if ! Add new scattering angle
              end if ! Need to check previous angles

            end do ! Input azimuth
          end do ! Input polar angle
        end do ! Output azimuth
      end do ! Output polar angle

      ! Allocate scattering angles
      allocate(Geom%V_CScatt(Geom%nScatt))
      allocate(Geom%V_SScatt(Geom%nScatt))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_CScatt)
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_SScatt)

      ! Initialize runner
      boxr => box
      ii = 0

      ! Copy angles
      do while (.True.)

        ! Store
        ii = ii + 1
        Geom%V_CScatt(ii) = boxr%val

        ! If no more boxes, leave
        if (.not.associated(boxr%next)) exit

        ! Shift
        boxr => boxr%next

      end do ! Until everything is copied

      !
      ! Clean boxes
      !

      ! While there is data in boxes
      do while (associated(box%next))

        ! Initialize pointer
        boxr => box%next

        ! If more boxes to the right
        if (associated(boxr%next)) then

          ! Point to second to last
          do while (.True.)

            ! Second-to-last found, leave
            if (.not.associated(boxr%next%next)) exit

            ! Shift one
            boxr => boxr%next

          end do ! Until reached second to last

          ! Remove tail
          deallocate(boxr%next)
          nullify(boxr%next)

        ! No more to the right
        else

          ! Deallocate
          deallocate(boxr)
          nullify(box%next)
          nullify(boxr)

        end if ! More to the right of the box

      end do ! Until the box is clean

      ! Order scattering angles
      call QsortC(Geom%V_CScatt)

      !
      ! Index the directions
      !

      ! Allocate indexing for quadrature
      allocate(Geom%i_scatt(Geom%nPh2,Geom%nTh, &
                            Geom%nPh*Geom%nTh))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Geom%i_scatt)

      ! Initialize running output direction
      mm = 0

      ! For each quadrature polar angle (output)
      do ii=1,Geom%nTh

        ! Get trigonometry
        Co = Geom%V_mu(ii)
        So = sin(Geom%V_theta(ii))

        ! For each quadrature azimuth angle (output)
        do jj=1,Geom%nPh

          ! Advance output direction
          mm = mm + 1

          ! For each quadrature polar angle (input)
          do kk=1,Geom%nTh

            ! Get trigonometry
            Coi = Co*Geom%V_mu(kk)
            Soi = So*sin(Geom%V_theta(kk))

            ! For each quadrature azimuth angle (input)
            do ll=1,Geom%nPh2

              ! Cosine scattering angle
              Ctheta = Coi + Soi*cos(Geom%V_phi(jj) - &
                                     Geom%V_phi(ll))

              ! If cosine exact 1 or overflown
              if (Ctheta.ge.1.0) then

                ! Angle is zero
                theta = 0d0

              ! If cosine exact -1 or overflown
              else if (Ctheta.le.-1.0) then

                ! Angle is 180
                theta = PI

              ! Intermediate values
              else

                ! Arccos
                theta = acos(Ctheta)

              end if ! Control invalid cosines

              !
              ! Find angle position in vector
              !

              ! For every scattering angle
              do nn=1,Geom%nScatt

                ! If same angle
                if (abs(theta-Geom%V_CScatt(nn)).le.TINYA) then

                  ! Save index and leave
                  Geom%i_scatt(ll,kk,mm) = nn
                  exit

                end if ! Same angle

              end do ! Scattering angle
            end do ! Quadrature azimuth (input)
          end do ! Quadrature polar (input)
        end do ! Quadrature azimuth (output)
      end do ! Quadrature polar (output)

      ! Compute cosines and sines for scattering angles
      Geom%V_SScatt = sin(Geom%V_CScatt)
      Geom%V_CScatt = cos(Geom%V_CScatt)

      ! And done
      return

      end subroutine get_scattering

!#####################################################################
!#####################################################################
!#####################################################################

      !> Generate an array of unique scattering angles for the a given
      !! line of sight\n
      !!  Geom(Geometry_class): Structure with geometric data
      !!          ith(integer): Index of current polar LOS\n
      !!          iph(integer): Index of current azimuth LOS
      subroutine get_scattering_los(Geom,ith,iph)

      ! I/O

      type(Geometry_class), intent(inout):: Geom
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
      if (allocated(Geom%i_scatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%i_scatt)
        deallocate(Geom%i_scatt)
      end if
      if (allocated(Geom%V_CScatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%V_CScatt)
        deallocate(Geom%V_CScatt)
      end if
      if (allocated(Geom%V_SScatt)) then
        MRAMc = MRAMc - 1d-6*sizeof(Geom%V_SScatt)
        deallocate(Geom%V_SScatt)
      end if
      Geom%nScatt = 0

      ! Get trigonometry for LOS
      Co = Geom%L_mu(ith)
      So = sin(Geom%L_theta(ith))

      ! For each quadrature polar angle (input)
      do kk=1,Geom%nTh

        ! Get trigonomet
        Coi = Co*Geom%V_mu(kk)
        Soi = So*sin(Geom%V_theta(kk))

        ! For each quadrature azimuth angle (input)
        do ll=1,Geom%nPh2

          ! Cosine scattering angle
          Ctheta = Coi + Soi*cos(Geom%L_phi(iph) - &
                                 Geom%V_phi(ll))

          ! If cosine exact 1 or overflown
          if (Ctheta.ge.1.0) then

            ! Angle is zero
            theta = 0d0

          ! If cosine exact -1 or overflown
          else if (Ctheta.le.-1.0) then

            ! Angle is 180
            theta = PI

          ! Intermediate values
          else

            ! Arccos
            theta = acos(Ctheta)

          end if ! Control invalid cosines

          ! If not started
          if (box%val.lt.0d0) then

            ! Add to first
            box%val = theta
            Geom%nScatt = 1

          ! List already started
          else

            ! Initialize value finder
            nofound = .True.
            boxr => box

            ! Check till done
            do while (.True.)

              ! If the new angle is close enough to the one
              ! in this box
              if (abs(theta-boxr%val).lt.TINYA) then

                ! Signal found and leave
                nofound = .False.
                exit

              end if ! Close angle

              ! If last one, exit loop
              if (.not.associated(boxr%next)) exit

              ! Point to next
              boxr => boxr%next

            end do ! Till done with the list

            ! If not found, add
            if (nofound) then

              ! Add one new angle
              Geom%nScatt = Geom%nScatt + 1

              ! Allocate next
              allocate(boxr%next)
              boxr => boxr%next
              nullify(boxr%next)

              ! Save value
              boxr%val = theta

            end if ! Add new scattering angle
          end if ! Need to check previous angles

        end do ! Input azimuth
      end do ! Input polar angle

      ! Allocate scattering angles
      allocate(Geom%V_CScatt(Geom%nScatt))
      allocate(Geom%V_SScatt(Geom%nScatt))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_CScatt)
      MRAMc = MRAMc + 1d-6*sizeof(Geom%V_SScatt)

      ! Initialize runner
      boxr => box
      ii = 0

      ! Copy angles
      do while (.True.)

        ! Store
        ii = ii + 1
        Geom%V_CScatt(ii) = boxr%val

        ! If no more boxes, leave
        if (.not.associated(boxr%next)) exit

        ! Shift
        boxr => boxr%next

      end do ! Until everything is copied

      !
      ! Clean boxes
      !

      ! While there is data in boxes
      do while (associated(box%next))

        ! Initialize
        boxr => box%next

        ! If more
        if (associated(boxr%next)) then

          ! Point to second to last
          do while (.True.)

            ! Second-to-last found, leave
            if (.not.associated(boxr%next%next)) exit

            ! Shift one
            boxr => boxr%next

          end do ! Until reached second to last

          ! Remove tail
          deallocate(boxr%next)
          nullify(boxr%next)

        ! No more
        else

          ! Deallocate
          deallocate(boxr)
          nullify(box%next)
          nullify(boxr)

        end if ! More to the right of the box

      end do ! Until the box is clean

      ! Order scattering angles
      call QsortC(Geom%V_CScatt)

      !
      ! Index the directions
      !

      ! Allocate indexing
      allocate(Geom%i_scatt(Geom%nPh2,Geom%nTh,1))

      ! Memory count
      MRAMc = MRAMc + 1d-6*sizeof(Geom%i_scatt)

      ! Initialize running output direction
      mm = 0

      ! Get trigonometry for line of sight
      Co = Geom%L_mu(ith)
      So = sin(Geom%L_theta(ith))

      ! For each quadrature polar angle (input)
      do kk=1,Geom%nTh

        ! Get trigonometry
        Coi = Co*Geom%V_mu(kk)
        Soi = So*sin(Geom%V_theta(kk))

        ! For each quadrature azimuth angle (input)
        do ll=1,Geom%nPh2

          ! Cosine scattering angle
          Ctheta = Coi + Soi*cos(Geom%L_phi(iph) - &
                                 Geom%V_phi(ll))

          ! If cosine exact 1 or overflown
          if (Ctheta.ge.1.0) then

            ! Angle is zero
            theta = 0d0

          ! If cosine exact -1 or overflown
          else if (Ctheta.le.-1.0) then

            ! Angle is 180
            theta = PI

          ! Intermediate values
          else

            ! Arccos
            theta = acos(Ctheta)

          end if ! Control invalid cosines

          !
          ! Find angle position in vector
          !

          ! For every scattering angle
          do nn=1,Geom%nScatt

            ! If same angle
            if (abs(theta-Geom%V_CScatt(nn)).le.TINYA) then

              ! Save index and leave
              Geom%i_scatt(ll,kk,1) = nn
              exit

            end if ! Same angle

          end do ! Scattering angle
        end do ! Quadrature azimuth (input)
      end do ! Quadrature polar (input)

      ! Compute cosines and sines for scattering angles
      Geom%V_SScatt = sin(Geom%V_CScatt)
      Geom%V_CScatt = cos(Geom%V_CScatt)

      ! And done
      return

      end subroutine get_scattering_los

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the angle between two directions\n
      !!   th1(double): Polar angle direction 1\n
      !!   ph1(double): Azimuth angle direction 1\n
      !!   th2(double): Polar angle direction 2\n
      !!   ph2(double): Azimuth angle direction 2
      double precision function atom2lab(th1,ph1,th2,ph2)

      ! I/O

      double precision,intent(in):: th1,ph1,th2,ph2

      ! Local

      double precision:: CTheta


      ! Get cosine between directions
      CTheta = cos(th1)*cos(th2) + &
               sin(th1)*sin(th2)*cos(ph2-ph1)

      ! If over 1 by numerical noise
      if (CTheta.ge.1d0) then

        ! Angle is zero
        atom2lab = 0d0

      ! If below -1 by numerical noise
      elseif (CTheta.le.-1d0) then

        ! Angle is 180º
        atom2lab = pi

      ! Normal case
      else

        ! Get arccos
        atom2lab = acos(CTheta)

      endif

      end function atom2lab

!#####################################################################
!#####################################################################
!#####################################################################

      end module fieldb_mod
