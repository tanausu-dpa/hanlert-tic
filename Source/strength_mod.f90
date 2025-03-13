      !> Transition strengths
      module strength_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     18/04/2017
!  Last version:
!     20/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     20/12/2024:    V4.0.0 - Always use the dipole strength in the
!                             energy eigenbasis (TdPA)
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
!  strength
!    Calculate the electric dipole transition strength for a given
!  atom
!
!  strength_ev
!    Calculate the electric dipole transition strength in the energy
!  eigenbasis for a given atom
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use aborted_mod
      use funnj_mod
      use parameters_mod , only : TINYB , TINYEV
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the electric dipole transition strength for a given
      !! atom\n
      !!    Atom(Atom_class): Structure with atomic data\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      subroutine strength(Atom,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: itran,iterm,iterm1,iJ,iJ1,nM,nM1,iM,iM1,iq

      double precision:: rL,rL1,S,rJ,rJ1,rJmax,rJ1max
      double precision:: rM,rM1,q,pL,pL1,pJ,CCJ


      ! Routine name
      urou = 'strength'

      ! If already allocated, skip
      if (allocated(Atom%rdip)) return

      ! Allocate first level of dipole strength
      allocate(Atom%rdip(Atom%ntran))

      ! For each lower term
      do iterm=1,Atom%nMulti-1

        ! Orbitan and spin momenta
        rL = Atom%rLval(iterm)
        S = Atom%Sval(iterm)

        ! Maximum J value
        rJmax = rL + S

        ! Number of magnetic sublevels
        nM = nint(2d0*rJmax + 1d0)

        ! Orbital momentum factor
        pL = sqrt(2d0*rL + 1d0)

        ! For each upper term
        do iterm1=iterm+1,Atom%nMulti

          ! Check if there is a transition
          itran = Atom%irad(iterm,iterm1)

          ! Skip if no transition
          if (itran.eq.0) cycle

          ! Get orbitan angular momentum
          rL1 = Atom%rLval(iterm1)

          ! Get maximum J value
          rJ1max = rL1 + S

          ! Get number of magnetic sublevels
          nM1 = nint(2d0*rJ1max + 1d0)

          ! Orbital momentum factor
          pL1 = sqrt(2d0*rL1 + 1d0)

          ! Add to memory
          MRAMc = MRAMc + 1d-6*sizeof(Atom%rdip(itran))

          ! Allocate and initialize second level of dipole array
          allocate(Atom%rdip(itran)%rdip(-1:1,nM1,nM, &
                                         Atom%nJ(iterm1), &
                                         Atom%nJ(iterm)))
          MRAMc = MRAMc + 1d-6*sizeof(Atom%rdip(itran)%rdip)
          Atom%rdip(itran)%rdip = 0d0

          !
          ! LS configurations
          !

          ! For each level of the lower term
          do iJ=1,Atom%nJ(iterm)

            ! Get total angular momentum
            rJ = Atom%rJval(iJ,iterm)

            ! Total angular momentum factor
            pJ = sqrt(2d0*rJ + 1d0)

            ! For each level of the upper term
            do iJ1=1,Atom%nJ(iterm1)

              ! Get total angular momentum
              rJ1 = Atom%rJval(iJ1,iterm1)

              ! Check electric dipole rule
              if (abs(rJ1-rJ).gt.1.1d0) cycle

              ! Get 6J
              CCJ = pJ*sqrt(2d0*rJ1 + 1d0)* &
                    fun6j(rL,rL1,1d0,rJ1,rJ,S,Flgsg)

              !
              ! Azimuthal components
              !

              ! For each magnetic component of the lower level
              do iM=1,nM

                ! Get magnetic quantum number
                rM = -rJmax + dble(iM-1)

                ! For each magnetic component of the upper level
                do iM1=1,nM1

                  ! Get magnetic quantum number
                  rM1 = -rJ1max + dble(iM1-1)

                  ! Get difference in magnetic quantum number
                  q = rM - rM1
                  iq = nint(q)

                  ! Selection rule
                  if (abs(iq).gt.1) cycle

                  ! Make use of conjugation properties of matrix
                  ! elements
                  Atom%rdip(itran)%rdip(-iq,iM1,iM,iJ1,iJ) = &
                                   Flgsg%sg(nint(rL1+S-rM1))*CCJ* &
                                   fun3j(rJ1,rJ,1d0,-rM1,rM,-q,Flgsg)

                end do ! iM1
              end do ! iM
            end do ! iJ1
          end do ! iJ
        end do ! iterm1
      end do ! iterm

      ! Check if everything is fine
      call control

      return

      end subroutine strength

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the electric dipole transition strength in the
      !! energy eigenbasis for a given atom\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!  Bfield(Bfield_class): Structure with magnetic field data
      subroutine strength_ev(Atom,Bfield)

      ! I/O

      type(Bfield_class), intent(in):: Bfield
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iz,itran,iterm,iterm1,maxnMu,maxnMl,maxnju,maxnjl
      integer:: iJ,iJ1,nM,nM1,iM,iM1,iq,jM,jM1,kM,kM1

      double precision:: rL,rL1,S,rJ,rJ1,rJmax,rJ1max,rM,rM1,q,cM,cM1


      ! Routine name
      urou = 'strength_ev'

      ! Allocate first level of dipole strength
      allocate(Atom%rdipev(Rz0:Rz1))

      ! For each height
      do iz=Rz0,Rz1

        ! Count memory
        MRAMc = MRAMc + 1d-6*sizeof(Atom%rdipev(iz))

        ! If no magnetic field, skip
        if (Bfield%Bstrength(iz).le.TINYB) cycle

        ! Allocate for transitions
        allocate(Atom%rdipev(iz)%rdipev(Atom%ntran))

        ! For each lower term
        do iterm=1,Atom%nMulti-1

          ! Get orbital and spin momenta
          rL = Atom%rLval(iterm)
          S = Atom%Sval(iterm)

          ! Get maximum value of total angular momentum
          rJmax = rL + S

          ! Get number of magnetic components
          nM = nint(2d0*rJmax + 1d0)

          ! For each upper term
          do iterm1=iterm+1,Atom%nMulti

            ! Check if there is a transition
            itran = Atom%irad(iterm,iterm1)

            ! Skip if no transition
            if (itran.eq.0) cycle

            ! Count memory
            MRAMc = MRAMc + 1d-6*sizeof(Atom%rdipev(iz)%rdipev(itran))

            ! Get orbital momentum
            rL1 = Atom%rLval(iterm1)

            ! Get maximum total angular momentum
            rJ1max = rL1 + S

            ! Get number of magnetic components
            nM1 = nint(2d0*rJ1max + 1d0)

            ! Maximum sizes
            maxnMu = nM1
            maxnMl = nM
            maxnju = maxval(Atom%nblk(1:nM1,iterm1))
            maxnjl = maxval(Atom%nblk(1:nM,iterm))

            ! Allocate second level of dipole array
            allocate(Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(-1:1,maxnju,maxnjl,maxnMu,maxnMl))
            MRAMc = MRAMc + 1d-6*sizeof(Atom%rdipev(iz)% &
                                             rdipev(itran)%rdip)

            ! For each magnetic component of the lower level
            do iM=1,nM

              ! Get magnetic quantum number
              rM = -rJmax + dble(iM-1)

              ! For each magnetic component of the upper level
              do iM1=1,nM1

                ! Get magnetic quantum number
                rM1 = -rJ1max + dble(iM1-1)

                ! Get change in magnetic quantum number
                q = rM - rM1
                iq = nint(q)

                ! Selection rule
                if (abs(iq).gt.1) cycle

                ! For each lower level state
                do jM=1,Atom%nblk(iM,iterm)

                  ! For each upper level state
                  do jM1=1,Atom%nblk(iM1,iterm1)

                    ! Initialize
                    Atom%rdipev(iz)%rdipev(itran)% &
                         rdip(-iq,jM1,jM,iM1,iM) = 0d0

                    ! Run over lower level block
                    do kM=1,Atom%nblk(iM,iterm)

                      ! Get eigenvector
                      cM = Atom%evec(kM,jM,iM,iterm,iz)

                      ! If coefficient too small, skip
                      if (abs(cM).lt.TINYEV) cycle

                      ! Get J info
                      iJ = Atom%iJval(kM,iM,iterm)
                      rJ = Atom%rJval(iJ,iterm)

                      ! Run over upper level block
                      do kM1=1,Atom%nblk(iM1,iterm1)

                        ! Get eigenvector
                        cM1 = Atom%evec(kM1,jM1,iM1,iterm1,iz)

                        ! If coefficient too small, skip
                        if (abs(cM1).lt.TINYEV) cycle

                        ! Get J info
                        iJ1 = Atom%iJval(kM1,iM1,iterm1)
                        rJ1 = Atom%rJval(iJ1,iterm1)

                        ! Selection rule
                        if (abs(rJ1-rJ).gt.1.1d0) cycle

                        ! Add to strength
                        Atom%rdipev(iz)%rdipev(itran)% &
                             rdip(-iq,jM1,jM,iM1,iM) = &
                                      Atom%rdipev(iz)%rdipev(itran)% &
                                           rdip(-iq,jM1,jM,iM1,iM) + &
                                      cM*cM1*Atom%rdip(itran)% &
                                               rdip(-iq,iM1,iM,iJ1,iJ)
                      end do ! kM1
                    end do ! kM
                  end do ! jM1
                end do ! jM
              end do ! iM1
            end do ! iM
          end do ! Upper term
        end do ! Lower term
      end do ! height

      ! Check if everything is fine
      call control
      return

      end subroutine strength_ev

!#####################################################################
!#####################################################################
!#####################################################################

      end module strength_mod
