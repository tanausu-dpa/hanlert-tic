      !> Transition strengths
      module strength_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/18/2017
!  Last version:
!     04/02/2024 V3.1.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     04/02/2024:    V3.1.1 - Limit the search in nblk to the valid
!                             ones to find the maximum size for all
!                             M (TdPA)
!                           - Bugfix: Wrong operatior in if clauses
!                             to check magnitude of eigenvectors in
!                             strength_ev (TdPA)
!
!     04/01/2024:    V3.1.0 - Removed unused conjugates in rdip,
!                             restricting the allocation to what is
!                             actually needed (TdPA)
!                           - Changed a comparison between integers
!                             with "nint" to one between floats (TdPA)
!                           - Added strength_ev routine (TdPA)
!
!     10/26/2022:    V3.0.1 - Changed the storage structure of the
!                             rdip variable (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     11/19/2019:    V1.1.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     04/18/2017:    V1.0.0 - First version (TdPA)
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
!  strength:
!    This subroutine calculates the component strengths
!  (electric-dipole elements) for all possible transitions
!  in the atom
!
!  strength_ev:
!    This subroutine calculates the component strengths
!  (electric-dipole elements) for all possible transitions
!  in the atom in the energy representation
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

      !> Computes transition strengths for the radiation transfer
      !! coefficients\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      subroutine strength(Atom,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: itran,iterm,iterm1
      integer:: iJ,iJ1,nM,nM1,iM,iM1,iq

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

        rL = Atom%rLval(iterm)
        S = Atom%Sval(iterm)

        rJmax = rL + S

        nM = nint(2d0*rJmax + 1d0)

        pL = sqrt(2d0*rL + 1d0)

        ! and each upper term
        do iterm1=iterm+1,Atom%nMulti

          ! Check if there is a transition
          itran = Atom%irad(iterm,iterm1)

          if (itran.eq.0) cycle

          rL1 = Atom%rLval(iterm1)

          rJ1max = rL1 + S

          nM1 = nint(2d0*rJ1max + 1d0)

          pL1=sqrt(2d0*rL1 + 1d0)

          ! Allocate second level of dipole array
          allocate(Atom%rdip(itran)%rdip(-1:1,nM1,nM, &
                                         Atom%nJ(iterm1), &
                                         Atom%nJ(iterm)))
          Atom%rdip(itran)%rdip = 0d0

          ! LS configurations
          ! For each level of the lower term
          do iJ=1,Atom%nJ(iterm)

            rJ = Atom%rJval(iJ,iterm)

            pJ=sqrt(2d0*rJ + 1d0)

            ! For each level of the upper term
            do iJ1=1,Atom%nJ(iterm1)

              rJ1 = Atom%rJval(iJ1,iterm1)

              if (abs(rJ1-rJ).gt.1.1d0) cycle

              CCJ = pJ*sqrt(2d0*rJ1 + 1d0)* &
                    fun6j(rL,rL1,1d0,rJ1,rJ,S,Flgsg)

              ! Azimuthal components
              ! For each magnetic component of the lower level
              do iM=1,nM

                rM = -rJmax + dble(iM-1)

                ! For each magnetic component of the upper level
                do iM1=1,nM1

                  rM1 = -rJ1max + dble(iM1-1)

                  q = rM - rM1
                  iq = nint(q)

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

#ifdef RDIPEV
!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes transition strengths for the radiation transfer
      !! coefficients in energy basis\n
      !!     Atom(Atom_class): Structure with the atomic data\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      !!      Bfield(Bfield_class): Structure with magnetic field
      subroutine strength_ev(Atom,Flgsg,Bfield)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Atom_class), intent(inout):: Atom

      ! Local

      integer:: iz,itran,iterm,iterm1
      integer:: maxnMu,maxnMl,maxnju,maxnjl
      integer:: iJ,iJ1,nM,nM1,iM,iM1,iq
      integer:: jM,jM1,kM,kM1

      double precision:: rL,rL1,S,rJ,rJ1,rJmax,rJ1max
      double precision:: rM,rM1,q,pJ,cM,cM1


      ! Routine name
      urou = 'strength_ev'

      ! Allocate first level of dipole strength
      allocate(Atom%rdipev(Rz0:Rz1))

      ! For each height
      do iz=Rz0,Rz1

        ! If no magnetic field, skip
        if (Bfield%Bstrength(iz).le.TINYB) cycle

        ! Allocate for transitions
        allocate(Atom%rdipev(iz)%rdipev(Atom%ntran))

        ! For each lower term
        do iterm=1,Atom%nMulti-1

          rL = Atom%rLval(iterm)
          S = Atom%Sval(iterm)

          rJmax = rL + S

          nM = nint(2d0*rJmax + 1d0)

          ! and each upper term
          do iterm1=iterm+1,Atom%nMulti

            ! Check if there is a transition
            itran = Atom%irad(iterm,iterm1)

            if (itran.eq.0) cycle

            rL1 = Atom%rLval(iterm1)

            rJ1max = rL1 + S

            nM1 = nint(2d0*rJ1max + 1d0)

            ! Maximum sizes
            maxnMu = nM1
            maxnMl = nM
            maxnju = maxval(Atom%nblk(1:nM1,iterm1))
            maxnjl = maxval(Atom%nblk(1:nM,iterm))

            ! Allocate second level of dipole array
            allocate(Atom%rdipev(iz)%rdipev(itran)% &
                          rdip(-1:1,maxnju,maxnjl,maxnMu,maxnMl))

            ! For each magnetic component of the lower level
            do iM=1,nM

              ! Get M
              rM = -rJmax + dble(iM-1)

              ! For each magnetic component of the upper level
              do iM1=1,nM1

                rM1 = -rJ1max + dble(iM1-1)

                q = rM - rM1
                iq = nint(q)

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

#endif
!#####################################################################
!#####################################################################
!#####################################################################

      end module strength_mod
