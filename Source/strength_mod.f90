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
!     10/26/2022 V3.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
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
!    This subroutine calculates the component strengths
!  (electric-dipole elements) for all possible transitions
!  in the atom
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use funnj_mod
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
      integer:: maxnM,maxnJ
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

          maxnM = max(nM,nM1)
          maxnJ = max(Atom%nJ(iterm),Atom%nJ(iterm1))

          ! Allocate second level of dipole array
          allocate(Atom%rdip(itran)%rdip(-1:1,maxnM,maxnM, &
                                         maxnJ,maxnJ))
          Atom%rdip(itran)%rdip = 0d0

          ! LS configurations
          ! For each level of the lower term
          do iJ=1,Atom%nJ(iterm)

            rJ = Atom%rJval(iJ,iterm)

            pJ=sqrt(2d0*rJ + 1d0)

            ! For each level of the upper term
            do iJ1=1,Atom%nJ(iterm1)

              rJ1 = Atom%rJval(iJ1,iterm1)

              if (nint(abs(rJ1-rJ)).gt.1) cycle

              CCJ = pJ*sqrt(2d0*rJ1 + 1d0)* &
                    fun6j(rL,rL1,1d0,rJ1,rJ,S,Flgsg)

              ! Azimuthal components
              ! For each magnetic component of the lower level
              do iM=1,nM

                rM = -rJmax + dble(iM-1)

                ! For each magnetic component of the upper level
                do iM1=1,nM1

                  rM1 = -rJ1max + dble(iM1-1)

                  q=rM-rM1
                  iq=nint(q)

                  if (abs(iq).gt.1) cycle

                  Atom%rdip(itran)%rdip(iq,iM,iM1,iJ,iJ1) = &
                                   Flgsg%sg(nint(rL+S-rM))*CCJ* &
                                   fun3j(rJ,rJ1,1d0,-rM,rM1,q,Flgsg)

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

      end module strength_mod
