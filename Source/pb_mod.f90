      !> Diagonalization of the Hamiltonian
      module pb_mod
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
!     08/08/2024 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/08/2024:    V3.0.2 - Removed unused variables (TdPA)
!
!     12/12/2023:    V3.0.1 - Added PB0 subroutine (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     12/17/2019:    V1.2.0 - Added options to deal in different ways
!                             with the Zeeman effect (TdPA)
!
!     10/30/2017:    V1.1.1 - Using multilevel flag instead of the S
!                             value (TdPA)
!
!     10/11/2017:    V1.1.0 - Takes into account multi-level if the
!                             spin is 0 (TdPA)
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
!  PB:
!    This subroutine calculates the energy eigenvalues and
!    eigenvectors for a multiplet, in the magnetic-field regime of the
!    incomplete Paschen-Back effect. Diagonality with respect to M is
!    exploited for block-diagonalization
!
!  PB0:
!    Initialize some quantities for zero field.
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use funnj_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Diagonalizes the atomic Hamiltonian.\n
      !!        iz(integer): Height index\n
      !!     iterm(integer): Term index\n
      !!     larmor(dfloat): Magnetic field in larmor frequency
      !!                     units\n
      !!      mode(integer): Type of Zeeman effect\n
      !! Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!   Atom(Atom_class): Structure with the atomic data
      subroutine PB(iz,iterm,larmor,mode,Flgsg,Atom)

      ! I/O

      type(Atom_class),intent(inout):: Atom
      type(Fctsg_class),intent(in):: Flgsg
      integer, intent(in):: iterm, iz, mode
      double precision, intent(in):: larmor

      ! Local

      integer:: nM,iM,iJ,iJ1,idJ,i,i1,j,j1,INFO

      double precision:: rL,S,rM,rJ,rJ1,rJm,rJmin,rJmax
      double precision:: pS,pJ,pJ1,comm
      double precision, dimension(Atom%nJmax,Atom%nJmax):: matr
      double precision, dimension(Atom%nJmax):: diag, odiag
      double precision, dimension(2*(Atom%nJmax-1)):: WORK

      ! Initialize variables
      diag = 0d0
      odiag = 0d0
      matr = 0d0

      ! Get term quantities
      rL = Atom%rLval(iterm)
      S = Atom%Sval(iterm)

      ! Multilevel
      if (Atom%ML) then

        rJ = Atom%rJval(1,iterm)
        nM = nint(2d0*rJ + 1d0)

        ! For each M
        do iM=1,nM

          rM = -rJ + dble(iM-1)

          ! Idenfity the J index for that M (only one)
          Atom%iJval(1,iM,iterm) = 1

          ! Size of the block of this M (just one)
          Atom%nblk(iM,iterm) = 1

          ! Magnetic shift
          Atom%eval(1,iM,iterm,iz) = larmor*rM*Atom%gL(iterm)

          ! Diagonal eigenvector
          Atom%evec(1,1,iM,iterm,iz) = 1d0

        end do

        ! If no splitting (2) or no zeeman (1)
        if (mode.eq.1.or.mode.eq.2) then

          ! For each M
          do iM=1,nM
            Atom%eval(1,iM,iterm,iz) = 0d0
          end do

        end if ! no splitting or no zeeman

      ! Multiterm
      else

        rJmin = abs(rL - S)
        rJmax = rL + S

        pS = S*(S+1d0)*(2d0*S + 1d0)

        nM = nint(2d0*rJmax + 1d0)

        ! Run over the magnetic block
        do iM=1,nM

          rM = -rJmax + dble(iM-1)

          rJm = max(abs(rM),rJmin)

          ! initialize the column index
          i=0

          ! If no Zeeman
          if (mode.eq.1) then

            ! Initialize
            Atom%evec(:,:,iM,iterm,iz) = 0d0

            ! Column loop
            do iJ = 1,Atom%nJ(iterm)

              rJ = Atom%rJval(iJ,iterm)

              if (rJ.ge.rJm) then

                pJ = 2d0*rJ + 1d0

                ! Increase the column index
                i = i + 1

                Atom%iJval(i,iM,iterm) = iJ

                Atom%eval(i,iM,iterm,iz) = Atom%FSfreq(iJ,iterm) - &
                                           Atom%TRfreq(iterm)

                Atom%evec(i,i,iM,iterm,iz) = 1d0

              end if ! Test J compatibility with M

            end do ! iJ

            ! Block dimension
            Atom%nblk(iM,iterm) = i

          ! If any part of Zeeman
          else

            ! Build the Hamiltonian matrix (diagonal and subdiagonal
            ! see DSTEV)
            ! Column loop
            do iJ = 1,Atom%nJ(iterm)

              rJ = Atom%rJval(iJ,iterm)

              if (rJ.ge.rJm) then

                pJ = 2d0*rJ + 1d0

                ! Initialize the row index
                i1 = i
                ! Increase the column index
                i = i + 1

                Atom%iJval(i,iM,iterm) = iJ

                ! Row loop
                do iJ1=iJ,Atom%nJ(iterm)

                  rJ1 = Atom%rJval(iJ1,iterm)

                  if (rJ1.ge.rJm) then

                    idJ = nint(rJ1-rJ)

                    if ((idJ.eq.0).or.(iabs(idJ).eq.1)) then

                      pJ1 = 2d0*rJ1 + 1d0

                      comm = Flgsg%sg(nint(rJmax+rJ+rJ1+rM))* &
                             sqrt(pJ*pJ1*pS)* &
                             fun3j(rJ1,rJ,1d0,-rM,rM,0d0,Flgsg)* &
                             fun6j(rJ1,rJ,1d0,S,S,rL,Flgsg)

                      ! Linear Zeeman effect
                      if (mode.eq.3) then
                        if (idJ.ne.0) comm=.0D0
                      end if

                      ! Increase row index
                      i1 = i1+1

                      if (i1.eq.i) then
                        diag(i) = (Atom%FSfreq(iJ,iterm) - &
                                   Atom%TRfreq(iterm)) + &
                                   larmor*(rM + comm)
                        matr(i,i) = 1d0

                        ! If no splitting (2) or no zeeman (1)
                        if (mode.eq.2) Atom%eval(i,iM,iterm,iz) = &
                                             Atom%FSfreq(iJ,iterm) - &
                                             Atom%TRfreq(iterm)
                      else
                        odiag(i1-1) = larmor*comm
                      end if

                    end if ! \DeltaJ <= 1
                  end if ! rJ > rJm

                end do ! iJ1

              end if ! Test J compatibility with M

            end do ! iJ

            ! Block dimension
            Atom%nblk(iM,iterm) = i

            ! Diagonalize
            if (i.gt.1) call DSTEV('V',i,diag,odiag,matr,Atom%nJmax, &
                                   WORK,INFO)

            ! Store eigenvalues and eigenvectors
            do j=1,i

              ! Only if no no splitting (2)
              if (mode.ne.2) &
                Atom%eval(j,iM,iterm,iz) = diag(j)

              do j1=1,i

                Atom%evec(j1,j,iM,iterm,iz) = matr(j1,j)

                matr(j1,j) = 0d0

              end do ! j1
            end do ! j

          end if ! If zeeman effect
        end do ! iM

      end if ! Multilevel or multiterm

      end subroutine PB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize sizes in atomic Hamiltonian for B=0\n
      !!     iterm(integer): Term index\n
      !!   Atom(Atom_class): Structure with the atomic data
      subroutine PB0(iterm,Atom)

      ! I/O

      type(Atom_class),intent(inout):: Atom
      integer, intent(in):: iterm

      ! Local

      integer:: nM,iM,iJ,i

      double precision:: rL,S,rM,rJ,rJm,rJmin,rJmax

      ! Get term quantities
      rL = Atom%rLval(iterm)
      S = Atom%Sval(iterm)

      ! Multilevel
      if (Atom%ML) then

        rJ = Atom%rJval(1,iterm)
        nM = nint(2d0*rJ + 1d0)

        ! For each M
        do iM=1,nM

          ! Idenfity the J index for that M (only one)
          Atom%iJval(1,iM,iterm) = 1

          ! Size of the block of this M (just one)
          Atom%nblk(iM,iterm) = 1

        end do

      ! Multiterm
      else

        rJmin = abs(rL - S)
        rJmax = rL + S

        nM = nint(2d0*rJmax + 1d0)

        ! Run over the magnetic block
        do iM=1,nM

          rM = -rJmax + dble(iM-1)

          rJm = max(abs(rM),rJmin)

          ! initialize the column index
          i=0

          ! Column loop
          do iJ = 1,Atom%nJ(iterm)

            rJ = Atom%rJval(iJ,iterm)

            if (rJ.ge.rJm) then

              ! Increase the column index
              i = i + 1

              Atom%iJval(i,iM,iterm) = iJ

            end if ! Test J compatibility with M

          end do ! iJ

          ! Block dimension
          Atom%nblk(iM,iterm) = i

        end do ! iM

      end if ! Multilevel or multiterm

      end subroutine PB0

!#####################################################################
!#####################################################################
!#####################################################################

      end module pb_mod
