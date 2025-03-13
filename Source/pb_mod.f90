      !> Diagonalization of the Hamiltonian
      module pb_mod
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
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Revised headers (TdPA)
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
!  PB
!    Diagonalize the atomic Hamiltonian
!
!  PB0
!    Initialize quantities related to diagonalization when there is no
!  magnetic field
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

      !> Diagonalize the atomic Hamiltonian\n
      !!         iz(integer): Current height index\n
      !!      iterm(integer): Current term index\n
      !!      larmor(dfloat): Magnetic field in larmor frequency
      !!                      units\n
      !!       mode(integer): Type of Zeeman effect\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!    Atom(Atom_class): Structure with atomic data
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

        ! Get angular momentum and number of magnetic sublevels
        rJ = Atom%rJval(1,iterm)
        nM = nint(2d0*rJ + 1d0)

        ! For each M
        do iM=1,nM

          ! Get M value
          rM = -rJ + dble(iM-1)

          ! Idenfity the J index for that M (only one)
          Atom%iJval(1,iM,iterm) = 1

          ! Size of the block of this M (just one)
          Atom%nblk(iM,iterm) = 1

          ! Magnetic shift from Linear Zeeman
          Atom%eval(1,iM,iterm,iz) = larmor*rM*Atom%gL(iterm)

          ! Diagonal eigenvector
          Atom%evec(1,1,iM,iterm,iz) = 1d0

        end do ! Magnetic sublevels

        ! If no splitting (2) or no zeeman (1)
        if (mode.eq.1.or.mode.eq.2) then

          ! For each M
          do iM=1,nM

            ! Make displacement equal to zero
            Atom%eval(1,iM,iterm,iz) = 0d0

          end do ! Magnetic sublevels

        end if ! no splitting or no zeeman

      ! Multiterm
      else

        ! Get angular momentum limits
        rJmin = abs(rL - S)
        rJmax = rL + S

        ! Get Spin factor
        pS = S*(S+1d0)*(2d0*S + 1d0)

        ! Number of unique magnetic quantum numbers
        nM = nint(2d0*rJmax + 1d0)

        ! Run over the magnetic number blocks
        do iM=1,nM

          ! Get value of magnetic quantum number
          rM = -rJmax + dble(iM-1)

          ! Get minimum possible value of angular momentum
          rJm = max(abs(rM),rJmin)

          ! initialize the column index
          i = 0

          ! If no Zeeman
          if (mode.eq.1) then

            ! Initialize eigenvectors to zero
            Atom%evec(:,:,iM,iterm,iz) = 0d0

            ! For each J level in the term
            do iJ=1,Atom%nJ(iterm)

              ! Get angular momentum
              rJ = Atom%rJval(iJ,iterm)

              ! If valid for the magnetic quantum number
              if (rJ.ge.rJm) then

                ! Get angular momentum factor
                pJ = 2d0*rJ + 1d0

                ! Increase the column index
                i = i + 1

                ! Save index for the position in magnetic block
                Atom%iJval(i,iM,iterm) = iJ

                ! No magnetic shift
                Atom%eval(i,iM,iterm,iz) = Atom%FSfreq(iJ,iterm) - &
                                           Atom%TRfreq(iterm)

                ! Make diagonal
                Atom%evec(i,i,iM,iterm,iz) = 1d0

              end if ! Test J compatibility with M

            end do ! J levels

            ! Block dimension
            Atom%nblk(iM,iterm) = i

          ! If any part of Zeeman
          else

            !
            ! Build the Hamiltonian matrix (diagonal and subdiagonal
            ! see DSTEV)
            !

            ! For each level in the term
            ! Column loop
            do iJ=1,Atom%nJ(iterm)

              ! Get angular momentum
              rJ = Atom%rJval(iJ,iterm)

              ! If valid for this magnetic quantum number
              if (rJ.ge.rJm) then

                ! Get angular momentum factor
                pJ = 2d0*rJ + 1d0

                ! Initialize the row index
                i1 = i

                ! Increase the column index
                i = i + 1

                ! Save index for the position in magnetic block
                Atom%iJval(i,iM,iterm) = iJ

                ! For each level in the term
                ! Row loop
                do iJ1=iJ,Atom%nJ(iterm)

                  ! Get angular momentum
                  rJ1 = Atom%rJval(iJ1,iterm)

                  ! If valid for this magnetic quantum number
                  if (rJ1.ge.rJm) then

                    ! J difference
                    idJ = nint(rJ1-rJ)

                    ! If difference would allow electric dipole
                    if ((idJ.eq.0).or.(iabs(idJ).eq.1)) then

                      ! Get angular momentum factor
                      pJ1 = 2d0*rJ1 + 1d0

                      ! Get Hamiltonian element
                      comm = Flgsg%sg(nint(rJmax+rJ+rJ1+rM))* &
                             sqrt(pJ*pJ1*pS)* &
                             fun3j(rJ1,rJ,1d0,-rM,rM,0d0,Flgsg)* &
                             fun6j(rJ1,rJ,1d0,S,S,rL,Flgsg)

                      ! If linear zeeman effect and interference,
                      ! neglect it
                      if (mode.eq.3.and.iDj.ne.0) comm=.0D0

                      ! Increase row index
                      i1 = i1 + 1

                      ! If diagonal
                      if (i1.eq.i) then

                        ! Store in diagonal
                        diag(i) = (Atom%FSfreq(iJ,iterm) - &
                                   Atom%TRfreq(iterm)) + &
                                   larmor*(rM + comm)

                        ! Initialize
                        matr(i,i) = 1d0

                        ! If no splitting (2) or no zeeman (1), get
                        ! eigenvalue with no magnetic shift
                        if (mode.eq.2) Atom%eval(i,iM,iterm,iz) = &
                                             Atom%FSfreq(iJ,iterm) - &
                                             Atom%TRfreq(iterm)

                      ! Not diagonal
                      else

                        ! Store off-diagonal
                        odiag(i1-1) = larmor*comm

                      end if ! Diagonal element
                    end if ! \DeltaJ <= 1
                  end if ! rJ > rJm

                end do ! Levels in term (iJ1)

              end if ! Test J compatibility with M

            end do ! Levels in term (iJ)

            ! Block dimension
            Atom%nblk(iM,iterm) = i

            ! Diagonalize
            if (i.gt.1) &
              call DSTEV('V',i,diag,odiag,matr,Atom%nJmax,WORK,INFO)

            !
            ! Store eigenvalues and eigenvectors
            !

            ! For each element in the block
            do j=1,i

              ! Only if no no-splitting (2), get eigenvalue
              if (mode.ne.2) &
                Atom%eval(j,iM,iterm,iz) = diag(j)

              ! For each element in the block
              do j1=1,i

                ! Get eigenvector
                Atom%evec(j1,j,iM,iterm,iz) = matr(j1,j)

                ! And reset value
                matr(j1,j) = 0d0

              end do ! Element in the block (j1)
            end do ! Element in the block (j)

          end if ! If zeeman effect

        end do ! iM

      end if ! Multilevel or multiterm

      end subroutine PB

!#####################################################################
!#####################################################################
!#####################################################################

      !> Initialize quantities related to diagonalization when there
      !! is no magnetic field\n
      !!    iterm(integer): Current term index\n
      !!  Atom(Atom_class): Structure with atomic data
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

        ! Get angular momentum
        rJ = Atom%rJval(1,iterm)

        ! Get number of magnetic sublevels
        nM = nint(2d0*rJ + 1d0)

        ! For each M
        do iM=1,nM

          ! Idenfity the J index for that M (only one)
          Atom%iJval(1,iM,iterm) = 1

          ! Size of the block of this M (just one)
          Atom%nblk(iM,iterm) = 1

        end do ! Magnetic sublevels

      ! Multiterm
      else

        ! Get limits of angular momentum
        rJmin = abs(rL - S)
        rJmax = rL + S

        ! Get number of magnetic components
        nM = nint(2d0*rJmax + 1d0)

        ! Run over the magnetic blocks
        do iM=1,nM

          ! Get magnetic quantum number
          rM = -rJmax + dble(iM-1)

          ! Get minimum valid angular momentum
          rJm = max(abs(rM),rJmin)

          ! initialize the column index
          i = 0

          ! For every level in the term
          ! Column loop
          do iJ=1,Atom%nJ(iterm)

            ! Get angular momentum
            rJ = Atom%rJval(iJ,iterm)

            ! If valid for this magnetic quantum number
            if (rJ.ge.rJm) then

              ! Increase the column index
              i = i + 1

              ! Save index
              Atom%iJval(i,iM,iterm) = iJ

            end if ! Test J compatibility with M

          end do ! Levels (iJ)

          ! Block dimension
          Atom%nblk(iM,iterm) = i

        end do ! Magnetic blocks (iM)

      end if ! Multilevel or multiterm

      end subroutine PB0

!#####################################################################
!#####################################################################
!#####################################################################

      end module pb_mod
