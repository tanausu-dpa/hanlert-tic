      !> Radiation field tensors
      module jcalc_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     26/04/2017
!  Last version:
!     15/05/2024 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.1 - Generalized declarations of Atom to
!                             allow for empty arrays for any of
!                             them (TdPA)
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
!  Jcalc
!    Calculate the JKQ, JKQS, J00P, JKQC in serial mode
!
!  FInt_line
!    Add contribution to the integrals of the mean radiation field in
!  MPI for bound-bound transitions
!
!  FInt_rest
!    Add contribution to the integrals of the frequency dependent mean
!  radiation field and the mean intensity for bound-free transitions
!  in MPI
!
!  addJKQasym
!    Add ad-hoc radiation field tensors to the total JKQ tensors
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use fieldb_mod
      use math_mod
      use parameters_mod , only : cZero, c2, convF, pi, cSaha, kb, &
                                  fktoJ, TINYB
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the JKQ, JKQS, J00P, JKQC in serial mode\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!        omega(double(:)): Frequency array\n
      !!        Wfreq(double(:)): Frequency trapezoidal weights\n
      !!            lf0(integer): First frequency index for
      !!                          bound-bound transitions\n
      !!            lf1(integer): Last frequency index for bound-bound
      !!                          transitions\n
      !!            pf0(integer): First frequency index for bound-free
      !!                          transitions\n
      !!            pf1(integer): Last frequency index for bound-free
      !!                          transitions\n
      !!               T(double): Temperature\n
      !!              ne(double): Electron number density\n
      !!              WA(double): Angular weight\n
      !!     Stokes(double(:,:)): Stokes parameters\n
      !!       Prof(double(:,:)): Bound-bound normalized line
      !!                          profiles\n
      !!     TSo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                          vertical reference frame\n
      !!    TKQo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                          suitable reference frame\n
      !!    JKQ(dcomplex(:,:,:)): Radiation field tensors
      !!                          integrated over the absorption
      !!                          profile\n
      !!   JKQS(dcomplex(:,:,:)): Radiation field tensors
      !!                          integrated over the emission
      !!                          profile\n
      !!       J00P(double(:,:)): Intensity integrals in the
      !!                          photoionization rates\n
      !!   JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                          frequency dependence\n
      !!         iexu(double(:)): Pre-computed frequency exponential
      subroutine Jcalc(Atom,omega,Wfreq,lf0,lf1,pf0,pf1,T,ne,WA,Stk, &
                       Prof,TSo,TKQo,JKQ,JKQS,J00P,JKQC,iexu)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      integer, intent(in):: lf0,lf1,pf0,pf1
      double precision, intent(in):: T,ne,WA
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), pointer, intent(in):: iexu
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:),intent(inout):: J00P
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TSo
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQo
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout)::JKQS
      complex(kind=8), dimension(0:2,0:2,nfreq), intent(inout):: JKQC

      ! Local

      integer:: Kmax,lKmax
      integer:: ifreq,if0,if1,ia,itran,jtran,K,iQ,iil,jjl,nf

      double precision:: c0,c1,c3,Saha,arg,WF,WFS,W0,W1
      double precision, dimension(:), pointer:: exu

      complex(kind=8), dimension(0:2,0:2,lf0:lf1):: integr


      !
      ! Initializations
      !

      ! Limits for multipoles
      Kmax = min(Krad, 2)
      lKmax = min(Kradl, 2)

      ! Saha non-line-dependent part
      Saha = cSaha*ne/(T**(1.5d0))
      arg = fktoJ/kb/T


      !
      ! Calculate JKQC, the argument of JKQ, and the frequency
      ! exponential
      !

      ! If pre-computed exponential
      if (PIRAM.and.pf1.ge.pf0) then

        ! Point to data
        exu(pf0:pf1) => iexu

      ! Not precomputed but frequencies to calculate
      else if (pf1.ge.pf0) then

        ! Allocate
        allocate(exu(pf0:pf1))

      end if ! Precalculated exponentials

      ! If no pre-computed
      if (.not.PIRAM.or.pf1.lt.pf0) then

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)

          ! Calculate exponential
          exu(ifreq) = diexp(WF)

        end do ! frequencies

      end if ! pre-computed

      ! For each frequency
      do ifreq=1,nfreq

        ! For each K
        do K=0,Kmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if (iQ.ne.0.and.axial) cycle

            ! Integrate JKQC
            JKQC(iQ,K,ifreq) = JKQC(iQ,K,ifreq) + &
                               WA*sum(Stk(:,ifreq)* &
                                      TSo(:,iQ,K))

          end do ! Q
        end do ! K
      end do ! frequencies

      ! For each frequency with lines
      do ifreq=lf0,lf1

        ! For each K
        do K=0,lKmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if (iQ.ne.0.and.axial) cycle

            ! Compute the sum TKQ*Stokes
            integr(iQ,K,ifreq) = sum(Stk(:,ifreq)*TKQo(:,iQ,K))

          end do ! Q
        end do ! K
      end do ! frequencies


      !
      ! Calculate J00 for transitions (b-b and b-f)
      !

      ! Variables for photoion.
      c0 = WA*4d-8*pi/convF

      ! Initialize profile index
      iil = 1

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! Get continuous index
          jtran = itran + Atom(ia)%tshift

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%if0(itran)
          W0 = Atom(ia)%W0(itran)
          if1 = Atom(ia)%if1(itran)
          W1 = Atom(ia)%W1(itran)
          nf = if1 - if0

          ! For each frequency
          do ifreq=if0,if1

            ! If left boundary
            if (ifreq.eq.if0) then

              ! Index
              jjl = iil

              ! Weight
              WF = W0*Prof(jjl,1)*WA

              ! Stimulated
              if (stm) WFS = W0*Prof(jjl,2)*WA

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Index
              jjl = iil + nf

              ! Weight
              WF = W1*Prof(jjl,1)*WA

              ! Stimulated
              if (stm) WFS = W1*Prof(jjl,2)*WA

            ! Normal
            else

              ! Index
              jjl = iil + ifreq - if0

              ! Weight
              WF = Wfreq(ifreq)*Prof(jjl,1)*WA

              ! Stimulated
              if (stm) WFS = Wfreq(ifreq)*Prof(jjl,2)*WA

            end if ! Extrema

            ! For each K
            do K=0,Atom(ia)%Krad(itran)

              ! For each Q
              do iQ=0,K

                ! If axial atmosphere and Q!=0, skip
                if (iQ.ne.0.and.axial) cycle

                ! Add the contribution to the JKQ integral
                JKQ(iQ,K,jtran) = JKQ(iQ,K,jtran) + &
                                  WF*integr(iQ,K,ifreq)

                ! Stimulated
                if(stm) &
                  JKQS(iQ,K,jtran) = JKQS(iQ,K,jtran) + &
                                     WFS*integr(iQ,K,ifreq)

              end do ! Q
            end do ! K
          end do ! frequencies

          ! Advance iil
          iil = iil + nf + 1

        end do ! b-b transitions

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

          ! Continuous index
          jtran = itran + Atom(ia)%pshift

          ! Saha factor
          c3 = Saha*exp(Atom(ia)%phot(itran)%edge*arg)* &
               Atom(ia)%phot(itran)%glu

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%phot(itran)%if0
          W0 = Atom(ia)%phot(itran)%W0
          if1 = Atom(ia)%phot(itran)%if1
          W1 = Atom(ia)%phot(itran)%W1
          nf = if1 - if0

          ! For each frequency
          do ifreq=if0,if1

            ! If left boundary
            if (ifreq.eq.if0) then

              ! Weight
              c1 = W0

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Weight
              c1 = W1

            ! Normal
            else

              ! Weight
              c1 = Wfreq(ifreq)

            end if

            ! Weight with cross section and constants
            c1 = c0*c1*Atom(ia)%phot(itran)%alpha(ifreq)/omega(ifreq)

            ! Weight for J00
            WF = c1*Stk(0,ifreq)

            ! Add the contribution to the integral J00
            J00P(jtran,1) = J00P(jtran,1) + WF
            J00P(jtran,2) = J00P(jtran,2) + WF*exu(ifreq)*c3

          end do ! frequencies
        end do ! b-f transitions
      end do ! atoms

      ! If precomputed exponential
      if (PIRAM.and.pf1.ge.pf0) then

        ! Free pointer
        nullify(exu)

      ! If local allocation
      else if (pf1.ge.pf0) then

        ! Free memory
        deallocate(exu)
        nullify(exu)

      end if ! Exponential allocations

      end subroutine Jcalc

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add contribution to the integrals of the mean radiation field
      !! in MPI for bound-bound transitions\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!         Wfreq(double(:)): Frequency trapezoidal weights\n
      !!             lf0(integer): First frequency index for
      !!                           bound-bound transitions\n
      !!             lf1(integer): Last frequency index for
      !!                           bound-bound transitions\n
      !!            if0p(integer): Lower limit of processor\n
      !!            proc(integer): CPU ID\n
      !!      Stokes(double(:,:)): Stokes parameters\n
      !!        Prof(double(:,:)): Bound-bound normalized line
      !!                           profiles\n
      !!     TKQo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                           suitable reference frame\n
      !!        Norm(double(:,:)): Normalization factor for the
      !!                           line profiles\n
      !!  Bstk(dcomplex(:,:,:,:)): Radiation field tensors partial
      !!                           integrals
      subroutine FInt_line(Atom,Wfreq,lf0,lf1,if0p,proc, &
                           Stk,Prof,TKQo,Norm,Bstk)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      integer, intent(in):: proc,lf0,lf1,if0p
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:), intent(inout):: Norm
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQo
      complex(kind=8), dimension(0:2,0:2,2,nxtran), &
                       intent(inout):: BStk

      ! Local

      integer:: ifreq,ifreqs,ia,itran,jtran,K,iQ,lKmax
      integer:: if0,if1,if0s,if1s,iil

      double precision:: WF,WFS,W0,W1

      complex(kind=8), dimension(0:2,0:2,lf0:lf1):: integr


      !
      ! Initializations
      !

      ! Maximum multipole
      lKmax = min(Kradl, 2)

      ! For each frequency with lines
      do ifreq=lf0,lf1

        ! Shift index
        ifreqs = ifreq - if0p

        ! For each K
        do K=0,lKmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if(iQ.ne.0.and.axial)cycle

            ! Compute the sum TKQ*Stokes
            integr(iQ,K,ifreq) = sum(Stk(:,ifreqs)*TKQo(:,iQ,K))

          end do ! Q
        end do ! K
      end do ! frequencies

      !
      ! Calculate J00 for transitions b-b
      !

      ! Initialize index for profiles
      iil = 0

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(itran)%Mabsent(proc)) cycle

          ! Continuous index
          jtran = itran + Atom(ia)%tshift

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%Mif0(itran,proc)
          if0s = if0 - if0p
          W0 = Atom(ia)%MW0(itran,proc)
          if1 = Atom(ia)%Mif1(itran,proc)
          if1s = if1 - if0p
          W1 = Atom(ia)%MW1(itran,proc)

          !
          ! Boundary points
          !

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W0*Prof(iil,1)
          Norm(1,jtran) = Norm(1,jtran) + WF

          ! Stimulated
          if (stm) then
            WFS = W0*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! For each K
          do K=0,Atom(ia)%Krad(itran)

            ! For each Q
            do iQ=0,K

              ! If axial atmosphere and Q!=0, skip
              if(iQ.ne.0.and.axial)cycle

              ! Add the contribution to the JKQ integral
              BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                   WF*integr(iQ,K,if0)

              ! Stimulated
              if(stm) &
                BStk(iQ,K,2,jtran) = BStk(iQ,K,2,jtran) + &
                                     WFS*integr(iQ,K,if0)

            end do ! Q
          end do ! K

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ! Shift index
            ifreqs = ifreq - if0p

            ! Add the profile of the line to the weights and
            ! the norm integral
            iil = iil + 1
            WF = Wfreq(ifreq)*Prof(iil,1)
            Norm(1,jtran) = Norm(1,jtran) + WF

            ! Stimulated
            if (stm) then
              WFS = Wfreq(ifreq)*Prof(iil,2)
              Norm(2,jtran) = Norm(2,jtran) + WFS
            end if

            ! For each K
            do K=0,Atom(ia)%Krad(itran)

              ! For each Q
              do iQ=0,K

                ! If axial atmosphere and Q!=0, skip
                if(iQ.ne.0.and.axial)cycle

                ! Add the contribution to the JKQ integral
                BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                     WF*integr(iQ,K,ifreq)

                ! Stimulated
                if(stm) BStk(iQ,K,2,jtran) = BStk(iQ,K,2,jtran) + &
                                             WFS*integr(iQ,K,ifreq)

              end do ! Q
            end do ! K
          end do ! frequencies

          ! If single frequency, skip
          if (if1s.le.if0s) cycle

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W1*Prof(iil,1)
          Norm(1,jtran) = Norm(1,jtran) + WF

          ! Stimulated
          if (stm) then
            WFS = W1*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! For each K
          do K=0,Atom(ia)%Krad(itran)

            ! For each Q
            do iQ=0,K

              ! If axial atmosphere and Q!=0, skip
              if(iQ.ne.0.and.axial)cycle

              ! Add the contribution to the JKQ integral
              BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                   WF*integr(iQ,K,if1)

              ! Stimulated
              if(stm) BStk(iQ,K,2,jtran) = BStk(iQ,K,2,jtran) + &
                                           WFS*integr(iQ,K,if1)

            end do ! Q
          end do ! K
        end do ! b-b transitions
      end do ! atoms

      end subroutine FInt_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add contribution to the integrals of the frequency dependent
      !! mean radiation field and the mean intensity for bound-free
      !! transitions in MPI\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       omega(double(:)): Frequency array\n
      !!       Wfreq(double(:)): Frequency trapezoidal weights\n
      !!           pf0(integer): First frequency index for bound-free
      !!                         transitions\n
      !!           pf1(integer): Last frequency index for bound-free
      !!                         transitions\n
      !!          if0l(integer): First frequency index for process\n
      !!          if1l(integer): Last frequency index for process\n
      !!          if0p(integer): Lower limit of processor\n
      !!              T(double): Temperature\n
      !!          proc(integer): CPU ID\n
      !!       Stk(double(:,:)): Stokes parameters\n
      !!    TSo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                         vertical reference frame\n
      !!      J00P(double(:,:)): Intensity integrals in the
      !!                         photoionization rates\n
      !!  JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                         frequency dependence\n
      !!        iexu(double(:)): Pre-computed frequency exponential
      subroutine FInt_rest(Atom,omega,Wfreq,pf0,pf1,if0l,if1l,if0p, &
                           T,proc,WA,Stk,TSo,J00P,JKQC,iexu)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      integer, intent(in):: proc,if0p,pf0,if0l,if1l,pf1
      double precision, intent(in):: T,WA
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), intent(in):: iexu
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:),intent(inout):: J00P
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TSo
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(inout):: JKQC

      ! Local

      integer:: ifreq,ifreqs,ia,itran,jtran,K,iQ,Kmax
      integer:: if0,if1,if0s,if1s

      double precision:: WF,c0,c1,W0,W1
      double precision, dimension(:), allocatable:: exu


      !
      ! Initializations
      !

      ! Maximum multipole
      Kmax = min(Krad, 2)


      !
      ! Calculate JKQC and frequency exponential
      !

      ! Allocate exu
      if (pf1.ge.pf0) allocate(exu(pf0:pf1))

      ! If pre-computed
      if (PIRAM.and.pf1.ge.pf0) then

        ! Copy exponential
        exu = iexu

      ! If no pre-computed
      else

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)

          ! Calculate exponential
          exu(ifreq) = diexp(WF)

        end do ! frequencies

      end if ! pre-computed

      ! For each frequency
      do ifreq=if0l,if1l

        ! Shift index
        ifreqs = ifreq - if0p

        ! For each K
        do K=0,Kmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if(iQ.ne.0.and.axial)cycle

            ! Integrate JKQC
            JKQC(iQ,K,ifreq) = JKQC(iQ,K,ifreq) + &
                               WA*sum(Stk(:,ifreqs)*TSo(:,iQ,K))

          end do ! Q
        end do ! K
      end do ! frequencies


      !
      ! Calculate J00 for transitions (b-b and b-f)
      !

      ! Variables for photoion.
      c0 = WA*4d-8*pi/convF

      ! For each atom
      do ia=1,nA

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(itran)%Mabsent(proc)) cycle

          ! Get continuous index
          jtran = itran + Atom(ia)%pshift

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%phot(itran)%Mif0(proc)
          if0s = if0 - if0p
          W0 = Atom(ia)%phot(itran)%MW0(proc)
          if1 = Atom(ia)%phot(itran)%Mif1(proc)
          if1s = if1 - if0p
          W1 = Atom(ia)%phot(itran)%MW1(proc)

          !
          ! Boundary points
          !

          ! Weight with cross section and constants
          c1 = c0*W0*Atom(ia)%phot(itran)%alpha(if0)/ &
               omega(if0)

          ! Weight for J00
          WF = c1*Stk(0,if0s)

          ! Add the contribution to the integral J00
          J00P(jtran,1) = J00P(jtran,1) + WF
          J00P(jtran,2) = J00P(jtran,2) + WF*exu(if0)

          ! Weight with cross section and constants
          c1 = c0*W1*Atom(ia)%phot(itran)%alpha(if1)/ &
               omega(if1)

          ! Weight for J00
          WF = c1*Stk(0,if1s)

          ! Add the contribution to the integral J00
          J00P(jtran,1) = J00P(jtran,1) + WF
          J00P(jtran,2) = J00P(jtran,2) + WF*exu(if1)

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ! shift index
            ifreqs = ifreq - if0p

            ! Weight with cross section and constants
            c1 = c0*Wfreq(ifreq)*Atom(ia)%phot(itran)%alpha(ifreq)/ &
                 omega(ifreq)

            ! Weight for J00
            WF = c1*Stk(0,ifreqs)

            ! Add the contribution to the integral J00
            J00P(jtran,1) = J00P(jtran,1) + WF
            J00P(jtran,2) = J00P(jtran,2) + WF*exu(ifreq)

          end do ! frequencies
        end do ! b-f transitions
      end do ! atoms

      ! Free
      if (allocated(exu)) deallocate(exu)

      end subroutine FInt_rest

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add ad-hoc radiation field tensors to the total JKQ tensors\n
      !!   Bfield(Bfield_class): Structure with magnetic field data\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                         and J-symbols\n
      !!  JKQa(dcomplex(:,:,:)): Extra asymmetry for the radiation
      !!                         field tensors\n
      !!   JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                         over the absorption profile\n
      !!  JKQS(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                         over the emission profile\n
      !!  JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                         frequency dependence
      subroutine addJKQasym(Bfield,Flgsg,JKQa,JKQ,JKQS,JKQC)

      ! I/O

      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(in):: Flgsg
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(:,:,:,:), intent(inout):: JKQ
      complex(kind=8), dimension(:,:,:,:), intent(inout):: JKQS
      complex(kind=8), dimension(:,:,:,:), intent(inout):: JKQC

      ! Local

      integer:: ifreq,iz,ii,jz


      !
      ! First the continuum, which is in the vertical reference
      ! frame
      !

      ! If forcing the asymmetry
      if (force_asym) then

        ! For each height
        do iz=1,Rnz

          ! For each frequency
          do ifreq=1,nfreq

            ! Add contribution [(3,1,ifreq,iz) is (0,0,ifreq,iz)]
            JKQC(:,2:3,ifreq,iz) = JKQa(:,:,iz)* &
                                   dble(JKQC(3,1,ifreq,iz))

          end do ! frequencies
        end do ! heights

      ! Additive instead
      else

        ! For each height
        do iz=1,Rnz
          ! For each frequency
          do ifreq=1,nfreq

            ! Add contribution [(3,1,ifreq,iz) is (0,0,ifreq,iz)]
            JKQC(:,2:3,ifreq,iz) = JKQC(:,2:3,ifreq,iz) + &
                                   JKQa(:,:,iz)* &
                                   dble(JKQC(3,1,ifreq,iz))

          end do ! frequencies
        end do ! heights

      end if ! Forced or additive

      !
      ! Now the lines, this may need to be rotated
      !
      do iz=1,Rnz

        ! Shift index
        jz = iz + Rz0 - 1

        ! If magnetic field
        if (Bfield%Bstrength(jz).gt.TINYB) then

          ! Rotate to the vertical
          call fieldB(JKQ(:,:,:,iz),nxtran,Flgsg, &
                      -Bfield%Btheta(jz), &
                      -Bfield%Bphi(jz),-1)
        end if

        ! Forcing asymmetry
        if (force_asym) then

          ! For each transition
          do ii=1,nxtran

            ! Add contribution
            JKQ(:,2:3,ii,iz) = JKQa(:,:,iz)*dble(JKQ(3,1,ii,iz))

          end do ! Transitions

        ! Additive instead
        else

          ! For each transition
          do ii=1,nxtran

            ! Add contribution
            JKQ(:,2:3,ii,iz) = JKQ(:,2:3,ii,iz) + &
                               JKQa(:,:,iz)*dble(JKQ(3,1,ii,iz))

          end do ! Transitions

        end if ! Forced or additive

        ! If magnetic field
        if (Bfield%Bstrength(jz).gt.TINYB) then

          ! Rotate back to magnetic field reference frame
          call fieldB(JKQ(:,:,:,iz),nxtran,Flgsg, &
                      Bfield%Btheta(jz),Bfield%Bphi(jz),1)

        end if ! If magnetic field

        ! If stimulated emission
        if (stm) then

          ! If magnetic field
          if (Bfield%Bstrength(jz).gt.TINYB) then

            ! Rotate to the vertical
            call fieldB(JKQS(:,:,:,iz),nxtran,Flgsg, &
                        -Bfield%Btheta(jz), &
                        -Bfield%Bphi(jz),-1)

          end if ! If magnetic field

          ! Forcing asymmetry
          if (force_asym) then

            ! For each transition
            do ii=1,nxtran

              ! Add contribution
              JKQS(:,2:3,ii,iz) = JKQS(:,2:3,ii,iz) + &
                                  JKQa(:,:,iz)*dble(JKQS(3,1,ii,iz))

            end do ! Transitions

          ! Additive instead
          else

            ! For each transition
            do ii=1,nxtran

              ! Add contribution
              JKQS(:,2:3,ii,iz) = JKQS(:,2:3,ii,iz) + &
                                  JKQa(:,:,iz)*dble(JKQS(3,1,ii,iz))

            end do ! Transitions

          end if ! Forced or additive

          ! If magnetic field
          if (Bfield%Bstrength(jz).gt.TINYB) then

            ! Rotate back to magnetic field reference frame
            call fieldB(JKQS(:,:,:,iz),nxtran,Flgsg, &
                        Bfield%Btheta(iz),Bfield%Bphi(iz),1)

          end if ! Magnetic field
        end if ! stimulated emission

      end do ! Heights

      end subroutine addJKQasym

!#####################################################################
!#####################################################################
!#####################################################################

      end module jcalc_mod
