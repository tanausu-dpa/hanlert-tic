      !> Mean intensities
      module jcalci_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     26/04/2017
!  Last version:
!     18/03/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/03/2025:    V4.0.1 - Added the option to skip the calculation
!                             of the lambda operator for bound-free
!                             transitions (TdPA)
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
!  JcalcI
!    Calculate the J00, J00S, J00P, J00C, and lambda operators in
!  serial mode
!
!  FIntI_line
!    Add contribution to the integrals of the mean intensity in MPI
!  for bound-bound transitions
!
!  FIntI_rest
!    Add contribution to the integrals of the frequency dependent mean
!  intensity and the mean intensity and lambda operator for bound-free
!  transitions in MPI
!
!  JcalcJ
!    Calculate the frequency dependent mean intensity in serial mode
!
!  FIntJ
!    Calculate the contribution to the frequency dependent mean
!  intensity in MPI
!
!  Jgen
!    Calculate the K=Q=0 component of the integrated radiation field
!  tensors in serial mode
!
!  FJgInt
!    Calculate the contribution to K=Q=0 component of the integrated
!  radiation field tensors in MPI
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use math_mod
      use parameters_mod , only : cZero, c2, convF, pi, cSaha, kb, &
                                  fktoJ
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the J00, J00S, J00P, J00C, and lambda operators in
      !! serial mode\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!        omega(double(:)): Frequency array\n
      !!        Wfreq(double(:)): Frequency trapezoidal weights\n
      !!            pf0(integer): First frequency index for bound-free
      !!                          transitions\n
      !!            pf1(integer): Last frequency index for bound-free
      !!                          transitions\n
      !!               T(double): Temperature\n
      !!              ne(double): Electron number density\n
      !!            iph(integer): Direction azimuth index\n
      !!            ith(integer): Direction polar index\n
      !!          Stk(double(:)): Intensity\n
      !!        rLine(double(:)): Bound-bound transition strength\n
      !!        rPhot(double(:)): Bound-free transition strength\n
      !!       Prof(double(:,:)): Bound-bound normalized line
      !!                          profiles\n
      !!          J00(double(:)): Mean intensity integrated over the
      !!                          absorption profile\n
      !!         J00S(double(:)): Mean intensity integrated over the
      !!                          emission profile\n
      !!       J00P(double(:,:)): Intensity integrals in the
      !!                          photoionization rates\n
      !!         J00C(double(:)): Mean intensity with frequency
      !!                          dependence\n
      !!    LambdaL(double(:,:)): Bound-bound Lambda operator\n
      !!  LambdaP(double(:,:,:)): Bound-free Lambda operator\n
      !!            ALI(logical): If computing Lambda operator\n
      !!           ALIp(logical): If computing Lambda operator for
      !!                          photoionizations\n
      !!         iexu(double(:)): Pre-computed frequency exponential
      subroutine JcalcI(Atom,Geom,omega,Wfreq,pf0,pf1,T,ne,iph,ith, &
                        Stk,rLine,rPhot,Prof,J00,J00S,J00P,J00C, &
                        LambdaL,LambdaP,ALI,ALIp,iexu)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      logical, intent(in):: ALI,ALIp
      integer, intent(in):: ith,iph,pf0,pf1
      double precision, intent(in):: T, ne
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(inout):: J00C
      double precision, dimension(:), intent(in), pointer:: iexu
      double precision, dimension(:), intent(inout):: J00
      double precision, dimension(:), intent(inout):: J00S
      double precision, dimension(:), intent(in):: rLine
      double precision, dimension(:), intent(in):: rPhot
      double precision, dimension(:,:),intent(inout):: J00P
      double precision, dimension(:,:), intent(inout):: LambdaL
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:,:),intent(inout):: LambdaP

      ! Local

      logical:: nALI,nALIp

      integer:: iil,iip,jjl,jjp,nf
      integer:: ifreq,if0,if1,ia,itran,jtran,ftran,jftran

      double precision:: c0,c1,c3,Saha,arg,WA,WF,WFS,W0,W1

      ! Pointers

      double precision, dimension(:), pointer:: exu


      !
      ! Initializations
      !

      ! Pointer
      nullify(exu)

      ! Saha non-line-dependent part
      Saha = cSaha*ne/(T**(1.5d0))
      arg = fktoJ/kb/T

      ! Angular weight
      WA = Geom%W_mu(ith)*Geom%W_mux(iph)

      ! Not ALI
      nALI = .not.ALI
      nALIp = .not.(ALI.and.ALIp)

      ! If pre-computed exponentials
      if (PIRAM.and.pf1.ge.pf0) then

        ! Point
        exu(pf0:pf1) => iexu

      ! If non-allocated but needed
      else if (pf1.ge.pf0) then

        ! Allocate
        allocate(exu(pf0:pf1))

      end if ! Allocated exponential


      !
      ! Calculate J00C and frequency exponential
      !

      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreq)

      end do ! Frequencies

      ! If no pre-computed
      if (.not.PIRAM.or.pf1.lt.pf0) then

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)

          ! Inverse exponential
          exu(ifreq) = diexp(WF)

        end do ! frequencies

      end if ! pre-computed


      !
      ! Calculate J00 for transitions (b-b and b-f)
      !

      ! Variables for photoion.
      c0 = WA*4d-8*pi/convF

      ! Initialize indexes for profiles
      iil = 1
      iip = 1

      ! For each atom
      do ia=1,nA

        ! For each FS transition
        do ftran=1,Atom(ia)%nftran

          ! Transition indexes
          itran = Atom(ia)%ifst(ftran)
          jftran = ftran + Atom(ia)%tfshift

          ! Get limits data, limit weights, and size
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
              WF = W0*Prof(iil,1)*WA

             !! Stimulated emission
             !if (stm) WFS = W0*Prof(jjl,2)*WA

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Index
              jjl = iil + nf

              ! Weight
              WF = W1*Prof(jjl,1)*WA

             !! Stimulated emission
             !if (stm) WFS = W1*Prof(jjl,2)*WA

            ! Normal
            else

              ! Index
              jjl = iil + ifreq - if0

              ! Weight
              WF = Wfreq(ifreq)*Prof(jjl,1)*WA

             !! Stimulated emission
             !if (stm) WFS = Wfreq(ifreq)*Prof(jjl,2)*WA

            end if

            ! Add the contribution to the integrals J00
            J00(jftran) = J00(jftran) + WF*Stk(ifreq)

           !! Stimulated emission
           !if(stm) J00S(jftran) = J00S(jftran) + WFS*Stk(ifreq)

            ! No ALI in iteration
            if (nALI) cycle

            ! Add contribution to Lambda operator
            LambdaL(1,jftran) = LambdaL(1,jftran) + WF*rLine(jjl)

          end do ! Frequencies

          ! Advance iil
          iil = iil + nf + 1

        end do ! b-b transitions

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

          ! Get continuous index
          jtran = itran + Atom(ia)%pshift

          ! Saha factor
          c3 = Saha*exp(Atom(ia)%phot(itran)%edge*arg)* &
               Atom(ia)%phot(itran)%glu

          ! Get limits data, limit weights, and size
          if0 = Atom(ia)%phot(itran)%if0
          W0 = Atom(ia)%phot(itran)%W0
          if1 = Atom(ia)%phot(itran)%if1
          W1 = Atom(ia)%phot(itran)%W1
          nf = if1 - if0

          ! For each frequency
          do ifreq=if0,if1

            ! If left boundary
            if (ifreq.eq.if0) then

              ! Index
              jjp = iip

              ! Weight
              c1 = W0

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Index
              jjp = iip + nf

              ! Weight
              c1 = W1

            ! Normal
            else

              ! Index
              jjp = iip + ifreq - if0

              ! Weight
              c1 = Wfreq(ifreq)

            end if

            ! Weight with cross section and constants
            c1 = c0*c1*Atom(ia)%phot(itran)%alpha(ifreq)/ &
                 omega(ifreq)

            ! Weight for J00
            WF = c1*Stk(ifreq)

            ! Add the contribution to the integral J00
            J00P(jtran,1) = J00P(jtran,1) + WF
            J00P(jtran,2) = J00P(jtran,2) + WF*exu(ifreq)*c3

            ! Check ALI
            if (nALIp) cycle

            ! Weight for Lambda operator
            WFS = c1*rPhot(jjp)

            ! Add the contribution to the Lambda operator
            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + &
                                 WFS*exu(ifreq)*c3

          end do ! Frequencies

          ! Advance index
          iip = iip + nf + 1

        end do ! b-f transitions
      end do ! atoms

      ! If pre-computed exponential
      if (PIRAM.and.pf1.ge.pf0) then

        ! Nullify local pointer
        nullify(exu)

      ! If locally allocated
      else if (pf1.ge.pf0) then

        ! Free exponential
        deallocate(exu)
        nullify(exu)

      end if ! Pre-computed exponential

      return

      ! Deceive the compiler
      c0 = J00S(1)

      end subroutine JcalcI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add contribution to the integrals of the mean intensity in
      !! MPI for bound-bound transitions\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!      MPID(MPI_class): Structure with MPI data\n
      !!     Wfreq(double(:)): Frequency trapezoidal weights\n
      !!        proc(integer): CPU ID\n
      !!       Stk(double(:)): Intensity\n
      !!     rLine(double(:)): Bound-bound transition strength\n
      !!    Prof(double(:,:)): Bound-bound normalized line profiles\n
      !!    Norm(double(:,:)): Normalization factor for the line
      !!                       profiles\n
      !!    Bstk(double(:,:)): Mean intensity partial integrals\n
      !!    BLam(double(:,:)): Bound-bound Lambda operator partial
      !!                       integral\n
      !!         ALI(logical): If computing Lambda operator
      subroutine FIntI_line(Atom,MPID,Wfreq,proc,Stk,rLine,Prof, &
                            Norm,Bstk,BLam,ALI)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: ALI
      integer, intent(in):: proc
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(in):: rLine
      double precision, dimension(:,:), intent(inout):: Norm
      double precision, dimension(:,:), intent(inout):: BLam
      double precision, dimension(:,:), intent(inout):: BStk
      double precision, dimension(:,:), intent(in):: Prof

      ! Local

      logical:: nALI

      integer:: ifreq,ifreqs,ia,ftran,itran,jftran
      integer:: if0,if1,if0p,if0s,if1s,iil

      double precision:: WF,W0,W1!,WFS


      !
      ! Initializations
      !

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      ! Not calculating ALI
      nALI = .not.ALI


      !
      ! Calculate J00 for transitions (b-b and b-f)
      !

      ! Initialize indexes for profiles
      iil = 0

      ! For each atom
      do ia=1,nA

        ! For each FS transition
        do ftran=1,Atom(ia)%nftran

          ! Transition index
          itran = Atom(ia)%ifst(ftran)

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(itran)%Mabsent(proc)) cycle

          ! Continuous transition index
          jftran = ftran + Atom(ia)%tfshift

          ! Get limits data
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
          Norm(1,jftran) = Norm(1,jftran) + WF

         !! Stimulated emission
         !if (stm) then
         !  WFS = W0*Prof(iil,2)
         !  Norm(2,jftran) = Norm(2,jftran) + WFS
         !end if

          ! Add the contribution to the integrals J00 in this
          ! direction
          BStk(1,jftran) = BStk(1,jftran) + WF*Stk(if0s)

         !! Stimulated emission
         !if(stm) &
         !  BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(if0s)

          ! If calculating lambda operator
          if (ALI) &
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ! Shifted frequency
            ifreqs = ifreq - if0p

            ! Add the profile of the line to the weights and
            ! the norm integral
            iil = iil + 1
            WF = Wfreq(ifreq)*Prof(iil,1)
            Norm(1,jftran) = Norm(1,jftran) + WF

           !! Stimulated emission
           !if (stm) then
           !  WFS = Wfreq(ifreq)*Prof(iil,2)
           !  Norm(2,jftran) = Norm(2,jftran) + WFS
           !end if

            ! Add the contribution to the integrals J00 in this
            ! direction
            BStk(1,jftran) = BStk(1,jftran) + WF*Stk(ifreqs)

           !! Stimulated emission
           !if(stm) BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(ifreqs)

            ! If not computing lambda operator, skip
            if (nALI) cycle

            ! Add contribution to lambda operator
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

          end do ! frequencies

          ! If no more than 1 frequency
          if (if1s.le.if0s) cycle

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W1*Prof(iil,1)
          Norm(1,jftran) = Norm(1,jftran) + WF

         !! Stimulated emission
         !if (stm) then
         !  WFS = W1*Prof(iil,2)
         !  Norm(2,jftran) = Norm(2,jftran) + WFS
         !end if

          ! Add the contribution to the integrals J00 in this
          ! direction
          BStk(1,jftran) = BStk(1,jftran) + WF*Stk(if1s)

         !! Stimulated emission
         !if(stm) BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(if1s)

          ! If computing lambda operator
          if (ALI) &
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

        end do ! b-b transitions
      end do ! atoms

      end subroutine FIntI_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Add contribution to the integrals of the frequency dependent
      !! mean intensity and the mean intensity and lambda operator for
      !! bound-free transitions in MPI\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!         MPID(MPI_class): Structure with MPI data\n
      !!        omega(double(:)): Frequency array\n
      !!        Wfreq(double(:)): Frequency trapezoidal weights\n
      !!            pf0(integer): First frequency index for bound-free
      !!                          transitions\n
      !!            pf1(integer): Last frequency index for bound-free
      !!                          transitions\n
      !!               T(double): Temperature\n
      !!           proc(integer): CPU ID\n
      !!              WA(double): Angular integral weight\n
      !!          Stk(double(:)): Intensity\n
      !!        rPhot(double(:)): Bound-free transition strength\n
      !!       J00P(double(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!         J00C(double(:)): Mean intensity with frequency
      !!                            dependence\n
      !!  LambdaP(double(:,:,:)): Bound-free Lambda operator\n
      !!            ALI(logical): If computing Lambda operator\n
      !!           ALIp(logical): If computing Lambda operator for
      !!                          photoionizations\n
      !!         iexu(double(:)): Pre-computed frequency exponential
      subroutine FIntI_rest(Atom,MPID,omega,Wfreq,pf0,pf1,T,proc,WA, &
                            Stk,rPhot,J00P,J00C,LambdaP,ALI,ALIp,iexu)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: ALI,ALIP
      integer, intent(in):: proc
      integer, intent(in):: pf0,pf1
      double precision, intent(in):: WA
      double precision, intent(in):: T
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(in):: iexu
      double precision, dimension(:), intent(inout):: J00C
      double precision, dimension(:), intent(in):: rPhot
      double precision, dimension(:,:),intent(inout):: J00P
      double precision, dimension(:,:,:),intent(inout):: LambdaP

      ! Local

      logical:: nALI,yALI

      integer:: ifreq,ifreqs,ia,itran,jtran,if0,if1,if0p,if0s,if1s,iip

      double precision:: WF,WFS,c0,c1,W0,W1
      double precision, dimension(:), allocatable:: exu


      !
      ! Initializations
      !

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      ! Not ALI
      yALI = ALI.and.ALIp
      nALI = .not.yALI


      !
      ! Calculate J00C and frequency exponential
      !

      ! For each frequency
      do ifreq=MPID%if0(proc),MPID%if1(proc)

        ! Get shifted frequency
        ifreqs = ifreq - if0p

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreqs)

      end do ! Frequencies

      ! Allocate exu if frequencies
      if (pf1.ge.pf0) allocate(exu(pf0:pf1))

      ! If pre-computed
      if (PIRAM.and.pf1.ge.pf0) then

        ! Copy from input
        exu = iexu

      ! If no pre-computed
      else

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)

          ! Get inverse exponential
          exu(ifreq) = diexp(WF)

        end do ! frequencies

      end if ! pre-computed


      !
      ! Calculate J00 for transitions (b-f)
      !

      ! Variables for photoion.
      c0 = WA*4d-8*pi/convF

      ! Initialize indexes for profiles
      iip = 0

      ! For each atom
      do ia=1,nA

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(itran)%Mabsent(proc)) cycle

          ! Get running transition index
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
          WF = c1*Stk(if0s)

          ! Add the contribution to the integral J00
          J00P(jtran,1) = J00P(jtran,1) + WF
          J00P(jtran,2) = J00P(jtran,2) + WF*exu(if0)

          ! If computing lambda operator
          if (yALI) then

            ! Weight for Lambda operator
            iip = iip + 1
            WFS = c1*rPhot(iip)

            ! Add the contribution to the Lambda operator
            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + WFS*exu(if0)

          end if

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ! Shift index
            ifreqs = ifreq - if0p

            ! Weight with cross section and constants
            c1 = c0*Wfreq(ifreq)*Atom(ia)%phot(itran)%alpha(ifreq)/ &
                 omega(ifreq)

            ! Weight for J00
            WF = c1*Stk(ifreqs)

            ! Add the contribution to the integral J00
            J00P(jtran,1) = J00P(jtran,1) + WF
            J00P(jtran,2) = J00P(jtran,2) + WF*exu(ifreq)

            ! ALI
            if (nALI) cycle

            ! Weight for Lambda operator
            iip = iip + 1
            WFS = c1*rPhot(iip)

            ! Add the contribution to the Lambda operator
            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + WFS*exu(ifreq)

          end do ! frequencies

          if (if1s.le.if0s) cycle

          ! Weight with cross section and constants
          c1 = c0*W1*Atom(ia)%phot(itran)%alpha(if1)/ &
               omega(if1)

          ! Weight for J00
          WF = c1*Stk(if1s)

          ! Add the contribution to the integral J00
          J00P(jtran,1) = J00P(jtran,1) + WF
          J00P(jtran,2) = J00P(jtran,2) + WF*exu(if1)

          ! Add the contribution to the Lambda operator
          if (yALI) then

            ! Weight for Lambda operator
            iip = iip + 1
            WFS = c1*rPhot(iip)

            ! Contribution to lambda operator
            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + WFS*exu(if1)

          end if ! Lambda operator

        end do ! b-f transitions
      end do ! atoms

      ! Free
      if (allocated(exu)) deallocate(exu)

      end subroutine FIntI_rest

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the contribution to the frequency dependent mean
      !! intensity in serial mode\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!          iph(integer): Direction azimuth index\n
      !!          ith(integer): Direction polar index\n
      !!        Stk(double(:)): Intensity\n
      !!       J00C(double(:)): Frequency dependent mean intensity
      subroutine JcalcJ(Geom,iph,ith,Stk,J00C)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      integer, intent(in):: ith,iph
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(inout):: J00C

      ! Local

      integer:: ifreq

      double precision:: WA

      !
      ! Initializations
      !

      ! Angular weight
      WA = Geom%W_mu(ith)*Geom%W_mux(iph)

      !
      ! Calculate J00C and frequency exponential
      !

      ! For each frequency
      do ifreq=1,nfreq

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreq)

      end do ! frequencies

      end subroutine JcalcJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the contribution to the frequency dependent mean
      !! intensity in MPI\n
      !!  MPID(MPI_class): Structure with MPI data\n
      !!       WA(double): Angular weight\n
      !!    proc(integer): CPU ID\n
      !!   Stk(double(:)): Intensity\n
      !!  J00C(double(:)): Frequency dependent mean intensity
      subroutine FIntJ(MPID,WA,proc,Stk,J00C)

      ! I/O

      type(MPI_class), intent(in):: MPID
      integer, intent(in):: proc
      double precision, intent(in):: WA
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(inout):: J00C

      ! Local

      integer:: ifreq,ifreqs,if0p


      !
      ! Initializations
      !

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      !
      ! Calculate J00C
      !

      ! For each frequency
      do ifreq=MPID%if0(proc),MPID%if1(proc)

        ! Shift index
        ifreqs = ifreq - if0p

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreqs)

      end do ! frequencies

      end subroutine FIntJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the K=Q=0 component of the integrated radiation
      !! field tensors in serial mode\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!         Wfreq(double(:)): Frequency trapezoidal weights\n
      !!               WA(double): Angular weight\n
      !!          inpt(double(:)): Quantity to integrate\n
      !!        Prof(double(:,:)): Bound-bound normalized line
      !!                           profiles\n
      !!   JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                           over the absorption profile\n
      !!  JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                           over the emission profile\n
      subroutine Jgen(Atom,Wfreq,WA,inpt,Prof,JKQ,JKQS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      double precision, intent(in):: WA
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(nfreq), intent(in):: inpt
      double precision, dimension(:,:), intent(in):: Prof
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran), &
                       intent(inout):: JKQS

      ! Local

      integer:: ifreq,if0,if1,ia,itran,jtran,iil,jjl,nf

      double precision:: WF,WFS,W0,W1


      !
      ! Calculate J00 for transitions (b-b)
      !

      ! Initialize index for profiles
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
              WF = W0*Prof(iil,1)*WA

              ! Stimulated
              if (stm) WFS = W0*Prof(iil,2)*WA

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

            end if ! Extreme of not

            ! Add the contribution to the JKQ integral
            JKQ(0,0,jtran) = JKQ(0,0,jtran) + &
                             WF*dcmplx(inpt(ifreq),0d0)

            ! Stimulated
            if(stm) &
              JKQS(0,0,jtran) = JKQS(0,0,jtran) + &
                                WFS*dcmplx(inpt(ifreq),0d0)

          end do ! frequencies

          ! Advance index
          iil = iil + nf + 1

        end do ! b-b transitions
      end do ! atoms

      end subroutine Jgen

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the contribution to K=Q=0 component of the
      !! integrated radiation field tensors in MPI\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!          MPID(MPI_class): Structure with MPI data\n
      !!         Wfreq(double(:)): Frequency trapezoidal weights\n
      !!            proc(integer): CPU ID\n
      !!          inpt(double(:)): Quantity to integrate\n
      !!        Prof(double(:,:)): Bound-bound normalized line
      !!                           profiles\n
      !!      Norm(double(:,:,:)): Normalization factor for the
      !!                           line profiles\n
      !!    Bstk(dcomplex(:,:,:)): Radiation field tensors partial
      !!                           integrals\n
      subroutine FJgInt(Atom,MPID,Wfreq,proc,inpt,Prof,Norm,Bstk)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      integer, intent(in):: proc
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(:), intent(in):: inpt
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:), intent(inout):: Norm
      double precision, dimension(:,:), intent(inout):: BStk

      ! Local

      integer:: ifreq,ia,itran,jtran,if0,if1,if0p,iil

      double precision:: WF,WFS,W0,W1


      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1


      !
      ! Calculate J00 for transitions (b-b)
      !

      ! Initialize index for profiles
      iil = 0

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(itran)%Mabsent(proc)) cycle

          jtran = itran + Atom(ia)%tshift

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%Mif0(itran,proc)
          W0 = Atom(ia)%MW0(itran,proc)
          if1 = Atom(ia)%Mif1(itran,proc)
          W1 = Atom(ia)%MW1(itran,proc)

          !
          ! Boundary points
          !

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W0*Prof(iil,1)
          Norm(1,jtran) = Norm(1,jtran) + WF

          ! Simulated
          if (stm) then
            WFS = W0*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! Add the contribution to the JKQ integral
          BStk(1,jtran) = BStk(1,jtran) + WF*inpt(if0)

          ! Stimulated
          if(stm) &
            BStk(2,jtran) = BStk(2,jtran) + WFS*inpt(if0)

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

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

            ! Add the contribution to the JKQ integral
            BStk(1,jtran) = BStk(1,jtran) + WF*inpt(ifreq)

            ! Stimulated
            if(stm) BStk(2,jtran) = BStk(2,jtran) + WFS*inpt(ifreq)

          end do ! frequencies

          ! If less than two frequencies, skip
          if (if1.le.if0) cycle

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

          ! Add the contribution to the JKQ integral
          BStk(1,jtran) = BStk(1,jtran) + WF*inpt(if1)

          ! Stimulated
          if(stm) &
            BStk(2,jtran) = BStk(2,jtran) + WFS*inpt(if1)

        end do ! b-b transitions
      end do ! atoms

      end subroutine FJgInt

!#####################################################################
!#####################################################################
!#####################################################################

      end module jcalci_mod
