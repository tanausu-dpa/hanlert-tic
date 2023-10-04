      !> Mean intensities
      module jcalci_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/26/2017
!  Last version:
!     10/25/2022 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/25/2022:    V3.0.2 - Changed iexu and exu to pointers to
!                             correctly manage the data regardless of
!                             the size and allocation status, as well
!                             as avoiding the copy of allocated data
!                             internally (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     09/30/2021:    V2.0.1 - Added missing nowait (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Split the Fint routines in two
!                             different routines to facilitate the
!                             use of OpenMP (TdPA)
!                           - Added OpenMP to the serial JcalcI,
!                             JcalcJ, and Jgen subroutines (TdPA)
!                           - Changed some parts of the source to
!                             make easier using OpenMP (TdPA)
!
!     07/31/2019:    V1.3.2 - The variable exu is allocated and not
!                             given a fixed size (TdPA)
!
!     12/10/2019:    V1.3.1 - The photoionization RAM storing has its
!                             own flag now, PIRAM (TdPA)
!
!     05/31/2019:    V1.3.0 - Changed the dimensionality of the
!                             profile and ratio variables. Now it
!                             runs sequentially on atoms, transitions
!                             and frequencies to save memory and
!                             reduce the size of data shared through
!                             MPI messages (TdPA)
!
!     05/08/2019:    V1.2.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!
!     04/08/2019:    V1.1.2 - Bugfix: There was an out-of-bounds when
!                             storing profiles, but without explicit
!                             photoionizations (TdPA)
!
!     03/18/2019:    V1.1.1 - Added option of precomputing the
!                             exponentials (TdPA)
!                           - Can skip ALI computations (TdPA)
!
!     02/20/2019:    V1.1.0 - Exponentials are done with diexp (TdPA)
!
!     08/31/2017:    V1.0.8 - Bugfix: The input does not change in
!                             terms of shift in FJgIn (TdPA)
!
!     08/30/2017:    V1.0.7 - Some variables now need shifts because
!                             of changes in the master calls (TdPA)
!
!     08/22/2017:    V1.0.6 - Changed inputs in Fint to comply with
!                             changes in solver (TdPA)
!
!     06/22/2017:    V1.0.5 - Added Jgen and FJgInt (TdPA)
!
!     06/14/2017:    V1.0.4 - exu allocated in a range that is
!                             actually used (TdPA)
!
!     06/13/2017:    V1.0.3 - Limit the exponentials to the region
!                             with photoionizations (TdPA)
!
!     06/12/2017:    V1.0.2 - Take advantage of knowing the true
!                             limits of transitions (TdPA)
!
!     06/09/2017:    V1.0.1 - Added JcalcJ and FintJ (TdPA)
!                           - Removed debugging variable in FintI
!                             (TdPA)
!
!     04/26/2017:    V1.0.0 - First version (TdPA)
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
!  JcalcI
!    Calculate the contribution to J00, J00S, J00P, J00C, and Lambda
!    operators (serial)
!  FIntI_line
!    Calculate partially the integral of J00, J00S and b-b Lambda
!    operator (mpi)
!  FIntI_rest
!    Calculate partially the integral of J00P, J00C, and b-f lambda
!    operator (mpi)
!  JcalcJ
!    Calculate the contribution to J00C (serial)
!  FIntJ
!    Calculate partially the integral of J00C (mpi)
!  Jgen
!    Calculate the contribution to J00 and J00, complex (serial)
!  FJgInt
!    Calculate partially the integral of J00 and J00, complex (mpi)
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use math_mod
      use omp_mod
      use parameters_mod , only : cZero, c2, convF, pi, cSaha, kb, &
                                  fktoJ
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the integrals of the mean intensity
      !! when only one CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!              pf0(integer): First frequency index of this
      !!                            CPU for bound-free transitions\n
      !!              pf1(integer): Last frequency index of this
      !!                            CPU for bound-free transitions\n
      !!                 T(dfloat): Temperature\n
      !!                ne(dfloat): Electron density\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!          rLine(dfloat(:)): Bound-bound transition strength\n
      !!          rPhot(dfloat(:)): Bound-free transition strength\n
      !!         Prof(dfloat(:,:)): Bound-bound normalized line
      !!                            profiles\n
      !!            J00(dfloat(:)): Mean intensity integrated
      !!                            over absorption profile\n
      !!           J00S(dfloat(:)): Mean intensity integrated
      !!                            over emission profile\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!           J00C(dfloat(:)): Mean intensity with frequency
      !!                            dependence\n
      !!      LambdaL(dfloat(:,:)): Bound-bound Lambda operator\n
      !!    LambdaP(dfloat(:,:,:)): Bound-free Lambda operator\n
      !!              ALI(logical): Compute for ALI\n
      !!           iexu(dfloat(:)): Pre-computed exponentials
      subroutine JcalcI(Atom,Geom,omega,Wfreq,pf0,pf1,T,ne,iph,ith, &
                        Stk,rLine,rPhot,Prof,J00,J00S,J00P,J00C, &
                        LambdaL,LambdaP,ALI,iexu)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      logical, intent(in):: ALI
      integer, intent(in):: ith,iph,pf0,pf1
      double precision, intent(in):: T, ne
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), intent(in):: Stk
      double precision, dimension(:), intent(inout):: J00C
      double precision, dimension(:), pointer, intent(in):: iexu
      double precision, dimension(:), intent(inout):: J00
      double precision, dimension(:), intent(inout):: J00S
      double precision, dimension(:), intent(in):: rLine
      double precision, dimension(:), intent(in):: rPhot
      double precision, dimension(:,:),intent(inout):: J00P
      double precision, dimension(:,:), intent(inout):: LambdaL
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:,:),intent(inout):: LambdaP

      ! Local

      logical:: nALI

      integer:: ifreq,if0,if1,ia,itran,jtran,ftran,jftran
      integer:: iil,iip,jjl,jjp,nf

      double precision:: c0,c1,c3,Saha,arg,WA,WF,WFS,W0,W1
      double precision, dimension(:), pointer:: exu


      !
      ! Initializations
      !

      ! Saha non-line-dependent part
      Saha = cSaha*ne/(T**(1.5d0))
      arg = fktoJ/kb/T

      ! Angular weight
      WA = Geom%W_mu(ith)*Geom%W_mux(iph)

      ! Not ALI
      nALI = .not.ALI

      ! Allocate or point exu
      if (PIRAM.and.pf1.ge.pf0) then
        exu(pf0:pf1) => iexu
      else if (pf1.ge.pf0) then
        allocate(exu(pf0:pf1))
      end if

!$omp parallel default(none) &
!$omp private(c0,ifreq,WF,WFS,ia,ftran,itran,jtran,jftran) &
!$omp private(iil,iip,jjl,jjp,nf,if0,if1,c1,c3) &
!$omp private(W0,W1) &
!$omp shared(Saha,arg,WA,nALI,Stk,pf0,pf1,PIRAM,T,iexu,omega) &
!$omp shared(Prof,rLine,rPhot,J00C,exu) &
!$omp shared(Wfreq,Atom,nfreq,na,stm) &
!$omp reduction(+ : J00,LambdaL,J00P,LambdaP)


      !
      ! Calculate J00C and frequency exponential
      !

      ! For each frequency
!$omp do
      do ifreq=1,nfreq

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreq)

      end do
!$omp end do nowait

      ! If no pre-computed
      if (.not.PIRAM.or.pf1.lt.pf0) then

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
!$omp do
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)
          exu(ifreq) = diexp(WF)

        end do ! frequencies
!$omp end do

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

          itran = Atom(ia)%ifst(ftran)
          jftran = ftran + Atom(ia)%tfshift

          ! Get limits data
          if0 = Atom(ia)%if0(itran)
          W0 = Atom(ia)%W0(itran)
          if1 = Atom(ia)%if1(itran)
          W1 = Atom(ia)%W1(itran)
          nf = if1 - if0
!$omp do
          ! For each frequency
          do ifreq=if0,if1

            ! If left boundary
            if (ifreq.eq.if0) then

              ! Index
              jjl = iil

              ! Weight
              WF = W0*Prof(iil,1)*WA
             !if (stm) WFS = W0*Prof(jjl,2)*WA

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Index
              jjl = iil + nf

              ! Weight
              WF = W1*Prof(jjl,1)*WA
             !if (stm) WFS = W1*Prof(jjl,2)*WA

            ! Normal
            else

              ! Index
              jjl = iil + ifreq - if0

              ! Weight
              WF = Wfreq(ifreq)*Prof(jjl,1)*WA
             !if (stm) WFS = Wfreq(ifreq)*Prof(jjl,2)*WA

            end if

            ! Add the contribution to the integrals J00
            J00(jftran) = J00(jftran) + WF*Stk(ifreq)
           !if(stm) J00S(jftran) = J00S(jftran) + WFS*Stk(ifreq)

            ! No ALI in iteration
            if (nALI) cycle

            ! Add contribution to Lambda operator
            LambdaL(1,jftran) = LambdaL(1,jftran) + WF*rLine(jjl)

          end do
!$omp end do nowait

          ! Advance iil
          iil = iil + nf + 1

        end do ! b-b transitions

        ! For each b-f transition
        do itran=1,Atom(ia)%nphot

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
!$omp do
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
            if (nALI) cycle

            ! Weight for Lambda operator
            WFS = c1*rPhot(jjp)

            ! Add the contribution to the Lambda operator
            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + &
                                 WFS*exu(ifreq)*c3
          end do
!$omp end do nowait

          ! Advance index
          iip = iip + nf + 1

        end do ! b-f transitions
      end do ! atoms

!$omp end parallel

      ! Free exu
      if (PIRAM.and.pf1.ge.pf0) then
        nullify(exu)
      else if (pf1.ge.pf0) then
        deallocate(exu)
        nullify(exu)
      end if

      end subroutine JcalcI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the integrals of the mean intensity
      !! with several CPU for b-b transitions.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!             proc(integer): CPU ID\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!          rLine(dfloat(:)): Bound-bound transition strength\n
      !!         Prof(dfloat(:,:)): Bound-bound non-normalized line
      !!                            profiles\n
      !!         Norm(dfloat(:,:)): Normalization factor for the
      !!                            line profiles\n
      !!         Bstk(dfloat(:,:)): Mean intensity partial integrals\n
      !!         BLam(dfloat(:,:)): Bound-bound Lambda operator
      !!                            partial integral\n
      !!              ALI(logical): Compute for ALI
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

      double precision:: WF,W0,W1
     !double precision:: WFS


      !
      ! Initializations
      !

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      ! Not ALI
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

          itran = Atom(ia)%ifst(ftran)

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(itran)%Mabsent(proc)) cycle

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
         !if (stm) then
         !  WFS = W0*Prof(iil,2)
         !  Norm(2,jftran) = Norm(2,jftran) + WFS
         !end if

          ! Add the contribution to the integrals J00 in this
          ! direction
          BStk(1,jftran) = BStk(1,jftran) + WF*Stk(if0s)
         !if(stm) &
         !  BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(if0s)

          if (ALI) &
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ifreqs = ifreq - if0p

            ! Add the profile of the line to the weights and
            ! the norm integral
            iil = iil + 1
            WF = Wfreq(ifreq)*Prof(iil,1)
            Norm(1,jftran) = Norm(1,jftran) + WF
           !if (stm) then
           !  WFS = Wfreq(ifreq)*Prof(iil,2)
           !  Norm(2,jftran) = Norm(2,jftran) + WFS
           !end if

            ! Add the contribution to the integrals J00 in this
            ! direction
            BStk(1,jftran) = BStk(1,jftran) + WF*Stk(ifreqs)
           !if(stm) BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(ifreqs)

            if (nALI) cycle
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

          end do ! frequencies

          if (if1s.le.if0s) cycle

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W1*Prof(iil,1)
          Norm(1,jftran) = Norm(1,jftran) + WF
         !if (stm) then
         !  WFS = W1*Prof(iil,2)
         !  Norm(2,jftran) = Norm(2,jftran) + WFS
         !end if

          ! Add the contribution to the integrals J00 in this
          ! direction
          BStk(1,jftran) = BStk(1,jftran) + WF*Stk(if1s)
         !if(stm) BStk(2,jftran) = BStk(2,jftran) + WFS*Stk(if1s)

          if (ALI) &
            BLam(1,jftran) = BLam(1,jftran) + WF*rLine(iil)

        end do ! b-b transitions
      end do ! atoms

      end subroutine FIntI_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the integrals of the mean intensity
      !! with several CPU for the photoionizations and continuum.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!              pf0(integer): First frequency index of this
      !!                            CPU for bound-free transitions\n
      !!              pf1(integer): Last frequency index of this
      !!                            CPU for bound-free transitions\n
      !!                 T(dfloat): Temperature\n
      !!             proc(integer): CPU ID\n
      !!                WA(dfloat): Angular integral weight\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!          rPhot(dfloat(:)): Bound-free transition strength\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!           J00C(dfloat(:)): Mean intensity with frequency
      !!                            dependence\n
      !!    LambdaP(dfloat(:,:,:)): Bound-free Lambda operator\n
      !!              ALI(logical): Compute for ALI\n
      !!           iexu(dfloat(:)): Pre-computed exponentials
      subroutine FIntI_rest(Atom,MPID,omega,Wfreq,pf0,pf1,T,proc,WA, &
                            Stk,rPhot,J00P,J00C,LambdaP,ALI,iexu)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: ALI
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

      logical:: nALI

      integer:: ifreq,ifreqs,ia,itran,jtran
      integer:: if0,if1,if0p,if0s,if1s,iip

      double precision:: WF,WFS,c0,c1,W0,W1
      double precision, dimension(:), allocatable:: exu


      !
      ! Initializations
      !

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      ! Not ALI
      nALI = .not.ALI


      !
      ! Calculate J00C and frequency exponential
      !

      ! For each frequency
      do ifreq=MPID%if0(proc),MPID%if1(proc)

        ifreqs = ifreq - if0p

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreqs)

      end do

      ! Allocate exu
      if (pf1.ge.pf0) allocate(exu(pf0:pf1))

      ! If pre-computed
      if (PIRAM.and.pf1.ge.pf0) then

        exu = iexu

      ! If no pre-computed
      else

        ! Initialize exponential argument constant
        c0 = c2*1d4/T

        ! For each frequency with photoionization
        do ifreq=pf0,pf1

          ! Argument frequency exponential
          WF = c0*omega(ifreq)
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

          if (ALI) then

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
          if (ALI) then

            ! Weight for Lambda operator
            iip = iip + 1
            WFS = c1*rPhot(iip)

            LambdaP(1,jtran,1) = LambdaP(1,jtran,1) + WFS
            LambdaP(1,jtran,2) = LambdaP(1,jtran,2) + WFS*exu(if1)

          end if

        end do ! b-f transitions

      end do ! atoms

      end subroutine FIntI_rest

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the directional integral of the mean
      !! intensity when only one CPU.\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!           J00C(dfloat(:)): Mean intensity with frequency
      !!                            dependence
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
!$omp parallel do default(none) &
!$omp private(ifreq) shared(nfreq,J00C,WA,Stk)
      do ifreq=1,nfreq

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreq)

      end do ! frequencies
!$omp end parallel do

      end subroutine JcalcJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the directional integral of the mean
      !! intensity with several CPU.\n
      !!     MPID(MPI_class): Structure with MPI data\n
      !!          WA(dfloat): Angular weight\n
      !!       proc(integer): CPU ID\n
      !!    Stk(dfloat(:,:)): Stokes parameters\n
      !!     J00C(dfloat(:)): Mean intensity with frequency dependence
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

        ifreqs = ifreq - if0p

        ! Calculate continuum mean intensity
        J00C(ifreq) = J00C(ifreq) + WA*Stk(ifreqs)

      end do ! frequencies

      end subroutine FIntJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds intensity contribution to the integrals of the radiation
      !! field tensors when only one CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!                WA(dfloat): Angular weight\n
      !!           inpt(dfloat(:)): Quantity to integrate\n
      !!         Prof(dfloat(:,:)): Bound-bound normalized line
      !!                            profiles\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQS(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      subroutine Jgen(Atom,Wfreq,WA,inpt,Prof,JKQ,JKQS)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      double precision, intent(in):: WA
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(nfreq), intent(in):: inpt
      double precision, dimension(:,:), intent(in):: Prof
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout)::JKQS

      ! Local

      integer:: ifreq,if0,if1,ia,itran,jtran,iil,jjl,nf

      double precision:: WF,WFS,W0,W1


      !
      ! Calculate J00 for transitions (b-b)
      !

!$omp parallel default(none) &
!$omp private(iil,jjl,ia,itran,jtran,if0,if1,W0,W1,nf,ifreq,WF,WFS) &
!$omp shared(nA,Atom,WA,Prof,inpt,Wfreq,stm) &
!$omp reduction(+ : JKQ,JKQS)

      ! Initialize index for profiles
      iil = 1

      ! For each atom
      do ia=1,nA

        ! For each b-b transition
        do itran=1,Atom(ia)%ntran

          jtran = itran + Atom(ia)%tshift

          ! Limits in frequency index of this transition
          if0 = Atom(ia)%if0(itran)
          W0 = Atom(ia)%W0(itran)
          if1 = Atom(ia)%if1(itran)
          W1 = Atom(ia)%W1(itran)
          nf = if1 - if0

!$omp do
          ! For each frequency
          do ifreq=if0,if1

            ! If left boundary
            if (ifreq.eq.if0) then

              ! Index
              jjl = iil

              ! Weight
              WF = W0*Prof(iil,1)*WA
              if (stm) WFS = W0*Prof(iil,2)*WA

            ! If right boundary
            else if (ifreq.eq.if1) then

              ! Index
              jjl = iil + nf

              ! Weight
              WF = W1*Prof(jjl,1)*WA
              if (stm) WFS = W1*Prof(jjl,2)*WA

            ! Normal
            else

                ! Index
                jjl = iil + ifreq - if0

                ! Weight
                WF = Wfreq(ifreq)*Prof(jjl,1)*WA
                if (stm) WFS = Wfreq(ifreq)*Prof(jjl,2)*WA

            end if

            ! Add the contribution to the JKQ integral
            JKQ(0,0,jtran) = JKQ(0,0,jtran) + &
                             WF*dcmplx(inpt(ifreq),0d0)
            if(stm) &
              JKQS(0,0,jtran) = JKQS(0,0,jtran) + &
                                WFS*dcmplx(inpt(ifreq),0d0)

          end do ! frequencies
!$omp end do nowait

          ! Advance index
          iil = iil + nf + 1

        end do ! b-b transitions
      end do ! atoms
!$omp end parallel

      end subroutine Jgen

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds intensity contribution to the integrals of the radiation
      !! field tensors with several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!             proc(integer): CPU ID\n
      !!           inpt(dfloat(:)): Quantity to integrate\n
      !!         Prof(dfloat(:,:)): Bound-bound non-normalized line
      !!                            profiles\n
      !!       Norm(dfloat(:,:,:)): Normalization factor for the
      !!                            line profiles\n
      !!     Bstk(dcomplex(:,:,:)): Radiation field tensors partial
      !!                            integrals\n
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
          if (stm) then
            WFS = W0*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! Add the contribution to the JKQ integral
          BStk(1,jtran) = BStk(1,jtran) + WF*inpt(if0)
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
            if (stm) then
              WFS = Wfreq(ifreq)*Prof(iil,2)
              Norm(2,jtran) = Norm(2,jtran) + WFS
            end if

            ! Add the contribution to the JKQ integral
            BStk(1,jtran) = BStk(1,jtran) + WF*inpt(ifreq)
            if(stm) BStk(2,jtran) = BStk(2,jtran) + WFS*inpt(ifreq)

          end do ! frequencies

          if (if1.le.if0) cycle

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W1*Prof(iil,1)
          Norm(1,jtran) = Norm(1,jtran) + WF
          if (stm) then
            WFS = W1*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! Add the contribution to the JKQ integral
          BStk(1,jtran) = BStk(1,jtran) + WF*inpt(if1)
          if(stm) &
            BStk(2,jtran) = BStk(2,jtran) + WFS*inpt(if1)

        end do ! b-b transitions
      end do ! atoms

      end subroutine FJgInt

!#####################################################################
!#####################################################################
!#####################################################################

      end module jcalci_mod
