      !> Radiation field tensors
      module jcalc_mod
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
!     09/29/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/29/2023:    V3.0.5 - Updated to term- and transition-wise
!                             K cut limits (TdPA)
!                           - Avoid computing tensor components
!                             with negative Q values (TdPA)
!
!     11/10/2022:    V3.0.4 - Added the possibility of forcing the
!                             axial asymmetry to be fully given by
!                             the input, ignoring any other
!                             contribution to it (TdPA)
!
!     10/25/2022:    V3.0.3 - Changed iexu and exu to pointers to
!                             correctly manage the data regardless of
!                             the size and allocation status, as well
!                             as avoiding the copy of allocated data
!                             internally (TdPA)
!                           - Implemented restriction of the height
!                             axis range (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!
!     07/08/2022:    V3.0.1 - Removed debugging print (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     09/30/2021:    V2.0.1 - Missing nowaits in jcalc (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Split the Fint routines in two
!                             different routines to facilitate the
!                             use of OpenMP (TdPA)
!                           - Added OpenMP to the serial Jcalc
!                             subroutines (TdPA)
!                           - Changed some parts of the source to
!                             make easier using OpenMP (TdPA)
!
!     01/13/2021:    V1.3.4 - Added addJKQasym subroutine (TdPA)
!
!     07/31/2020:    V1.3.3 - The variable exu is allocated and not
!                             given a fixed size (TdPA)
!
!     12/10/2019:    V1.3.2 - The photoionization RAM storing has its
!                             own flag now, PIRAM (TdPA)
!
!     08/19/2019:    V1.3.1 - Ignoring J00S integrals, commented the
!                             lines (TdPA)
!
!     05/31/2019:    V1.3.0 - Changed the dimensionality of the
!                             profile variable. Now it runs
!                             sequentially on atoms, transitions and
!                             frequencies to save memory and reduce
!                             the size of data shared through MPI
!                             messages (TdPA)
!
!     05/08/2019:    V1.2.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!
!     03/18/2019:    V1.1.1 - Added option of precomputing the
!                             exponentials (TdPA)
!
!     02/20/2019:    V1.1.0 - Exponentials are done with diexp (TdPA)
!
!     08/03/2018:    V1.0.8 - Kcut is now Krad for JKQ (TdPA)
!
!     07/31/2018:    V1.0.7 - Kcut also limits the radiation tensors
!                             that are actually computed (TdPA)
!
!     08/30/2017:    V1.0.6 - Some variables now need shifts because
!                             of changes in the master calls (TdPA)
!
!     08/21/2017:    V1.0.5 - Changed inputs in Fint to comply with
!                             changes in solver (TdPA)
!
!     06/16/2017:    V1.0.4 - Limit the integr to the region
!                             with lines (TdPA)
!
!     06/14/2017:    V1.0.3 - exu allocated in a range that is
!                             actually used (TdPA)
!
!     06/13/2017:    V1.0.2 - Limit the exponentials to the region
!                             with photoionizations (TdPA)
!
!     06/12/2017:    V1.0.1 - Take advantage of knowing the true
!                             limits of transitions (TdPA)
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
!  Jcalc
!    Calculate the contribution to JKQ, JKQS, JKQC, and J00P (serial)
!  FInt_line
!    Calculate partially the integral of JKQ and JKQS
!  FInt_rest
!    Calculate partially the integral of J00P and JKQC (mpi)
!  addJKQasym
!    Adds ad-hoc values to the JKQ tensors
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

      !> Adds contribution to the integrals of the radiation field
      !! tensors when only one CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!              lf0(integer): First frequency index of this
      !!                            CPU for bound-bound transitions\n
      !!              lf1(integer): Last frequency index of this
      !!                            CPU for bound-bound transitions\n
      !!              pf0(integer): First frequency index of this
      !!                            CPU for bound-free transitions\n
      !!              pf1(integer): Last frequency index of this
      !!                            CPU for bound-free transitions\n
      !!                 T(dfloat): Temperature\n
      !!                ne(dfloat): Electron density\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!               iz(integer): Height index\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!         Prof(dfloat(:,:)): Bound-bound normalized line
      !!                            profiles\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQS(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!           iexu(dfloat(:)): Pre-computed exponentials
      subroutine Jcalc(Atom,Geom,omega,Wfreq,lf0,lf1,pf0,pf1,T,ne, &
                       iph,ith,iz,Stk,Prof,JKQ,JKQS,J00P,JKQC, &
                       iexu)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      integer, intent(in):: ith,iph,iz,lf0,lf1,pf0,pf1
      double precision, intent(in):: T, ne
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), pointer, intent(in):: iexu
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:),intent(inout):: J00P
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(inout)::JKQS
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(inout):: JKQC

      ! Local

      integer:: ifreq,if0,if1,ia,itran,jtran,K,iQ,iil,jjl,nf
      integer:: Kmax,lKmax

      double precision:: c0,c1,c3,Saha,arg
      double precision:: WA,WF,WFS,W0,W1
      double precision, dimension(:), pointer:: exu

      complex(kind=8), dimension(-2:2,0:2,lf0:lf1):: integr


      !
      ! Initializations
      !

      Kmax = min(Krad, 2)
      lKmax = min(Kradl, 2)

      ! Saha non-line-dependent part
      Saha = cSaha*ne/(T**(1.5d0))
      arg = fktoJ/kb/T

      ! Angular weight
      WA = Geom%W_mu(ith)*Geom%W_mux(iph)


      !
      ! Calculate JKQC, the argument of JKQ, and the frequency
      ! exponential
      !

      ! Allocate or point exu
      if (PIRAM.and.pf1.ge.pf0) then
        exu(pf0:pf1) => iexu
      else if (pf1.ge.pf0) then
        allocate(exu(pf0:pf1))
      end if

!$omp parallel default(none) &
!$omp private(c0,ifreq,WF,K,iQ,ia,iil,itran,jtran,if0,if1,W0,W1,nf) &
!$omp private(c3,c1,WFS,jjl) &
!$omp shared(Kmax,lKmax,Saha,arg,WA,pf1,pf0,exu,iexu,T,PIRAM,nfreq) &
!$omp shared(JKQC,Stk,Geom,lf0,lf1,axial,integr,Atom,Wfreq,iz,ith) &
!$omp shared(iph,omega,na,Prof,stm) &
!$omp reduction(+ : JKQ,JKQS,J00P)

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

      ! For each frequency
!$omp do
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
                                      Geom%TS(:,iQ,K,iph,ith))

          end do ! Q
        end do ! K
      end do ! frequencies
!$omp end do

      ! For each frequency with lines
!$omp do
      do ifreq=lf0,lf1

        ! For each K
        do K=0,lKmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if (iQ.ne.0.and.axial) cycle

            ! Compute the sum TKQ*Stokes
            integr(iQ,K,ifreq) = sum(Stk(:,ifreq)* &
                                     Geom%TB(:,iQ,K,iph,ith,iz))

          end do ! Q
        end do ! K
      end do ! frequencies
!$omp end do


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
              WF = W0*Prof(jjl,1)*WA
              if (stm) WFS = W0*Prof(jjl,2)*WA

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

            ! For each K
            do K=0,Atom(ia)%Krad(itran)

              ! For each Q
              do iQ=0,K

                ! If axial atmosphere and Q!=0, skip
                if (iQ.ne.0.and.axial) cycle

                ! Add the contribution to the JKQ integral
                JKQ(iQ,K,jtran) = JKQ(iQ,K,jtran) + &
                                  WF*integr(iQ,K,ifreq)
                if(stm) &
                JKQS(iQ,K,jtran) = JKQS(iQ,K,jtran) + &
                                   WFS*integr(iQ,K,ifreq)

              end do ! Q
            end do ! K
          end do ! frequencies
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
!$omp end do nowait
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

      end subroutine Jcalc

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution to the integrals of the radiation field
      !! tensors with several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!              lf0(integer): First frequency index of this
      !!                            CPU for bound-bound transitions\n
      !!              lf1(integer): Last frequency index of this
      !!                            CPU for bound-bound transitions\n
      !!             proc(integer): CPU ID\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!               iz(integer): Height index\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!         Prof(dfloat(:,:)): Bound-bound non-normalized line
      !!                            profiles\n
      !!         Norm(dfloat(:,:)): Normalization factor for the
      !!                            line profiles\n
      !!   Bstk(dcomplex(:,:,:,:)): Radiation field tensors partial
      !!                            integrals\n
      subroutine FInt_line(Atom,MPID,Geom,Wfreq,lf0,lf1,proc, &
                           iph,ith,iz,Stk,Prof,Norm,Bstk)

      ! I/O
      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      integer, intent(in):: proc,lf0,lf1,iph,ith,iz
      double precision, dimension(:), intent(in):: Wfreq
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:), intent(in):: Prof
      double precision, dimension(:,:), intent(inout):: Norm
      complex(kind=8), dimension(0:2,0:2,2,nxtran), &
                                                  intent(inout):: BStk

      ! Local

      integer:: ifreq,ifreqs,ia,itran,jtran,K,iQ,lKmax
      integer:: if0,if1,if0p,if0s,if1s,iil

      double precision:: WF,WFS,W0,W1

      complex(kind=8), dimension(0:2,0:2,lf0:lf1):: integr


      !
      ! Initializations
      !

      lKmax = min(Kradl, 2)

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1

      ! For each frequency with lines
      do ifreq=lf0,lf1

        ifreqs = ifreq - if0p

        ! For each K
        do K=0,lKmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if(iQ.ne.0.and.axial)cycle

            ! Compute the sum TKQ*Stokes
            integr(iQ,K,ifreq) = sum(Stk(:,ifreqs)* &
                                     Geom%TB(:,iQ,K,iph,ith,iz))

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
          if (stm) then
            WFS = W0*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! For each K
          do K=0,Atom(ia)%Krad(jtran)

            ! For each Q
            do iQ=0,K

              ! If axial atmosphere and Q!=0, skip
              if(iQ.ne.0.and.axial)cycle

              ! Add the contribution to the JKQ integral
              BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                   WF*integr(iQ,K,if0)
              if(stm) &
              BStk(iQ,K,2,jtran) = BStk(iQ,K,2,jtran) + &
                                   WFS*integr(iQ,K,if0)

            end do ! Q
          end do ! K

          !
          ! For each non-boundary frequency
          !
          do ifreq=if0+1,if1-1

            ifreqs = ifreq - if0p

            ! Add the profile of the line to the weights and
            ! the norm integral
            iil = iil + 1
            WF = Wfreq(ifreq)*Prof(iil,1)
            Norm(1,jtran) = Norm(1,jtran) + WF
            if (stm) then
              WFS = Wfreq(ifreq)*Prof(iil,2)
              Norm(2,jtran) = Norm(2,jtran) + WFS
            end if

            ! For each K
            do K=0,Atom(ia)%Krad(jtran)

              ! For each Q
              do iQ=0,K

                ! If axial atmosphere and Q!=0, skip
                if(iQ.ne.0.and.axial)cycle

                ! Add the contribution to the JKQ integral
                BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                     WF*integr(iQ,K,ifreq)
                if(stm) BStk(iQ,K,2,jtran) = BStk(iQ,K,2,jtran) + &
                                             WFS*integr(iQ,K,ifreq)

              end do ! Q
            end do ! K
          end do ! frequencies

          if (if1s.le.if0s) cycle

          ! Add the profile of the line to the weights and
          ! the norm integral
          iil = iil + 1
          WF = W1*Prof(iil,1)
          Norm(1,jtran) = Norm(1,jtran) + WF
          if (stm) then
            WFS = W1*Prof(iil,2)
            Norm(2,jtran) = Norm(2,jtran) + WFS
          end if

          ! For each K
          do K=0,Atom(ia)%Krad(jtran)

            ! For each Q
            do iQ=0,K

              ! If axial atmosphere and Q!=0, skip
              if(iQ.ne.0.and.axial)cycle

              ! Add the contribution to the JKQ integral
              BStk(iQ,K,1,jtran) = BStk(iQ,K,1,jtran) + &
                                   WF*integr(iQ,K,if1)
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

      !> Adds contribution to the integrals of the radiation field
      !! tensors with several CPU.\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!          omega(dfloat(:)): Frequency array\n
      !!          Wfreq(dfloat(:)): Frequency trapezoidal weights\n
      !!              pf0(integer): First frequency index of this
      !!                            CPU for bound-free transitions\n
      !!              pf1(integer): Last frequency index of this
      !!                            CPU for bound-free transitions\n
      !!                 T(dfloat): Temperature\n
      !!             proc(integer): CPU ID\n
      !!              iph(integer): Output direction azimuth index\n
      !!              ith(integer): Output direction polar index\n
      !!          Stk(dfloat(:,:)): Stokes parameters\n
      !!         J00P(dfloat(:,:)): Intensity integrals in the
      !!                            photoionization rates\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!           iexu(dfloat(:)): Pre-computed exponentials
      subroutine FInt_rest(Atom,MPID,Geom,omega,Wfreq,pf0,pf1,T, &
                           proc,iph,ith,WA,Stk,J00P,JKQC,iexu)

      ! I/O
      type(Atom_class), dimension(:), intent(in):: Atom
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      integer, intent(in):: proc,iph,ith
      integer, intent(in):: pf0,pf1
      double precision, intent(in):: T,WA
      double precision, dimension(:), intent(in):: Wfreq, omega
      double precision, dimension(:), intent(in):: iexu
      double precision, dimension(0:3,nfreq), intent(in):: Stk
      double precision, dimension(:,:),intent(inout):: J00P
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(inout):: JKQC

      ! Local

      integer:: ifreq,ifreqs,ia,itran,jtran,K,iQ,Kmax
      integer:: if0,if1,if0p,if0s,if1s

      double precision:: WF,c0,c1,W0,W1
      double precision, dimension(:), allocatable:: exu


      !
      ! Initializations
      !

      Kmax = min(Krad, 2)

      ! Lower limit for processor
      if0p = MPID%if0(proc) - 1


      !
      ! Calculate JKQC and frequency exponential
      !

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

      ! For each frequency
      do ifreq=MPID%if0(proc),MPID%if1(proc)

        ifreqs = ifreq - if0p

        ! For each K
        do K=0,Kmax

          ! For each Q
          do iQ=0,K

            ! If axial atmosphere and Q!=0, skip
            if(iQ.ne.0.and.axial)cycle

            ! Integrate JKQC
            JKQC(iQ,K,ifreq) = JKQC(iQ,K,ifreq) + &
                               WA*sum(Stk(:,ifreqs)* &
                                      Geom%TS(:,iQ,K,iph,ith))

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

      end subroutine FInt_rest

!#####################################################################
!#####################################################################
!#####################################################################

      !> Adds contribution of asymmetric inputs to the JKQ\n
      !!   Bfield(Bfield_class): Structure with magnetic field data\n
      !!     Flgsg(Fctsg_class): Structure with factorials and signs\n
      !!   JKQa(dcomplex(:,:,:)): Radiation field tensors extra
      !!                          asymmetries\n
      !!  JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                          over absorption profile\n
      !! JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                          over emission profile\n
      !! JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                          frequency dependence
      subroutine addJKQasym(Bfield,Flgsg,JKQa,JKQ,JKQS,JKQC)

      ! I/O

      type(Bfield_class):: Bfield
      type(Fctsg_class):: Flgsg
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
            JKQC(:,2:3,ifreq,iz) = &
                               JKQa(:,:,iz)*dble(JKQC(3,1,ifreq,iz))
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
                               JKQa(:,:,iz)*dble(JKQC(3,1,ifreq,iz))

          end do ! frequencies
        end do ! heights

      end if

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
            JKQ(:,2:3,ii,iz) = JKQa(:,:,iz)*dble(JKQ(3,1,ii,iz))
          end do

        ! Additive instead
        else

          ! For each transition
          do ii=1,nxtran
            JKQ(:,2:3,ii,iz) = JKQ(:,2:3,ii,iz) + &
                               JKQa(:,:,iz)*dble(JKQ(3,1,ii,iz))
          end do

        end if

        ! If magnetic field
        if (Bfield%Bstrength(jz).gt.TINYB) then

          ! Rotate back to magnetic field reference frame
          call fieldB(JKQ(:,:,:,iz),nxtran,Flgsg, &
                      Bfield%Btheta(jz),Bfield%Bphi(jz),1)

        end if

        ! If stimulated emission
        if (stm) then

          ! If magnetic field
          if (Bfield%Bstrength(jz).gt.TINYB) then

            ! Rotate to the vertical
            call fieldB(JKQS(:,:,:,iz),nxtran,Flgsg, &
                        -Bfield%Btheta(jz), &
                        -Bfield%Bphi(jz),-1)
          end if

          ! Forcing asymmetry
          if (force_asym) then

            ! For each transition
            do ii=1,nxtran
              JKQS(:,2:3,ii,iz) = JKQS(:,2:3,ii,iz) + &
                                  JKQa(:,:,iz)*dble(JKQS(3,1,ii,iz))
            end do

          ! Additive instead
          else

            ! For each transition
            do ii=1,nxtran
              JKQS(:,2:3,ii,iz) = JKQS(:,2:3,ii,iz) + &
                                  JKQa(:,:,iz)*dble(JKQS(3,1,ii,iz))
            end do

          end if

          ! If magnetic field
          if (Bfield%Bstrength(jz).gt.TINYB) then

            ! Rotate back to magnetic field reference frame
            call fieldB(JKQS(:,:,:,iz),nxtran,Flgsg, &
                        Bfield%Btheta(iz),Bfield%Bphi(iz),1)

          end if
        end if ! stimulated emission

      end do ! Heights

      end subroutine addJKQasym

!#####################################################################
!#####################################################################
!#####################################################################

      end module jcalc_mod
