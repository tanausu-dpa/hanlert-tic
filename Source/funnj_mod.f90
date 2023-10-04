      !> n-J symbols
      module funnj_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     06/29/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added an extra check to ensure that
!                             we are not writing in the memoization
!                             if doing OpenMP actively (TdPA)
!
!     12/11/2019:    V1.1.0 - Now uses memoization (TdPA)
!
!     12/11/2019:    V1.1.0 - Now uses memoization (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  fn1:
!    Auxiliar for nJ-symbols
!
!  fn2:
!    Auxiliar for nJ-symbols
!
!  fun3j:
!    3J-symbol
!
!  fun6j:
!    6J-symbol
!
!  fun9j:
!    9J-symbol
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use memoization_mod
      use types_mod

      ! Functions for n-j symbols

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Used in 3j calculation\n
      !!           j1(dfloat): First value\n
      !!           j2(dfloat): Second value\n
      !!           j3(dfloat): Third value\n
      !!           m1(dfloat): Fourth value\n
      !!           m2(dfloat): Fifth value\n
      !!           m3(dfloat): Sixth value\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      double precision function fn1(j1,j2,j3,m1,m2,m3,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      double precision:: j1,j2,j3
      double precision:: m1,m2,m3

      ! Local

      integer:: l1,l2,l3,l4,l5,l6,l7,l8,l9,l10

      l1 = nint(j1+j2-j3)
      l2 = nint(j2+j3-j1)
      l3 = nint(j3+j1-j2)
      l4 = nint(j1+j2+j3)+1
      l5 = nint(j1+m1)
      l6 = nint(j1-m1)
      l7 = nint(j2+m2)
      l8 = nint(j2-m2)
      l9 = nint(j3+m3)
      l10 = nint(j3-m3)

      fn1 = .5d0*(Flgsg%flg(l1) + Flgsg%flg(l2) + Flgsg%flg(l3) - &
                  Flgsg%flg(l4) + Flgsg%flg(l5) + Flgsg%flg(l6) + &
                  Flgsg%flg(l7) + Flgsg%flg(l8) + Flgsg%flg(l9) + &
                  Flgsg%flg(l10))

      return

      end function fn1

!#####################################################################
!#####################################################################
!#####################################################################

      !> Used in 6j calculation\n
      !!          ij1(dfloat): First value\n
      !!          ij2(dfloat): Second value\n
      !!          ij3(dfloat): Third value\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      double precision function fn2(ij1,ij2,ij3,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      integer:: ij1,ij2,ij3

      ! Local

      integer:: l1,l2,l3,l4

      l1 = (ij1+ij2-ij3)/2
      l2 = (ij2+ij3-ij1)/2
      l3 = (ij3+ij1-ij2)/2
      l4 = (ij1+ij2+ij3)/2+1

      fn2 = .5d0*(Flgsg%flg(l1) + Flgsg%flg(l2) + Flgsg%flg(l3) - &
                  Flgsg%flg(l4))

      return

      end function fn2

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculates 3j value, each ji/mi is real\n
      !!           j1(dfloat): First value in 3J symbol\n
      !!           j2(dfloat): Second value in 3J symbol\n
      !!           j3(dfloat): Third value in 3J symbol\n
      !!           m1(dfloat): Fourth value in 3J symbol\n
      !!           m2(dfloat): Fifth value in 3J symbol\n
      !!           m3(dfloat): Sixth value in 3J symbol\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      double precision function fun3j(j1,j2,j3,m1,m2,m3,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      double precision:: j1,j2,j3
      double precision:: m1,m2,m3

      ! Local

      integer:: l1,l2,l3,l4,l5,l6
      integer:: ij1,ij2,ij3
      integer:: im1,im2,im3
      integer:: kmin,kmin1,kmin2
      integer:: kmax,kmax1,kmax2,kmax3
      integer:: i

      double precision:: term1,term2
      double precision:: sgn

      double precision, pointer:: val

      ! If doing memoization
      if (Flgsg%memo) then

        ! Convert to integers
        l1 = nint(2d0*j1)
        l2 = nint(2d0*j2)
        l3 = nint(2d0*j3)
        l4 = nint(2d0*m1)
        l5 = nint(2d0*m2)
        l6 = nint(2d0*m3)

        ! Look for existing element
        val => elem6D(l1,l2,l3,l4,l5,l6,Flgsg%J3)

        ! If found, return value
        if (associated(val)) then
          fun3j = val
          return
        end if ! Found value in jagged array
      end if ! Memoization

      fun3j = 0d0

      ij1 = nint(j1+j1)
      ij2 = nint(j2+j2)
      ij3 = nint(j3+j3)

      if(ij1+ij2-ij3.lt.0) return
      if(ij2+ij3-ij1.lt.0) return
      if(ij3+ij1-ij2.lt.0) return

      im1 = nint(m1+m1)
      im2 = nint(m2+m2)
      im3 = nint(m3+m3)

      if(im1+im2+im3.ne.0) return
      if(abs(im1).gt.ij1) return
      if(abs(im2).gt.ij2) return
      if(abs(im3).gt.ij3) return

      kmin = (ij3-ij1-im2)/2
      kmin1 = kmin
      kmin2 = (ij3-ij2+im1)/2

      kmin = max(-min(kmin,kmin2),0)

      kmax = nint(j1+j2-j3)
      kmax1 = kmax
      kmax2 = nint(j1-m1)
      kmax3 = nint(j2+m2)

      kmax = min(kmax,kmax2,kmax3)

      if(kmin.le.kmax) then

        term1 = fn1(j1,j2,j3,m1,m2,m3,Flgsg)

        sgn = Flgsg%sg((ij1-ij2-im3)/2)

        do i=kmin,kmax
          term2 = Flgsg%flg(i) + Flgsg%flg(kmin1+i) + &
                  Flgsg%flg(kmin2+i) + Flgsg%flg(kmax1-i) + &
                  Flgsg%flg(kmax2-i) + Flgsg%flg(kmax3-i)
          fun3j = Flgsg%sg(i)*exp(term1-term2) + fun3j
        end do

        fun3j = sgn*fun3j

      end if

      ! If memoization, insert
#ifdef _OPENMP
      if (Flgsg%memo.and.Flgsg%can_write) &
        call insert6D(l1,l2,l3,l4,l5,l6,fun3j,Flgsg%J3)
#else
      if (Flgsg%memo) call insert6D(l1,l2,l3,l4,l5,l6,fun3j,Flgsg%J3)
#endif

      return

      end function fun3j

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculates 6j value, each jii is real\n
      !!          j11(dfloat): First value in 6J symbol\n
      !!          j12(dfloat): Second value in 6J symbol\n
      !!          j13(dfloat): Third value in 6J symbol\n
      !!          j21(dfloat): Fourth value in 6J symbol\n
      !!          j22(dfloat): Fifth value in 6J symbol\n
      !!          j23(dfloat): Sixth value in 6J symbol\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      double precision function fun6j(j11,j12,j13,j21,j22,j23,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      double precision:: j11,j12,j13,j21,j22,j23

      ! Local

      integer:: l1,l2,l3,l4,l5,l6
      integer:: ij1,ij2,ij3,ij4,ij5,ij6
      integer:: ijm,ijm1,ijm2,ijm3,ijm4
      integer:: ijx,ijx1,ijx2,ijx3
      integer:: i

      double precision:: term1,term2

      double precision, pointer:: val


      ! If doing memoization
      if (Flgsg%memo) then

        ! Convert to integers
        l1 = nint(2d0*j11)
        l2 = nint(2d0*j12)
        l3 = nint(2d0*j13)
        l4 = nint(2d0*j21)
        l5 = nint(2d0*j22)
        l6 = nint(2d0*j23)

        ! Look for existing element
        val => elem6D(l1,l2,l3,l4,l5,l6,Flgsg%J6)

        ! If found, return value
        if (associated(val)) then
          fun6j = val
          return
        end if ! Found value in jagged array
      end if ! Memoization

      fun6j = 0d0

      ij1 = nint(j11+j11)
      ij2 = nint(j12+j12)
      ij3 = nint(j13+j13)
      ij4 = nint(j21+j21)
      ij5 = nint(j22+j22)
      ij6 = nint(j23+j23)

      ijm1 = (ij1+ij2+ij3)/2
      ijm2 = (ij1+ij5+ij6)/2
      ijm3 = (ij4+ij2+ij6)/2
      ijm4 = (ij4+ij5+ij3)/2

      ijm = ijm1

      ijm = max(ijm,ijm2,ijm3,ijm4)

      ijx1 = (ij1+ij2+ij4+ij5)/2
      ijx2 = (ij2+ij3+ij5+ij6)/2
      ijx3 = (ij3+ij1+ij6+ij4)/2

      ijx = ijx1

      ijx = min(ijx,ijx2,ijx3)

      if(ijm.le.ijx) then

        term1 = fn2(ij1,ij2,ij3,Flgsg) + fn2(ij1,ij5,ij6,Flgsg) + &
                fn2(ij4,ij2,ij6,Flgsg) + fn2(ij4,ij5,ij3,Flgsg)

        do i=ijm,ijx
          term2 = Flgsg%flg(i+1) - Flgsg%flg(i-ijm1) - &
                  Flgsg%flg(i-ijm2) - Flgsg%flg(i-ijm3) - &
                  Flgsg%flg(i-ijm4) - Flgsg%flg(ijx1-i) - &
                  Flgsg%flg(ijx2-i) - Flgsg%flg(ijx3-i)
          fun6j = Flgsg%sg(i)*exp(term1+term2) + fun6j
        end do

      end if

      ! If memoization, insert
#ifdef _OPENMP
      if (Flgsg%memo.and.Flgsg%can_write) &
        call insert6D(l1,l2,l3,l4,l5,l6,fun6j,Flgsg%J6)
#else
      if (Flgsg%memo) call insert6D(l1,l2,l3,l4,l5,l6,fun6j,Flgsg%J6)
#endif

      return

      end function fun6j

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculates 9j value, each jii is real\n
      !!          j11(dfloat): First value in 9J symbol\n
      !!          j12(dfloat): Second value in 9J symbol\n
      !!          j13(dfloat): Third value in 9J symbol\n
      !!          j21(dfloat): Fourth value in 9J symbol\n
      !!          j22(dfloat): Fifth value in 9J symbol\n
      !!          j23(dfloat): Sixth value in 9J symbol\n
      !!          j31(dfloat): Seventh value in 9J symbol\n
      !!          j32(dfloat): Eighth value in 9J symbol\n
      !!          j33(dfloat): Nineth value in 9J symbol\n
      !!   Flgsg(Fctsg_class): Structure with factorials and signs
      double precision function fun9j(j11,j12,j13,j21,j22,j23, &
                                      j31,j32,j33,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      double precision:: j11,j12,j13,j21,j22,j23,j31,j32,j33

      ! Local

      integer:: l1,l2,l3,l4,l5,l6,l7,l8,l9
      integer:: ij11,ij12,ij13,ij21,ij22,ij23,ij31,ij32,ij33
      integer:: kmin1,kmin2,kmin3,kmax1,kmax2,kmax3,k

      double precision:: hk

      double precision, pointer:: val

      ! If doing memoization
      if (Flgsg%memo) then

        ! Convert to integers
        l1 = nint(2d0*j11)
        l2 = nint(2d0*j12)
        l3 = nint(2d0*j13)
        l4 = nint(2d0*j21)
        l5 = nint(2d0*j22)
        l6 = nint(2d0*j23)
        l7 = nint(2d0*j31)
        l8 = nint(2d0*j32)
        l9 = nint(2d0*j33)

        ! Look for existing element
        val => elem9D(l1,l2,l3,l4,l5,l6,l7,l8,l9,Flgsg%J9)

        ! If found, return value
        if (associated(val)) then
          fun9j = val
          return
        end if ! Found value in jagged array
      end if ! Memoization

      fun9j = 0d0

      ij11 = nint(j11+j11)
      ij12 = nint(j12+j12)
      ij13 = nint(j13+j13)
      ij21 = nint(j21+j21)
      ij22 = nint(j22+j22)
      ij23 = nint(j23+j23)
      ij31 = nint(j31+j31)
      ij32 = nint(j32+j32)
      ij33 = nint(j33+j33)

      kmin1 = abs(ij11-ij33)
      kmin2 = abs(ij32-ij21)
      kmin3 = abs(ij23-ij12)

      kmin1 = max(kmin1,kmin2,kmin3)

      kmax1 = ij11+ij33
      kmax2 = ij32+ij21
      kmax3 = ij23+ij12

      kmax1 = min(kmax1,kmax2,kmax3)

      if(kmin1.le.kmax1) then

        do k=kmin1,kmax1,2
          hk = .5d0*dble(k)
          fun9j = Flgsg%sg(k)*dble(k+1)* &
                  fun6j(j11,j21,j31,j32,j33,hk,Flgsg)* &
                  fun6j(j12,j22,j32,j21,hk,j23,Flgsg)* &
                  fun6j(j13,j23,j33,hk,j11,j12,Flgsg) + &
                  fun9j
        end do

      end if

      ! If memoization, insert
#ifdef _OPENMP
      if (Flgsg%memo.and.Flgsg%can_write) &
        call insert9D(l1,l2,l3,l4,l5,l6,l7,l8,l9,fun9j,Flgsg%J9)
#else
      if (Flgsg%memo) &
        call insert9D(l1,l2,l3,l4,l5,l6,l7,l8,l9,fun9j,Flgsg%J9)
#endif

      return

      end function fun9j

!#####################################################################
!#####################################################################
!#####################################################################

      end module funnj_mod
