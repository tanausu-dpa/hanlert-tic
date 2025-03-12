      !> n-J symbols
      module funnj_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     12/03/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     12/03/2025:    V4.0.1 - Gave access to common_mod (TdPA)
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
!  fn1
!    Auxiliar function used in the calculation of J-symbols
!
!  fn2
!    Auxiliar function used in the calculation of J-symbols
!
!  fun3j
!    Calculate a 3J-symbol
!
!  fun6j
!    Calculate a 6J-symbol
!
!  fun9j
!    Calculate a 9J-symbol
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use memoization_mod
      use types_mod

      ! Functions for n-j symbols

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Auxiliar function used in the calculation of J-symbols\n
      !!          j1(double): First value\n
      !!          j2(double): Second value\n
      !!          j3(double): Third value\n
      !!          m1(double): Fourth value\n
      !!          m2(double): Fifth value\n
      !!          m3(double): Sixth value\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      double precision function fn1(j1,j2,j3,m1,m2,m3,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      double precision, intent(in):: j1,j2,j3,m1,m2,m3

      ! Local

      integer:: l1,l2,l3,l4,l5,l6,l7,l8,l9,l10


      ! Combinations of inputs
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

      ! Combine factorials
      fn1 = .5d0*(Flgsg%flg(l1) + Flgsg%flg(l2) + Flgsg%flg(l3) - &
                  Flgsg%flg(l4) + Flgsg%flg(l5) + Flgsg%flg(l6) + &
                  Flgsg%flg(l7) + Flgsg%flg(l8) + Flgsg%flg(l9) + &
                  Flgsg%flg(l10))

      return

      end function fn1

!#####################################################################
!#####################################################################
!#####################################################################

      !> Auxiliar function used in the calculation of J-symbols\n
      !!        ij1(integer): First value\n
      !!        ij2(integer): Second value\n
      !!        ij3(integer): Third value\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      double precision function fn2(ij1,ij2,ij3,Flgsg)

      ! I/O

      type(Fctsg_class),intent(in):: Flgsg
      integer, intent(in):: ij1,ij2,ij3

      ! Local

      integer:: l1,l2,l3,l4


      ! Combinations of inputs
      l1 = (ij1+ij2-ij3)/2
      l2 = (ij2+ij3-ij1)/2
      l3 = (ij3+ij1-ij2)/2
      l4 = (ij1+ij2+ij3)/2+1

      ! Combine factorials
      fn2 = .5d0*(Flgsg%flg(l1) + Flgsg%flg(l2) + Flgsg%flg(l3) - &
                  Flgsg%flg(l4))

      return

      end function fn2

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate a 3J-symbol\n
      !!          j1(double): First value in 3J symbol\n
      !!          j2(double): Second value in 3J symbol\n
      !!          j3(double): Third value in 3J symbol\n
      !!          m1(double): Fourth value in 3J symbol\n
      !!          m2(double): Fifth value in 3J symbol\n
      !!          m3(double): Sixth value in 3J symbol\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      double precision function fun3j(j1,j2,j3,m1,m2,m3,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: j1,j2,j3,m1,m2,m3

      ! Local

      integer:: l1,l2,l3,l4,l5,l6,ij1,ij2,ij3,im1,im2,im3
      integer:: kmin,kmin1,kmin2,kmax,kmax1,kmax2,kmax3,i

      double precision:: term1,term2,sgn
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

      !
      ! Not in memoization
      !

      ! Initialize
      fun3j = 0d0

      ! Combinations of upper row
      ij1 = nint(j1+j1)
      ij2 = nint(j2+j2)
      ij3 = nint(j3+j3)

      ! 3J rules
      if(ij1+ij2-ij3.lt.0) return
      if(ij2+ij3-ij1.lt.0) return
      if(ij3+ij1-ij2.lt.0) return

      ! Combinations of lower row
      im1 = nint(m1+m1)
      im2 = nint(m2+m2)
      im3 = nint(m3+m3)

      ! 3J rules
      if(im1+im2+im3.ne.0) return
      if(abs(im1).gt.ij1) return
      if(abs(im2).gt.ij2) return
      if(abs(im3).gt.ij3) return

      ! Get lower limit summation
      kmin = (ij3-ij1-im2)/2
      kmin1 = kmin
      kmin2 = (ij3-ij2+im1)/2
      kmin = max(-min(kmin,kmin2),0)

      ! Get upper limit summation
      kmax = nint(j1+j2-j3)
      kmax1 = kmax
      kmax2 = nint(j1-m1)
      kmax3 = nint(j2+m2)
      kmax = min(kmax,kmax2,kmax3)

      ! If valid number of terms
      if(kmin.le.kmax) then

        ! Call auxiliar
        term1 = fn1(j1,j2,j3,m1,m2,m3,Flgsg)

        ! Get sign
        sgn = Flgsg%sg((ij1-ij2-im3)/2)

        ! Perform sum
        do i=kmin,kmax

          ! Combination of factorials
          term2 = Flgsg%flg(i) + Flgsg%flg(kmin1+i) + &
                  Flgsg%flg(kmin2+i) + Flgsg%flg(kmax1-i) + &
                  Flgsg%flg(kmax2-i) + Flgsg%flg(kmax3-i)

          ! Add contribution
          fun3j = Flgsg%sg(i)*exp(term1-term2) + fun3j

        end do ! Sum

        ! Apply sign
        fun3j = sgn*fun3j

      end if ! Valid number of terms in summation

      ! If memoization
      if (Flgsg%memo) then
          
        ! Insert in memoization
        call insert6D(l1,l2,l3,l4,l5,l6,fun3j,Flgsg%J3)

        ! Count memory
        ERAMc = ERAMc + 8d-6

      end if ! Memoization

      return

      end function fun3j

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate a 6J-symbol\n
      !!         j11(double): First value in 6J symbol\n
      !!         j12(double): Second value in 6J symbol\n
      !!         j13(double): Third value in 6J symbol\n
      !!         j21(double): Fourth value in 6J symbol\n
      !!         j22(double): Fifth value in 6J symbol\n
      !!         j23(double): Sixth value in 6J symbol\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      double precision function fun6j(j11,j12,j13,j21,j22,j23,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: j11,j12,j13,j21,j22,j23

      ! Local

      integer:: l1,l2,l3,l4,l5,l6,ij1,ij2,ij3,ij4,ij5,ij6
      integer:: ijm,ijm1,ijm2,ijm3,ijm4,ijx,ijx1,ijx2,ijx3,i

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

      !
      ! Not in memoization
      !

      ! Initialize
      fun6j = 0d0

      ! Combinations of inputs
      ij1 = nint(j11+j11)
      ij2 = nint(j12+j12)
      ij3 = nint(j13+j13)
      ij4 = nint(j21+j21)
      ij5 = nint(j22+j22)
      ij6 = nint(j23+j23)

      ! Lower limit in summation
      ijm1 = (ij1+ij2+ij3)/2
      ijm2 = (ij1+ij5+ij6)/2
      ijm3 = (ij4+ij2+ij6)/2
      ijm4 = (ij4+ij5+ij3)/2
      ijm = ijm1
      ijm = max(ijm,ijm2,ijm3,ijm4)

      ! Upper limit in summation
      ijx1 = (ij1+ij2+ij4+ij5)/2
      ijx2 = (ij2+ij3+ij5+ij6)/2
      ijx3 = (ij3+ij1+ij6+ij4)/2
      ijx = ijx1
      ijx = min(ijx,ijx2,ijx3)

      ! If valid number of terms
      if(ijm.le.ijx) then

        ! Call auxiliar function
        term1 = fn2(ij1,ij2,ij3,Flgsg) + fn2(ij1,ij5,ij6,Flgsg) + &
                fn2(ij4,ij2,ij6,Flgsg) + fn2(ij4,ij5,ij3,Flgsg)

        ! Perform sum
        do i=ijm,ijx

          ! Combination of factorials
          term2 = Flgsg%flg(i+1) - Flgsg%flg(i-ijm1) - &
                  Flgsg%flg(i-ijm2) - Flgsg%flg(i-ijm3) - &
                  Flgsg%flg(i-ijm4) - Flgsg%flg(ijx1-i) - &
                  Flgsg%flg(ijx2-i) - Flgsg%flg(ijx3-i)

          ! Add contribution
          fun6j = Flgsg%sg(i)*exp(term1+term2) + fun6j

        end do ! Sum

      end if ! Valid number of terms in summation

      ! If memoization
      if (Flgsg%memo) then
          
        ! Insert in memoization
        call insert6D(l1,l2,l3,l4,l5,l6,fun6j,Flgsg%J6)

        ! Count memory
        ERAMc = ERAMc + 8d-6

      end if ! Memoization

      return

      end function fun6j

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate a 9J-symbol\n
      !!         j11(double): First value in 9J symbol\n
      !!         j12(double): Second value in 9J symbol\n
      !!         j13(double): Third value in 9J symbol\n
      !!         j21(double): Fourth value in 9J symbol\n
      !!         j22(double): Fifth value in 9J symbol\n
      !!         j23(double): Sixth value in 9J symbol\n
      !!         j31(double): Seventh value in 9J symbol\n
      !!         j32(double): Eighth value in 9J symbol\n
      !!         j33(double): Nineth value in 9J symbol\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols
      double precision function fun9j(j11,j12,j13,j21,j22,j23, &
                                      j31,j32,j33,Flgsg)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(in):: j11,j12,j13,j21,j22,j23, &
                                     j31,j32,j33

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

      !
      ! Not in memoization
      !

      ! Initialize
      fun9j = 0d0

      ! Combinations of inputs
      ij11 = nint(j11+j11)
      ij12 = nint(j12+j12)
      ij13 = nint(j13+j13)
      ij21 = nint(j21+j21)
      ij22 = nint(j22+j22)
      ij23 = nint(j23+j23)
      ij31 = nint(j31+j31)
      ij32 = nint(j32+j32)
      ij33 = nint(j33+j33)

      ! Lower limit in summation
      kmin1 = abs(ij11-ij33)
      kmin2 = abs(ij32-ij21)
      kmin3 = abs(ij23-ij12)
      kmin1 = max(kmin1,kmin2,kmin3)

      ! Upper limit in summation
      kmax1 = ij11+ij33
      kmax2 = ij32+ij21
      kmax3 = ij23+ij12
      kmax1 = min(kmax1,kmax2,kmax3)

      ! If valid number of terms
      if(kmin1.le.kmax1) then

        ! Perform sum
        do k=kmin1,kmax1,2

          ! Half K
          hk = .5d0*dble(k)

          ! Add contribution
          fun9j = Flgsg%sg(k)*dble(k+1)* &
                  fun6j(j11,j21,j31,j32,j33,hk,Flgsg)* &
                  fun6j(j12,j22,j32,j21,hk,j23,Flgsg)* &
                  fun6j(j13,j23,j33,hk,j11,j12,Flgsg) + &
                  fun9j

        end do ! Sum

      end if ! Valid number of terms in summation

      ! If memoization, insert
      if (Flgsg%memo) then

        ! Insert in memoization
        call insert9D(l1,l2,l3,l4,l5,l6,l7,l8,l9,fun9j,Flgsg%J9)

        ! Count memory
        ERAMc = ERAMc + 8d-6

      end if ! Memoization

      return

      end function fun9j

!#####################################################################
!#####################################################################
!#####################################################################

      end module funnj_mod
