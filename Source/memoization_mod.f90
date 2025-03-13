      !> 3J, 6J, and 9J RAM storage
      module memoization_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     \'Angel de Vicente (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     10/12/2019
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
!!$ * Copyright (C) 2019 by Angel de Vicente, angel.de.vicente@iac.es
!!$    https://github.com/angel-devicente/
!!$
!!$ ******************************************************************
!!$ * Angel de Vicente
!!$ *
!!$ * A simple program to perform memoization with dynamic jagged
!!$ * arrays
!!$ *
!!$ ******************************************************************
!!$ * TO DO:
!!$ *   + add proper description
!!$ *   + add error handling to code
!!$ ******************************************************************
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################


      !**************************************************************
      !
      ! Inserting functions for 1D-6D memoization jagged arrays
      !
      !***************************************************************

      subroutine insert_scalar(val, A)
      double precision :: val
      type(scalar) :: A

      if (associated(A%d)) then
        ! print*, "---> Repeated element. Ignoring: ", val
      else
        allocate(A%d)
        A%d = val
      end if
      end subroutine insert_scalar

      subroutine insert1D(l1, val, A)
      integer :: i, l1
      double precision :: val
      type(a1D) :: A
      type(scalar), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert_scalar(val,A%d(l1))
        else   ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert_scalar(val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert_scalar(val,A%d(l1))
      end if
      end subroutine insert1D

      subroutine insert2D(l1, l2, val, A)
      integer :: i, l1, l2
      double precision :: val
      type(a2D) :: A
      type(a1D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert1D(l2,val,A%d(l1))
        else   ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert1D(l2,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert1D(l2,val,A%d(l1))
      end if
      end subroutine insert2D

      subroutine insert3D(l1, l2, l3, val, A)
      integer :: i, l1,l2,l3
      double precision :: val
      type(a3D) :: A
      type(a2D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert2D(l2,l3,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert2D(l2,l3,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert2D(l2,l3,val,A%d(l1))
      end if
      end subroutine insert3D

      subroutine insert4D(l1, l2, l3, l4, val, A)
      integer :: i, l1,l2,l3,l4
      double precision :: val
      type(a4D) :: A
      type(a3D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert3D(l2,l3,l4,val,A%d(l1))
        else   ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert3D(l2,l3,l4,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert3D(l2,l3,l4,val,A%d(l1))
      end if
      end subroutine insert4D

      subroutine insert5D(l1, l2, l3, l4, l5, val, A)
      integer :: i, l1,l2,l3,l4,l5
      double precision :: val
      type(a5D) :: A
      type(a4D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert4D(l2,l3,l4,l5,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert4D(l2,l3,l4,l5,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert4D(l2,l3,l4,l5,val,A%d(l1))
      end if
      end subroutine insert5D

      subroutine insert6D(l1, l2, l3, l4, l5, l6,val, A)
      integer :: i, l1,l2,l3,l4,l5,l6
      double precision :: val
      type(a6D) :: A
      type(a5D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert5D(l2,l3,l4,l5,l6,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert5D(l2,l3,l4,l5,l6,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert5D(l2,l3,l4,l5,l6,val,A%d(l1))
      end if
      end subroutine insert6D

      subroutine insert7D(l1, l2, l3, l4, l5, l6, l7,val, A)
      integer :: i, l1,l2,l3,l4,l5,l6,l7
      double precision :: val
      type(a7D) :: A
      type(a6D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert6D(l2,l3,l4,l5,l6,l7,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert6D(l2,l3,l4,l5,l6,l7,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert6D(l2,l3,l4,l5,l6,l7,val,A%d(l1))
      end if
      end subroutine insert7D

      subroutine insert8D(l1, l2, l3, l4, l5, l6, l7, l8,val, A)
      integer :: i, l1,l2,l3,l4,l5,l6,l7,l8
      double precision :: val
      type(a8D) :: A
      type(a7D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert7D(l2,l3,l4,l5,l6,l7,l8,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert7D(l2,l3,l4,l5,l6,l7,l8,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert7D(l2,l3,l4,l5,l6,l7,l8,val,A%d(l1))
      end if
      end subroutine insert8D

      subroutine insert9D(l1, l2, l3, l4, l5, l6, l7, l8, l9,val, A)
      integer :: i, l1,l2,l3,l4,l5,l6,l7,l8,l9
      double precision :: val
      type(a9D) :: A
      type(a8D), dimension(:), pointer :: d_t

      if (associated(A%d)) then
        if (l1 >= lbound(A%d,1) .and. l1 <= ubound(A%d,1)) then
          call insert8D(l2,l3,l4,l5,l6,l7,l8,l9,val,A%d(l1))
        else ! not big enough: we need to realloc and place value
          if (l1 < lbound(A%d,1)) then ! we are extending the low part
            allocate(d_t(l1:ubound(A%d,1)))
            do i=l1,lbound(A%d,1)-1
              d_t(i)%d => null()
            end do
          else ! we are extending the high part
            allocate(d_t(lbound(A%d,1):l1))
            do i=ubound(A%d,1)+1,l1
              d_t(i)%d => null()
            end do
          end if
          do i=lbound(A%d,1),ubound(A%d,1)
            d_t(i)%d => A%d(i)%d
          end do
          deallocate(A%d)
          A%d => d_t
          call insert8D(l2,l3,l4,l5,l6,l7,l8,l9,val,A%d(l1))
        end if
      else ! dimension does not exist yet: create and place value
        if (l1 < 0) then
          allocate(A%d(l1:0))
        else
          allocate(A%d(0:l1))
        end if
        do i=lbound(A%d,1),ubound(A%d,1)
          A%d(i)%d => null()
        end do
        call insert8D(l2,l3,l4,l5,l6,l7,l8,l9,val,A%d(l1))
      end if
      end subroutine insert9D

      !**************************************************************
      !
      ! Search functions for 1D-6D memoization jagged arrays
      !
      ! Returns: NULL if element not in Jagged_Array.
      !          Pointer to value if present in Jagged_Array
      !
      !**************************************************************
      function elem1D(l1,A)
      integer :: l1
      type(a1D) :: A
      double precision, pointer :: elem1D
      nullify(elem1D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem1d => A%d(l1)%d
      end function elem1D

      function elem2D(l1,l2,A)
      integer :: l1,l2
      type(a2D) :: A
      double precision, pointer :: elem2D
      nullify(elem2D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem2D => elem1D(l2,A%d(l1))
      end function elem2D

      function elem3D(l1,l2,l3,A)
      integer :: l1,l2,l3
      type(a3D) :: A
      double precision, pointer :: elem3D
      nullify(elem3D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
      elem3D => elem2D(l2,l3,A%d(l1))
      end function elem3D

      function elem4D(l1,l2,l3,l4,A)
      integer :: l1,l2,l3,l4
      type(a4D) :: A
      double precision, pointer :: elem4D
      nullify(elem4D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem4D => elem3D(l2,l3,l4,A%d(l1))
      end function elem4D

      function elem5D(l1,l2,l3,l4,l5,A)
      integer :: l1,l2,l3,l4,l5
      type(a5D) :: A
      double precision, pointer :: elem5D
      nullify(elem5D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem5D => elem4D(l2,l3,l4,l5,A%d(l1))
      end function elem5D

      function elem6D(l1,l2,l3,l4,l5,l6,A)
      integer :: l1,l2,l3,l4,l5,l6
      type(a6D) :: A
      double precision, pointer :: elem6D
      nullify(elem6D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem6D => elem5D(l2,l3,l4,l5,l6,A%d(l1))
      end function elem6D

      function elem7D(l1,l2,l3,l4,l5,l6,l7,A)
      integer :: l1,l2,l3,l4,l5,l6,l7
      type(a7D) :: A
      double precision, pointer :: elem7D
      nullify(elem7D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem7D => elem6D(l2,l3,l4,l5,l6,l7,A%d(l1))
      end function elem7D

      function elem8D(l1,l2,l3,l4,l5,l6,l7,l8,A)
      integer :: l1,l2,l3,l4,l5,l6,l7,l8
      type(a8D) :: A
      double precision, pointer :: elem8D
      nullify(elem8D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem8D => elem7D(l2,l3,l4,l5,l6,l7,l8,A%d(l1))
      end function elem8D

      function elem9D(l1,l2,l3,l4,l5,l6,l7,l8,l9,A)
      integer :: l1,l2,l3,l4,l5,l6,l7,l8,l9
      type(a9D) :: A
      double precision, pointer :: elem9D
      nullify(elem9D)
      if (.not.associated(A%d)) return
      if ((l1 >= lbound(A%d,1)).and.(l1 <= ubound(A%d,1))) &
        elem9D => elem8D(l2,l3,l4,l5,l6,l7,l8,l9,A%d(l1))
      end function elem9D

!#####################################################################
!#####################################################################
!#####################################################################

      end module memoization_mod
