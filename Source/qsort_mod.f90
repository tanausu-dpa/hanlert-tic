      !> Sorting routine
      module qsort_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Made F conformant by Walt Brainerd
!
!  Modified:
!  Tanausú del Pino Alemán
!  Start:
!     04/19/2017
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
!
!     04/19/2016:    V1.0.0 - First version (TdPA)
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
!    Recursive Fortran 95 quicksort routine
!    sorts real numbers into ascending numerical order
!    Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
!    Based on algorithm from Cormen et al., Introduction to
!    Algorithms, 1997 printing
!
!#####################################################################
!#####################################################################
!#####################################################################

      public :: QsortC
      private :: Partition

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Main routine for sorting\n
      !!    A(dfloat(:)): Vector to sort and sorted vector
      recursive subroutine QsortC(A)

      ! I/O
      double precision, dimension(:),intent(inout):: A

      ! Local
      integer :: iq

      if (size(A) > 1) then
        call Partition(A,iq)
        call QsortC(A(:iq-1))
        call QsortC(A(iq:))
      endif

      end subroutine QsortC

!#####################################################################
!#####################################################################
!#####################################################################

      !> Partitions the vector to sort and marks the cut point\n
      !!       A(dfloat(:)): Vector to sort\n
      !!    marker(integer): Index where to cut
      subroutine Partition(A,marker)

      ! I/O

      integer, intent(out):: marker
      double precision, dimension(:),intent(inout):: A

      ! Local

      integer :: i, j
      double precision :: temp, x

      x = A(1)
      i = 0
      j = size(A) + 1

      do
        j = j - 1
        do
          if (A(j) <= x) exit
          j = j - 1
        end do
        i = i + 1
        do
          if (A(i) >= x) exit
          i = i + 1
        end do
        if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
        elseif (i == j) then
          marker = i + 1
          return
        else
          marker = i
          return
        endif
      end do

      end subroutine Partition

!#####################################################################
!#####################################################################
!#####################################################################

      end module qsort_mod
