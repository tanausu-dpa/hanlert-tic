      !> CHIANTI data manager
      module chianti_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     11/11/2022
!  Last version:
!     11/24/2022 V3.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     11/24/2022:    V3.0.0 - First version (TdPA)
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
!    Manages everything related to CHIANTI data
!
!    rCHIANTI
!      Read ionization fraction from the CHIANTI database
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use inter_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read CHIANTI data\n
      !!      Input(Input_class): Structure with settings data\n
      !!  chianti(chianti_class): Structure with the CHIANTI data
      subroutine rCHIANTI(Input,chianti)

      ! I/O

      type(Input_class), intent(in):: Input
      type(chianti_class), intent(out):: chianti

      ! Local

      integer:: ios,nT,nE,ndim,nels,i1,i2,iel,last
      real, dimension(:), allocatable:: buffer


      !
      ! If no input, exit
      !
      if (trim(Input%chianti_path).eq.'NONE') then

          ! This signals the lack of data
          chianti%nT = -1
          chianti%nE = -1
          return

      end if

      !
      ! Master reads and shared
      !
      if (gpid.eq.0) then

        ! Open file
        open(100, &
             file=trim(Input%chianti_path)//'ioneq/chianti.ioneq', &
             status='old', iostat=ios, err=1000)

        ! Read dimensions
        read(100,*,err=1100) nT,nE

        ! Number of entries and dimension of data
        nels = (nE+1)*nE/2 + nE
        ndim = (nels + 1)*nT

        ! Initialize buffer
        allocate(buffer(ndim))

        ! Read temperature
        read(100,*,err=1100) buffer(1:nT)

        ! Initialize last position
        last = nT

        ! Read elements
        do iel=1,nels

          ! Read into buffer
          read(100,*,err=1100) i1,i2,buffer(last+1:last+nT)

          ! Update last
          last = last + nT

        end do ! Read elements

        ! Close file
        close(100)

        ! Share with slaves
        call MPI_BCAST(nE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
        call MPI_BCAST(nT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
        call MPI_BCAST(buffer,ndim,MPI_REAL,0,MPI_COMM_WORLD,ios)

        ! Free data
        deallocate(buffer)

      !
      ! Slaves receive and store
      !
      else

        ! Receive dimensions
        call MPI_BCAST(chianti%nE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
        call MPI_BCAST(chianti%nT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)

        ! Compute dimension and allocate data
        chianti%Nioneq = ((chianti%nE+1)*chianti%nE/2 + &
                          chianti%nE + 1)*chianti%nT

        ! Allocate data and first level of pointers
        allocate(chianti%ioneq_data(chianti%Nioneq))
        allocate(chianti%ioneq(chianti%nE))

        ! Receive data
        call MPI_BCAST(chianti%ioneq_data,chianti%Nioneq, &
                       MPI_REAL,0,MPI_COMM_WORLD,ios)

        ! Point temperature
        chianti%ioneq_T => chianti%ioneq_data(1:chianti%nT)

        ! Initialize last position
        last = chianti%nT

        ! For each element
        do i1=1,chianti%nE

          ! Number of stages
          chianti%ioneq(i1)%nI = i1+1

          ! Allocate second level of pointers
          allocate(chianti%ioneq(i1)%stage(chianti%ioneq(i1)%nI))

          ! Prepare splines
          allocate(chianti%ioneq(i1)% &
                           b(chianti%nT,chianti%ioneq(i1)%nI))
          allocate(chianti%ioneq(i1)% &
                           c(chianti%nT,chianti%ioneq(i1)%nI))
          allocate(chianti%ioneq(i1)% &
                           d(chianti%nT,chianti%ioneq(i1)%nI))

          ! For each ionization
          do i2=1,chianti%ioneq(i1)%nI

            ! Point to ionization data
            chianti%ioneq(i1)%stage(i2)%p =>  &
                chianti%ioneq_data(last+1:last+chianti%nT)

            ! Get splines
            call spline(dble(chianti%ioneq_T), &
                        dble(chianti%ioneq(i1)%stage(i2)%p), &
                        chianti%ioneq(i1)%b(:,i2), &
                        chianti%ioneq(i1)%c(:,i2), &
                        chianti%ioneq(i1)%d(:,i2), &
                        chianti%nT)

            ! Update last
            last = last + chianti%nT

          end do ! Ion stage
        end do ! Element

      end if ! Master/slave

      return

1000  umsg = 'Error opening CHIANTI ioneq file'
      urou = 'rCHIANTI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error reading CHIANTI ioneq file'
      urou = 'rCHIANTI'
      close(100)
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine rCHIANTI

!#####################################################################
!#####################################################################
!#####################################################################

      end module chianti_mod
