      !> CHIANTI data manager
      module chianti_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     11/11/2022
!  Last version:
!     28/11/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     28/11/2024:    V4.0.0 - Updated calls to abortedS to not include
!                             thread information (TdPA)
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
!  rCHIANTI
!    Read the ionization fraction data from the CHIANTI database
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

      !> Read the ionization fraction data from the CHIANTI database\n
      !!      Input(Input_class): Structure with configuration data\n
      !!  chianti(chianti_class): Structure with the CHIANTI data
      subroutine rCHIANTI(Input,chianti)

      ! I/O

      type(Input_class), intent(in):: Input
      type(chianti_class), intent(out):: chianti

      ! Local

      integer:: ios,nT,nE,ndim,nels,i1,i2,iel,last
      real, dimension(:), allocatable:: buffer


      ! If no input
      if (trim(Input%chianti_path).eq.'NONE') then

          ! This signals the lack of data
          chianti%nT = -1
          chianti%nE = -1

          ! And exit
          return

      end if ! No Input

      !
      ! Master
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

        ! For each element
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
      ! Slaves
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

        ! Add memory count
        MRAMc = MRAMc + 1d-6*sizeof(chianti%ioneq_data)

        ! Receive data
        call MPI_BCAST(chianti%ioneq_data,chianti%Nioneq, &
                       MPI_REAL,0,MPI_COMM_WORLD,ios)

        ! Point temperature
        chianti%ioneq_T => chianti%ioneq_data(1:chianti%nT)

        ! Initialize last position
        last = chianti%nT

        ! For each element
        do i1=1,chianti%nE

          ! Add memory count
          MRAMc = MRAMc + 1d-6*sizeof(chianti%ioneq(i1))

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

          ! Add memory count
          MRAMc = MRAMc + 1d-6*sizeof(chianti%ioneq(i1)%b)
          MRAMc = MRAMc + 1d-6*sizeof(chianti%ioneq(i1)%c)
          MRAMc = MRAMc + 1d-6*sizeof(chianti%ioneq(i1)%d)

          ! For each ionization stage
          do i2=1,chianti%ioneq(i1)%nI

            ! Point to ionization data
            chianti%ioneq(i1)%stage(i2)%p =>  &
                chianti%ioneq_data(last+1:last+chianti%nT)

            ! Calculate coefficients for splines
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
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error reading CHIANTI ioneq file'
      urou = 'rCHIANTI'
      close(100)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine rCHIANTI

!#####################################################################
!#####################################################################
!#####################################################################

      end module chianti_mod
