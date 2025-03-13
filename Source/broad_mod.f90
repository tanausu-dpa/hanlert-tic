      !> Broadening of transition lines
      module broad_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
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
!  broad
!    Compute line broadening for every line for a given atom
!
!  broad_line
!    Compute line broadening for every line for a given LTE line
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use broadaux_mod
      use commons_mod
      use parameters_mod , only : c , PI , ccons
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute line broadening for every line for a given atom\n
      !!      Atom(Atom_class): Structure with atomic data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!  folder(character(:)): Path to the output folder\n
      !!       aparam(logical): If the equivalant parametric
      !!                        parameters need to be stored
      subroutine broad(Atom,Atmo,folder,aparam)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      character(len=500), intent(in):: folder
      logical, intent(in):: aparam

      ! Local

      logical:: laparam

      integer:: ios,iterml,itermu,itran

      double precision, dimension(NZ):: damp


      ! Allocate and count memory to store the damping
      ! parameter and the elastic rates
      allocate(Atom%ldamp(Atom%ntran,NZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%ldamp)
      allocate(Atom%qel(Atom%ntran,NZ))
      MRAMc = MRAMc + 1d-6*sizeof(Atom%qel)

      ! Initialize to zero
      Atom%ldamp = 0d0
      Atom%qel = 0d0

      ! If Gaussian profiles
      if (VOITY.eq.3) then

        ! If 1D synthesis
        if (run_mode.eq.0) then

          ! Remove memory count for inputs
          MRAMc = MRAMc - 1d-6*sizeof(Atom%broad_args)
          MRAMc = MRAMc - 1d-6*sizeof(Atom%broad_args)

          ! We can drop the inputs
          deallocate(Atom%broad_args)
          deallocate(Atom%broad_stark)

        end if ! 1D synthesis

        ! Get out
        return

      end if ! Gaussian profiles

      ! If the we are storing the parameters (only active atoms can
      ! have aparam = .True.) and the master
      if (aparam.and.pid.eq.0) then

        ! Flag for calls below
        laparam = .True.

        ! Open files
        open(600, file=trim(folder)//'/'//trim(Atom%file_label)// &
             '.avdwparam', iostat=ios, err=1000)
        open(700, file=trim(folder)//'/'//trim(Atom%file_label)// &
             '.astkparam', iostat=ios, err=1001)

        ! Prepare headers
        write(600,'(A)',err=1100) &
              'This file contains the parameters that you have '// &
              'to put in the model atom with the parameter '// &
              'option to mimic the current Van dew Waals '// &
              'approximation'
        write(600,'(A)',err=1100) 'Atom: '//trim(Atom%Element)
        write(600,'(A)',err=1100) &
              'Transition                             approxim.'// &
              '             A_H             B_H            '// &
              'A_He            B_He'
        write(700,'(A)',err=1101) &
              'This file contains the parameters that you have '// &
              'to put in the model atom with the parameter '// &
              'option to mimic the current Stark approximation'
        write(700,'(A)',err=1101) 'Atom: '//trim(Atom%Element)
        write(700,'(A)',err=1101) &
              'Transition                                         C'

      ! If not storing or it is a slave
      else

        ! Flag for calls below
        laparam = .False.

      end if ! Storing damping parameters and is master


      ! For each transition
      do itran=1,Atom%ntran

        ! Get term indexes
        itermu = Atom%fst(itran)%itermu
        iterml = Atom%fst(itran)%iterml

        ! Initialize the local variable
        damp = 0d0

        ! Van der Waals
        call broad_vdw(Atom,Atmo,iterml,itermu,itran,damp,laparam)

        ! Quadratic Stark
        call broad_stk(Atom,Atmo,iterml,itermu,itran,damp,laparam)

        ! Linear Stark
        call broad_lstk(Atom,Atmo,iterml,itermu,damp)

        ! Add to the line broadening
        Atom%ldamp(itran,:) = Atom%ldamp(itran,:) + &
                              1d-16*damp/c/(4d0*PI)

        ! Elastic collisions
        Atom%qel(itran,:) = Atom%qel(itran,:) + damp

      end do ! Transitions

      ! If 1D synthesis mode
      if (run_mode.eq.0) then

        ! Remove memory count for inputs
        MRAMc = MRAMc - 1d-6*sizeof(Atom%broad_args)
        MRAMc = MRAMc - 1d-6*sizeof(Atom%broad_stark)

        ! Deallocate used inputs
        deallocate(Atom%broad_args)
        deallocate(Atom%broad_stark)

      end if ! 1D synthesis mode

      ! Close aparam files if they were opened
      if (laparam) then
        close(600)
        close(700) 
      end if

      ! Control
      call control

      return

1000  umsg = 'Error opening avdwparam file'
      urou = 'broad'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing avdwparam file'
      close(600)
      close(700)
      urou = 'broad'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening astkparam file'
      urou = 'broad'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing astkparam file'
      close(600)
      close(700)
      urou = 'broad'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine broad

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute line broadening for every line for a given LTE line\n
      !!  line(LTEline_class): Structure with the LTE line data\n
      !!     Atmo(Atmo_class): Structure with atmospheric data
      subroutine broad_line(line,Atmo)

      ! I/O

      type(LTEline_class), intent(inout):: line
      type(Atmo_class), intent(in):: Atmo

      ! Local

      double precision, dimension(:), allocatable:: damp


      ! Allocate, count memory, and initialize the variable to store
      ! the damping parameter
      allocate(line%damp(NZ))
      MRAMc = MRAMc + 1d-6*sizeof(line%damp)
      line%damp = 0d0

      ! Collisional (inelastic) contribution
      if (line%f_c.gt.0d0) then

        ! Oscillator strength formulation
        ! *1d-8 Units of 10^8 s^-1
        ! *1d-8 Factor for damping to compensate c scaling
        line%damp = line%damp + 1d-16*ccons*line%f_c*Atmo%ne* &
                    (1d0 + line%nu/line%nl)/(2d0*line%Ju+1d0)/ &
                    sqrt(Atmo%T)/c/(4d0*PI)

      end if

      ! Radiative additional damping
      if (line%broad_rad.gt.0d0) then

        ! Collisional contribution
        line%damp = line%damp + 1d-16*line%broad_rad/c/(4d0*PI)

      end if

      ! If not Gaussian
      if (VOITY.ne.3) then

        ! Initialize the local variable
        allocate(damp(nz))
        damp = 0d0

        ! Van der Waals contribution
        call broad_vdw_LTE(line,Atmo,damp)

        ! Stark contribution
        call broad_stk_LTE(line,Atmo,damp)

        ! Linear Stark contribution
        call broad_lstk_LTE(line,Atmo,damp)

        ! Add to the line broadening
        line%damp = line%damp + 1d-16*damp/c/(4d0*PI)

        ! Free
        deallocate(damp)

      end if ! Non-Gaussian profiles

      return

      end subroutine broad_line

!#####################################################################
!#####################################################################
!#####################################################################

      end module broad_mod
