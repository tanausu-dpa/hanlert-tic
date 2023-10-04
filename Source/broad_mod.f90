      !> Broadening of transition lines
      module broad_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/19/2017
!  Last version:
!     09/25/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/25/2023:    V3.0.5 - Changed name of param files (TdPA)
!
!     08/07/2023:    V3.0.4 - Added broad_line routine (TdPA)
!
!     10/25/2022:    V3.0.3 - Added option for gaussian profiles
!                             which basically skips computing the
!                             damping parameter (TdPA)
!                           - Added a slash after output folder
!                             when saving file (TdPA)
!
!     07/13/2022:    V3.0.2 - The resource input parameter is no
!                             longer needed (TdPA)
!
!     07/08/2022:    V3.0.1 - Bugfix: Atom%broad_args and
!                             Atom%broad_stark can only be
!                             deallocated if we are doing a single 1D
!                             synthesis (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     11/19/2019:    V1.2.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     03/12/2019:    V1.2.0 - Can create a file with the broading
!                             parameters (TdPA)
!
!     02/20/2019:    V1.1.0 - New verbosity (TdPA)
!
!     09/15/2017:    V1.0.1 - Receiving Input%resource (TdPA)
!
!     04/19/2017:    V1.0.0 - First version (TdPA)
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
!
!  broad:
!    Calculates the damping parameters due to for the collisional
!  broadening of spectral lines
!
!  broad_line:
!    Calculate the damping parameter for LTE lines
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

      !> Calculates broadening contributions due to Van der Waals and
      !! Stark effects\n
      !!        Atom(Atom_class): Structure with the atomic data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    folder(character(:)): Path to the output folder\n
      !!         aparam(logical): Store VdW parametric quantities
      subroutine broad(Atom,Atmo,folder,aparam)

      ! I/O

      type(Atom_class), intent(inout):: Atom
      type(Atmo_class), intent(in):: Atmo
      character(len=500), intent(in):: folder
      logical, intent(in):: aparam

      ! Local

      logical:: laparam

      integer:: ios,iterm,iterm1,itran

      double precision, dimension(NZ):: damp


      ! Allocate the variable to store this damping parameter
      allocate(Atom%ldamp(Atom%ntran,NZ))
      Atom%ldamp = 0d0

      ! Gaussian?
      if (VOITY.eq.3) then
        ! Deallocate used inputs if we can free
        if (run_mode.eq.0) then
          deallocate(Atom%broad_args)
          deallocate(Atom%broad_stark)
        end if
        ! And get out
        return
      end if

      ! If the we are storing the parameters (only active can call
      ! with .True.), and only the master
      if (aparam.and.pid.eq.0) then

        laparam = .True.
        open(600, file=trim(folder)//'/'//trim(Atom%file_label)// &
             '.avdwparam', iostat=ios, err=1000)
        open(700, file=trim(folder)//'/'//trim(Atom%file_label)// &
             '.astkparam', iostat=ios, err=1001)

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

        laparam = .False.

      end if ! Storing damping parameters and is master


      ! For each pair of terms
      do iterm=1,Atom%nMulti-1
       do iterm1=iterm+1,Atom%nMulti

          ! Check if there is a transition
          itran = Atom%irad(iterm1,iterm)

          if (itran.le.0) cycle

          ! Initialize the local variable
          damp = 0d0

          ! Van der Waals contribution
          call broad_vdw(Atom,Atmo,iterm,iterm1,itran,damp, &
                         laparam)

          ! Stark contribution
          call broad_stk(Atom,Atmo,iterm,iterm1,itran,damp,laparam)

          ! Linear Stark contribution
          call broad_lstk(Atom,Atmo,iterm,iterm1,damp)

          ! Add to the line broadening
          Atom%ldamp(itran,:) = Atom%ldamp(itran,:) + &
                                1d-16*damp/c/(4d0*PI)

        end do
      end do


      ! Deallocate used inputs if we can free
      if (run_mode.eq.0) then
        deallocate(Atom%broad_args)
        deallocate(Atom%broad_stark)
      end if

      ! Close aparam file
      if (laparam) close(600)

      ! Control
      call control

      return

1000  umsg = 'Error opening avdwparam file'
      urou = 'broad'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1100  umsg = 'Error writing avdwparam file'
      close(600)
      close(700)
      urou = 'broad'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1001  umsg = 'Error opening astkparam file'
      urou = 'broad'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control
1101  umsg = 'Error writing astkparam file'
      close(600)
      close(700)
      urou = 'broad'
      call abortedS(umsg,urou,-1,.True.,.True.)
      call control

      end subroutine broad

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculates the total collisional broadening for a LTE line\n
      !!     line(LTEline_class): Structure with the LTE line data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data
      subroutine broad_line(line,Atmo)

      ! I/O

      type(LTEline_class):: line
      type(Atmo_class), intent(in):: Atmo

      ! Local
      double precision, dimension(:), allocatable:: damp

      ! Allocate the variable to store this damping parameter
      allocate(line%damp(NZ))
      line%damp = 0d0

      ! Collisional (inelastic)
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

      ! If not gaussian
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

      end if

      return

      end subroutine broad_line

!#####################################################################
!#####################################################################
!#####################################################################

      end module broad_mod
