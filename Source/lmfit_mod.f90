      !> Levenberg-Marquardt fit
      module lmfit_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/22/2023
!  Last version:
!     05/28/2024 V3.0.19
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     05/28/2024:   V3.0.19 - Moved the deallocation of
!                             Trial_Synthesis local variables before
!                             the return in case of error (TdPA)
!                           - If the forward solution fails during
!                             a trial in the backtracking, the
!                             algorithm will try to keep going with
!                             different lambda parameters. Have not
!                             checked the generality (TdPA)
!
!     05/20/2024:   V3.0.18 - Introduce an additional factor to
!                             scale regularization penalties (TdPA)
!
!     05/13/2024:   V3.0.17 - Completed the "shake-up" of the
!                             backtracking with parameters that can be
!                             changed in the input file (TdPA)
!
!     05/08/2024:   V3.0.16 - Testing new "shake-ups" to the
!                             backtracking (TdPA)
!
!     05/07/2024:   V3.0.15 - Added the possibility to track the
!                             lambda constant in the LM optimization
!                             to decide its initial value in the
!                             next iteration (TdPA)
!                           - Added the option to use 'emergency'
!                             measures when the backtracking gets
!                             stuck. This option is for experimenting
!                             and not ready to be documented for
!                             general use (TdPA)
!                           - Added predict_lambda and update_lambda
!                             subroutines (TdPA)
!
!     11/27/2023:   V3.0.14 - Wrong initialization of the variable
!                             to ensure that the error gets
!                             calculated (TdPA)
!
!     11/24/2023:   V3.0.13 - Added argument to LMFIT to skip the
!                             iterations for the pixel (TdPA)
!
!     10/04/2023:   V3.0.12 - Ensure that the '*' character before an
!                             iteration is also in the secondary
!                             verbosity file (TdPA)
!                           - If the solution is stored is decided by
!                             the new saving argument in LMFIT (TdPA)
!                           - Added new type of inversion (TdPA)
!
!     09/28/2023:   V3.0.11 - Add print priority to the '*' character
!                             before an iteration line (TdPA)
!                           - Conclusion of the bracketing is not
!                             printed with '+' symbol, similarly to
!                             the trials (TdPA)
!
!     09/26/2023:   V3.0.10 - Bugfix: There was no print of chi2
!                             in 15D mode (TdPA)
!
!     09/08/2023:    V3.0.9 - Verbosity update (TdPA)
!
!     08/17/2023:    V3.0.8 - Set penalty to 0 in backtracking when
!                             there is no regularization (TdPA)
!
!     07/31/2023:    V3.0.7 - Change the verbosity level in the
!                             inversion (HL)
!
!     07/03/2023:    V3.0.6 - Added reallocation of LM_Stru arrays
!                             at the beginning in case there was
!                             a failure doing a previous pixel (TdPA)
!                           - Initialize the synthesis solution as
!                             empty (TdPA)
!                           - No more file renaming, now set a given
!                             solution as best in RAM. To this end,
!                             the set_best routine has been added
!                             to the module (TdPA)
!                           - Added bound folding for velocity
!                             azimuth (TdPA)
!                           - If the iteration are exhausted, the
!                             check happens after checking for
!                             convergence (TdPA)
!                           - Only cycle if overflown Broyden if
!                             iterations are not exhausted (TdPA)
!                           - The results are written at the end
!                             of lmfit (TdPA)
!                           - Added verbosity of the model in
!                             Trial_synthesis, as well as a true
!                             copy (and freeing) of the model
!                             atmosphere (TdPA)
!
!     06/12/2023:    V3.0.5 - Rename the variable Inf_File (HL)
!                           - Keep penalty if the ratio is smaller
!                             than one (HL)
!
!     05/16/2023:    V3.0.4 - Added alternative method to compute the
!                             errors (TdPA)
!
!     04/25/2023:    V3.0.3 - Added the option to write the response
!                             functions into a file (TdPA)
!                           - Bugfix: The best chi2 was not being
!                             updated when there was an improvement
!                             after finding a worse chi2, stalling
!                             the backtracking and making it fail
!                             when it should not; this was also
!                             solved by HL in his own branch, merged
!                             before, but I changed it because it
!                             requires less lines (TdPA)
!
!     04/11/2023:    V3.0.2 - Verbose the lambda factor (HL)
!                           - Bugfix: update best index and chi2 (HL)
!
!     03/15/2023:    V3.0.1 - Made the limit to decide the initial
!                             penalty ratio parameters in module
!                             parameters_mod (TdPA)
!                           - Intpol_Atmo_all does not need the Flgsg
!                             argument (TdPA)
!                           - Cleaned not used variables (TdPA)
!                           - Removed unecessary broadcasts (TdPA)
!                           - Removed some commented lines (TdPA)
!                           - Removed some false branching
!                             conditionals (TdPA)
!                           - The Blos variables are in the same
!                             structure than the polar ones (TdPA)
!                           - The limits on the LM lambda, and the
!                             accepted and rejected factors are
!                             decided in the input, instead of
!                             hard-coded (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/23/2023:    V0.0.0 - Started from 12/12/2021
!                             TIC@lmfit_mod.f90 revision and
!                             05/12/2020 TIC@propose_mod.f90 revision
!                             from Hao (TdPA)
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
!    LMFIT:
!      Levenberg-Marquardt fit of Stokes parameters
!
!    Convergence_Check:
!      Check if LM problem is converged
!
!    Err:
!      Compute the error in the parameters
!
!    Trial_Synthesis:
!      Do a LM step and get new Stokes profiles
!
!    Backtracking:
!      Optimize lambda parameter with backtracking algorithm
!
!    Lambda_propose_fix:
!      Propose a lambda factor for the LM
!
!    predict_lambda:
!      Propose the value of lambda for the next iteration
!
!    update_lambda:
!      Update the tracked values of the lambda constant
!
!    Nodes_Modify:
!      Modify the node values according to the SVD solution
!
!    CheckLambda:
!      Check limits of lambda factor
!
!    set_best:
!      Sets the current solution as the best of the backtracking
!    or the global problem
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use bounds_mod
      use commons_mod
      use free_mod
      use hanlert_mod
      use initmodel_mod
      use inter_mod
      use iosolution_mod
      use iotic_mod
      use jacobian_mod
      use model_mod
      use parameters_mod , only: TINYPT, TINYREG
      use ratmo_mod
      use regul_mod
      use svd_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Levenberg-Marquardt fit of Stokes parameters
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!         Input(Input_class): Structure with settings data\n
      !!   Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                             data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!       LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                             other LM quantities\n
      !!             imask(integer): Indicate if this pixel is masked
      !!                             in the restart\n
      !!            saving(logical): If the result is to be stored
      subroutine LMFIT(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                       kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                       Inf_Nodes,Sol,LM_Stru,imask,saving)


      ! IO
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(LMFIT_class), intent(inout):: LM_Stru
      logical, intent(in):: saving
      integer, intent(in):: imask

      ! Local

      ! RT full solution
      type(Solution_F_class):: SolF

      ! Hessian computing
      double precision, dimension(:), allocatable:: Solution, Errors
      double precision, dimension(:,:), allocatable:: Stokes_Min
      double precision, dimension(:,:), allocatable:: Stokes_best

      ! Lambda boundary
      logical:: Flag_Convg, Flag_Jac

      integer:: indx_iter, indx_rej, i, Num_Broyden, max_iters

      double precision:: Chisq_old, Ratio
      double precision, dimension(:), allocatable:: Lam_track


      ! Check allocations (in case of previous failure)
      if (allocated(LM_Stru%Hessian)) then

        if (allocated(LM_Stru%ResidualI)) then
          deallocate(LM_Stru%ResidualI)
          deallocate(LM_Stru%WeightI)
          deallocate(LM_Stru%JacobianI)
        end if

        if (allocated(LM_Stru%Residual)) then
          deallocate(LM_Stru%Residual)
          deallocate(LM_Stru%Weight)
          deallocate(LM_Stru%Jacobian)
        end if

        ! Free memory
        deallocate(LM_Stru%Jacfvec,LM_Stru%diag)
        deallocate(LM_Stru%Hessian,LM_Stru%Hessian_og)
        deallocate(LM_Stru%Jacfvec_og)

        ! Free solutions
        call free_inv_solution(SolF)

      end if


      ! If Thermal inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Allocate arrays for just intensity
        allocate(LM_Stru%ResidualI(Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%WeightI(Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%JacobianI(Sol%Num_Wavelength,LM_Stru%Num))

      ! Magnetic or full inversion
      else

        ! Allocate full Stokes arrays
        allocate(LM_Stru%Residual(0:3,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%Weight(0:3,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%Jacobian(0:3,Sol%Num_Wavelength,LM_Stru%Num))

      end if ! Type of inversion

      ! Initialize ratio for the first iteration
      Ratio = 0.1d0

      ! If ratio larger than regularitaion limit, reduce it now
      if (Ratio.gt.Input%Regul_Limit) Ratio = Input%Regul_Limit

      ! Initialize as false
      LM_Stru%Flag_weight = .False.

      ! Flag solution buffer not initialized
      SolF%no_initialized = .True.

      ! Allocate solution and Hessian quantities
      allocate(Solution(LM_Stru%Num),Errors(LM_Stru%Num))
      allocate(LM_Stru%Jacfvec(LM_Stru%Num),LM_Stru%diag(LM_Stru%Num))
      allocate(LM_Stru%Hessian(LM_Stru%Num,LM_Stru%Num))
      allocate(LM_Stru%Hessian_og(LM_Stru%Num,LM_Stru%Num))
      allocate(LM_Stru%Jacfvec_og(LM_Stru%Num))

      ! Allocate lambda array
      if (Input%l_Lam_track) then
        allocate(Lam_track(Input%Lam_track))
        Lam_track = 0d0
      end if

      ! Master
      if (pid.eq.0) then

        ! Verbose
        if (gpid.eq.0) then
          umsg = '-------------------'
          call verbose
          umsg = '| First synthesis |'
          call verbose
          umsg = '-------------------'
          call verbose
        end if
        umsg = ' - First synthesis'
        call verboseI(3)

      end if


      ! Get synthesis results
      call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                      kurucz,MPID,Atmo,Bfield,Input,Sol,SolF,.False.)

      ! Set best
      call set_best(SolF,.True.,.False.)

      ! Compute merit function
      call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                          Inf_Nodes%Nodes_Type,LM_Stru)

      ! If regularizing
      if (Inf_Nodes%Regul_Flag) then

        ! Allocate space for regularization
        allocate(LM_Stru%Rgl%Regul_H(LM_Stru%Num,LM_Stru%Num))
        allocate(LM_Stru%Rgl%Regul_F(LM_Stru%Num))

        ! Get regularization
        call Get_Regl_all(Inf_Nodes,.False.,LM_Stru%Rgl)

        ! Factorize
        LM_Stru%Rgl%Penalty = LM_Stru%Rgl%Penalty*Input%Regul_factor

        ! If penalty and limit larger than the minimum penalty
        ! and the regularization limit is larger than the
        ! minimum
        if (LM_Stru%Rgl%Penalty.gt.TINYPT.and. &
            Input%Regul_Limit.gt.TINYREG) then

          ! Update ratio to scale penalty to 0.1*chi2
          LM_Stru%Rgl%Ratio = Ratio/LM_Stru%Rgl%Penalty* &
                              LM_Stru%Chisq*Input%Regul_factor

          ! Keep penalty if the ratio is smaller than 1
          if (LM_Stru%Rgl%Ratio.lt.1d0) LM_Stru%Rgl%Ratio = 1d0

        ! Too small penalty or regularization limit
        else

          ! Keep penalty as it is
          LM_Stru%Rgl%Ratio = 1.0

        end if

        ! Update chi2
        LM_Stru%Chisq = LM_Stru%Chisq + &
                        LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

      ! Not regularizing
      else

        ! No penalty
        LM_Stru%Rgl%Ratio = 0d0
        LM_Stru%Rgl%Penalty = 0d0

      end if ! Regularizing

      ! Store current chi2
      Chisq_old = LM_Stru%Chisq

      ! Master
      if (pid.eq.0) then

        ! Verbose merit function if in output
        umsg = ' * '
        call verboseI(0)
        if (vlevel.eq.0) call verboseI(3)

        ! Verbose merit function
        write(umsg,'(A,i4,3x,A,es15.4)')  &
           ' * Iteration = ',0, &
              'Chi2 = ',LM_Stru%Chisq

        if (gpid.eq.0) then
          call verboseI(0)
          call verboseI(4)
        else
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)
        end if

        ! Save first chi2
        LM_Stru%Chisq_0 = LM_Stru%Chisq

        ! If regularizing
        if (Inf_Nodes%Regul_Flag) then

          ! Verbosity
          write(umsg,'(A)') ' - Initial regularization:'
          call verboseI(3)
          write(umsg,'(A,es15.4,2x,A,es15.4,2x,A,es15.4)') &
            '   Scaled penalty = ', &
            LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio, &
            'Ratio = ',LM_Stru%Rgl%Ratio, &
            'Penalty = ',LM_Stru%Rgl%Penalty
          call verboseI(3)

        end if ! Regularizing

        ! Verbose the current nodes
        call Verbose_Model(Inf_Nodes,.False.)

      end if ! Master

      ! Initialize LM variables
      LM_Stru%diag = 1d0
      Flag_Jac = .True.
      Num_Broyden = 0
      LM_Stru%accepted = .False.
      Flag_Convg = .False.
      ! Type of LM method
      select case(Input%LM_Method)
        ! Classic
        case(0)
          LM_Stru%Lambda = 1d-3
        ! Backtracking
        case(1)
          LM_Stru%Lambda = 1d-1*Input%factoraccept
      end select
      LM_Stru%Lambda_bounds = Input%Lam_Range
      LM_Stru%Flag_weight = .False.
      LM_Stru%factorreject = Input%factorreject
      LM_Stru%factoraccept = Input%factoraccept
      LM_Stru%accepted = .True.
      Stokes_best = Sol%Stokes_out

      ! If masked
      if (imask.eq.1) then

        ! Don't iterate and assume not accepted
        max_iters = 0

      ! Not masked
      else

        ! True number of iterations
        max_iters = Input%Num_Iter

      end if ! Mask

      ! Iterate, in principle, to the maximum allowed
      do indx_iter=1,max_iters

        ! Predict lambda
        if (Input%l_Lam_track) &
          call predict_lambda(indx_iter,Lam_track, \
                              Input%Lam_track,LM_Stru)

        ! If after the first iteration the accepted flag is on,
        ! we are doing Broyden and we have not done three already
        if (Input%Broyden.and.LM_Stru%accepted.and. &
            Num_Broyden.lt.3.and.indx_iter.gt.1) then

          ! Call Broyden
          call Broyden_Rank1(Stokes_best, Stokes_Min, &
                             Inf_Stokes%Num_Wavelength, &
                             Solution, Inf_Nodes, LM_Stru)
          ! Deflag Jacobian
          Flag_Jac = .False.

        ! First iteration, not accepted LM, not doing Broyden,
        ! or already three Broyden
        else

          ! Get Jabocian
          call Jacobian_Compute(Input,Atom,Atomb,Mol,Geom, &
                                GeomI,Flgsg,Frec,fudge,kurucz, &
                                MPID,Atmo,Bfield,Inf_Nodes,Sol, &
                                SolF,LM_Stru)
          if (laborted) return

          ! Flag Jacobian
          Flag_Jac = .True.

        end if

        ! Compute Hessian
        call Hessian_Compute(Inf_Nodes%Nodes_Type, &
                             LM_Stru)

        ! If regularizing
        if (Inf_Nodes%Regul_Flag) then

          ! Get regularization contribution
          call Get_Regl_all(Inf_Nodes, .True., LM_Stru%Rgl)

          ! Add to Hessian
          LM_Stru%Hessian = LM_Stru%Hessian + &
                            LM_Stru%Rgl%Regul_H*LM_Stru%Rgl%Ratio
          ! Add to vector
          LM_Stru%Jacfvec = LM_Stru%Jacfvec + &
                            LM_Stru%Rgl%Regul_F*LM_Stru%Rgl%Ratio

        end if

        ! Diagonal from Hessian
        do i=1,LM_Stru%Num
          LM_Stru%diag(i) = LM_Stru%Hessian(i,i)
        end do

        ! Type of LM method
        select case(Input%LM_Method)

          ! Classic
          case(0)

            ! Run up to the maximum number of rejections
            do indx_rej=1,10

              ! Master
              if (pid.eq.0) then

                ! Verbose
                if (gpid.eq.0) then
                  umsg = '---------------------------------------'
                  call verbose
                  umsg = '| Trial synthesis for tradiational LM |'
                  call verbose
                  umsg = '---------------------------------------'
                  call verbose
                end if
                umsg = ' - Trial synthesis for tradiational LM'
                call verboseI(3)

              end if ! Master

              ! Get trial synthesis
              call Trial_Synthesis(LM_Stru, Inf_Nodes, &
                                   LM_Stru%Lambda, Solution, &
                                   Atmo, Bfield, Atom, Atomb, &
                                   Mol, Geom, GeomI, Flgsg, Frec, &
                                   fudge, kurucz, MPID, Input, &
                                   Sol, SolF)
              if (laborted) return

              ! Get merit function
              call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                                  Inf_Nodes%Nodes_Type, LM_Stru)

              ! If regularizing, add its penalty
              if (Inf_Nodes%Regul_Flag) &
                LM_Stru%Chisq = LM_Stru%Chisq + &
                                LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

              ! If improved the chi2
              if (LM_Stru%Chisq.lt.Chisq_old) then

                ! Accept the step
                LM_Stru%accepted = .True.
                Chisq_old = LM_Stru%Chisq
                Stokes_best = Sol%Stokes_out

                ! Modify nodes
                call Nodes_Modify(Solution, Inf_Nodes)

                ! Fold B azimuth
                i = Inf_Nodes%index_Bp
                call FoldBounds(Inf_Nodes%Node(i), &
                                Inf_Nodes%Num_Nodes(i))

                ! Fold v azimuth
                i = Inf_Nodes%index_vy
                if (Inf_Nodes%vtype.eq.1) &
                  call FoldBounds(Inf_Nodes%Node(i), &
                                  Inf_Nodes%Num_Nodes(i))

                ! Generate new model stratification
                call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield, &
                                     Atom,Atomb,Mol,Input,fudge)

                ! Set best solution
                call set_best(SolF,.True.,.False.)

                ! If Jacobian flagged
                if (Flag_Jac) then

                  ! Reset Broyden counter
                  Num_Broyden = 0

                ! Otherwise
                else

                  ! Add to counter
                  Num_Broyden = Num_Broyden + 1

                end if

                ! Leave loop
                exit

              ! Did not improve
              else

                ! Reject the step
                LM_Stru%accepted = .False.
                Sol%Stokes_out = Stokes_best

                ! If the linear method to stimate the Jacobian
                ! is not working, leave to calculate it
                if (.not.Flag_Jac.and.indx_rej.ge.3) exit

              end if ! Improvement of chi2

              ! Propose a fix to the LM lambda parameter, as we
              ! did not improve the chi2
              call Lambda_propose_fix(LM_Stru, LM_Stru%Lambda)

              ! If lambda below lower limit, keep it there
              if (LM_Stru%Lambda.le.LM_Stru%Lambda_bounds(1)) &
                LM_Stru%Lambda = LM_Stru%Lambda_bounds(1)

              ! If lambda above upper limit, leave
              if (LM_Stru%Lambda.ge.LM_Stru%Lambda_bounds(2)) exit

            end do ! Allowed rejections

          ! Backtracking
          case(1)

            ! Optimize lambda parameter with backtracking algorithm
            call Backtracking(LM_Stru, Atom, Atomb, Mol, Geom, &
                              GeomI, Flgsg, Frec, fudge, kurucz, &
                              MPID, Atmo, Bfield, Sol, SolF, &
                              Inf_Nodes, Input, Inf_Stokes, &
                              Solution, Stokes_Min)

            ! If improved the chi2
            if (LM_Stru%Chisq.lt.Chisq_old) then

              ! Set best solution
              call set_best(SolF,.True.,.True.)

              ! Accept the step
              LM_Stru%accepted = .True.
              Stokes_best = Stokes_Min

              ! Modify nodes
              call Nodes_Modify(Solution, Inf_Nodes)

              ! Fold B azimuth
              i = Inf_Nodes%index_Bp
              call FoldBounds(Inf_Nodes%Node(i), &
                              Inf_Nodes%Num_Nodes(i))

              ! Fold v azimuth
              i = Inf_Nodes%index_vy
              if (Inf_Nodes%vtype.eq.1) &
                call FoldBounds(Inf_Nodes%Node(i), &
                                Inf_Nodes%Num_Nodes(i))

              ! Get new model stratification
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield, &
                                   Atom,Atomb,Mol,Input,fudge)

              ! If Jacobian flagged
              if (Flag_Jac) then

                ! Reset Broyden counter
                Num_Broyden = 0

              ! Otherwise
              else

                ! Add to counter
                Num_Broyden = Num_Broyden + 1

              end if

            ! Did not improve
            else

              ! Reject
              LM_Stru%accepted = .False.

              ! Get new model stratification
              call Intpol_Atmo_all(Inf_Nodes,Atmo,Bfield, &
                                   Atom,Atomb,Mol,Input,fudge)


              ! If Jacobi flagged
              if (Flag_Jac) then

                ! Revert Stokes
                Sol%Stokes_out = Stokes_best

                ! Leave iterations
                exit

              ! Otherwise
              else

                ! Overflow Broyden index
                Num_Broyden = 10

              end if ! Jacobi flagged
            end if ! Improved chi2

        end select ! LM method

        ! Set Stokes to best
        Sol%Stokes_out = Stokes_best

        ! If last step was accepted
        if (LM_Stru%accepted) then

          ! Get chi2 without regularization
          LM_Stru%Chisq_og = LM_Stru%Chisq - &
                             LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

          ! Master
          if (pid.eq.0) then

            ! Verbose merit function if in output
            umsg = ' * '
            call verboseI(0)
            if (vlevel.eq.0) call verboseI(3)

            ! Verbose merit function
            write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Iteration = ',indx_iter, &
              'Chi2 = ',LM_Stru%Chisq, &
              'Chi2 (no regularization) = ',LM_Stru%Chisq_og

            if (gpid.eq.0) then
              call verboseI(0)
              call verboseI(4)
            else
              call verboseI(0)
              if (vlevel.eq.0) call verboseI(3)
            end if


            ! If regularizing
            if (Inf_Nodes%Regul_Flag) then

              ! Verbose
              write(umsg,'(A)') ' - Regularization:'
              call verboseI(3)
              write(umsg,'(A,es15.4,2(3x,A,es15.4))') &
                '   Scaled penalty = ', &
                LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio, &
                'Ratio = ',LM_Stru%Rgl%Ratio, &
                'Penalty = ',LM_Stru%Rgl%Penalty
              call verboseI(3)

            end if

            ! Call verbose
            call Verbose_Model(Inf_Nodes,.False.)

          end if

          ! Check convergence
          call Convergence_Check(Input, Chisq_old, &
                                 LM_Stru%Chisq, Flag_Convg)

          ! If converged
          if (Flag_Convg) then

            ! If Flagged Jacobian
            if (Flag_Jac) then

              ! Master
              if (pid.eq.0) then

                ! Verbose merit function if in output
                umsg = ' * '
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)

                ! Verbose
                write(umsg,'(A,es15.4)') &
                  ' * LM converged. Chi2 = ',LM_Stru%Chisq

                if (gpid.eq.0) then
                  call verboseI(0)
                  call verboseI(4)
                else
                  call verboseI(0)
                  if (vlevel.eq.0) call verboseI(3)
                end if

              end if ! Master

              ! Leave iterations, we are done
              exit

            ! Not flagged Jacobian
            else

              ! Overflow Broyden and do more
              Num_Broyden = 10
              if (indx_iter.lt.Input%Num_iter) cycle

            end if ! Flagged Jacobian
          end if ! Converged

          ! If last iteration, leave
          if (indx_iter.eq.Input%Num_Iter) then

            ! Master
            if (pid.eq.0) then

              ! Verbose
              umsg = ' - LM iterations exhausted'

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if ! Master

            ! Exit loop
            exit

          end if ! Master

          ! If regularizing and the regularization ratio is larger
          ! than 1
          if (Inf_Nodes%Regul_Flag.and.LM_Stru%Rgl%Ratio.gt.1d0) then

            ! Update ratio
            Ratio = (Input%Regul_Limit)/ &
                    LM_Stru%Rgl%Penalty*LM_Stru%Chisq* &
                    Input%Regul_factor

            ! Keep ratio at least 1
            if (Ratio.lt.1d0) Ratio = 1d0

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,es15.4)') ' - Current Ratio = ', Ratio
              call verboseI(3)

            end if

            ! If ratio is below the stored one and the resulting
            ! penalty is above the chosen limit
            if (Ratio.lt.LM_Stru%Rgl%Ratio) then

              ! Update ratio
              LM_Stru%Rgl%Ratio = Ratio

              ! Update chi2
              LM_Stru%Chisq = LM_Stru%Chisq_og + &
                              LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

              ! Master
              if (pid.eq.0) then

                ! Verbose
                write(umsg,'(A)') ' - Updated regularization ratio:'
                call verboseI(3)
                write(umsg,'(A,es15.4,4(3x,A,es15.4))') &
                  '   New scaled penalty = ', &
                  LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio,  &
                  'New Chi2 (total) = ', &
                  LM_Stru%Chisq, &
                  'Chi2 (no regulatization) = ', &
                  LM_Stru%Chisq_og, &
                  'Ratio = ',LM_Stru%Rgl%Ratio, &
                  'Penalty = ',LM_Stru%Rgl%Penalty
                call verboseI(3)

              end if ! Master
            end if ! Lower ratio, but acceptable
          end if ! Regularizing and ratio larger than 1

          ! Update chi2
          Chisq_old = LM_Stru%Chisq

        end if ! Accepted step

        ! Predict lambda
        if (Input%l_Lam_track) &
          call update_lambda(Lam_track,Input%Lam_track, &
                             LM_Stru%Lambda)

      end do ! Iterations

      ! If last iteration was accepted and the type of error
      ! is Hessian
      if (LM_Stru%accepted.and.Input%Err_Type.ge.1) then

        ! Get Merit function
        LM_Stru%Chisq = LM_Stru%Chisq - &
                        LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

        ! Get Jacobian
        call Jacobian_Compute(Input,Atom,Atomb,Mol,Geom, &
                              GeomI,Flgsg,Frec,fudge,kurucz, &
                              MPID,Atmo,Bfield,Inf_Nodes,Sol, &
                              SolF,LM_Stru)

        ! If Hessian or worst error
        if (Input%Err_Type.eq.1.or.Input%Err_Type.eq.3) then

          ! Get Hessian
          call Hessian_Compute(Inf_Nodes%Nodes_Type, LM_Stru)

        end if

      ! Not accepted last step or error not from Hessian
      else

        ! If regularizing
        if (Inf_Nodes%Regul_Flag) then

          ! Recover merit function
          LM_Stru%Chisq = LM_Stru%Chisq - &
                          LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

          ! Recover Hessian
          LM_Stru%Hessian = LM_Stru%Hessian - &
                            LM_Stru%Rgl%Regul_H*LM_Stru%Rgl%Ratio
          ! Recover Jacobian
          LM_Stru%Jacfvec = LM_Stru%Jacfvec - &
                            LM_Stru%Rgl%Regul_F*LM_Stru%Rgl%Ratio

        end if ! Regularizing
      end if ! Accepted step and hessian error

      ! Compute error
      call Err(LM_Stru,Input,Inf_Stokes,Sol,Inf_Nodes)

      !
      ! Write the results into file
      !

      ! Control
      call control

      !
      ! Prepare output model atmosphere
      call setup_Atmo_ouinv(Atom,Atomb,Mol,Atmo,MPID,Input,fudge)


      ! Not aborted
      if (.not.laborted) then

        ! If LOS field, convert to vertical
        if (Inf_Nodes%Btype.eq.1) &
          call Bconversion(Atmo%nz,Inf_Stokes%mu, &
                           Inf_Stokes%azimuth, &
                           Bfield%Blos,Bfield%Bpos, &
                           Bfield%Azimuth,Bfield%Bstrength, &
                           Bfield%Btheta,Bfield%Bphi)

        ! If LOS velocity, convert to cartesian
        if (Inf_Nodes%vtype.eq.1) &
          call vconversion(Atmo%nz,Inf_Stokes%mu, &
                           Inf_Stokes%azimuth, &
                           Atmo%vlos,Atmo%vpos, &
                           Atmo%vphi,Atmo%vx, &
                           Atmo%vy,Atmo%vz)

        ! Master and saving
        if (pid.eq.0.and.saving) then

          ! Fake infinite loop
          do while (.True.)

            ! Write column to result file
            call Write_Result(Inf_Stokes,Sol,LM_Stru, &
                              Inf_Nodes,Atmo,Bfield,Input)
            if (laborted) exit

            !
            ! Write synthesis results

            ! Thermal not from a full inversion
            if (Inf_Nodes%Nodes_type.eq.0.and. &
                Input%Type_inversion.ne.3.and. &
                Input%Type_inversion.ne.4) then

              ! Write Stokes
              call writestkI(Input%folder,0,0,Frec%omega,GeomI, &
                             SolF%e_Stk_b(0,:,1,1),Input%lim_stk)
              if (laborted) exit


              ! Output contribution
              if (Input%out_contr) then
                call writectrI_inv(Input%folder, &
                                   SolF%e_Ctr_b(0,:,:,1,1), &
                                   Input%lim_ctr)
                if (laborted) exit
              end if

            ! Non-thermal
            else if (Inf_Nodes%Nodes_type.ne.0) then

              ! Write Stokes
              call writestk(Input%folder,0,0,Frec%omega,GeomI, &
                            SolF%e_Stk_b(:,:,1,1),Input%lim_stk)
              if (laborted) exit


              ! Output contribution
              if (Input%out_contr) then
                call writectr_inv(Input%folder, &
                                  SolF%e_Ctr_b(:,:,:,1,1), &
                                  Input%lim_ctr)
                if (laborted) exit
              end if

            end if

            ! Thermal not from a full inversion or non-thermal
            if ((Inf_Nodes%Nodes_type.eq.0.and. &
                 Input%Type_inversion.ne.3.and. &
                 Input%Type_inversion.ne.4).or. &
                (Inf_Nodes%Nodes_type.ne.0)) then

              ! Output tau1
              if (Input%out_tau1) &
                call writetau_inv(Input%folder, &
                                  SolF%e_tau1_b(:,1,1), &
                                  Input%lim_tau)
            end if

            exit

          end do ! Fake infinite loop

          ! Call control if did not fail
          if (.not.laborted) call control

        ! Slave if safe
        else if (.not.laborted) then

          ! Control
          call control

        end if ! Master saving
      end if ! Not aborting

      ! If thermal inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Deallocate intensity quantities
        deallocate(LM_Stru%ResidualI)
        deallocate(LM_Stru%WeightI)
        deallocate(LM_Stru%JacobianI)

      ! If magnetic
      else

        ! Deallocate polarization quantities
        deallocate(LM_Stru%Residual)
        deallocate(LM_Stru%Weight)
        deallocate(LM_Stru%Jacobian)

      end if ! Type of inversion

      ! Free memory
      deallocate(Solution,Errors,LM_Stru%Jacfvec,LM_Stru%diag)
      deallocate(LM_Stru%Hessian,LM_Stru%Hessian_og)
      deallocate(LM_Stru%Jacfvec_og)
      if (Input%l_Lam_track) deallocate(Lam_track)

      ! Free solutions
      call free_inv_solution(SolF)

      return

      end subroutine LMFIT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check if LM problem is converged\n
      !!    Input(Input_class): Structure with settings data\n
      !!     Chisq_old(double): Previous chi2\n
      !!         Chisq(double): Current chi2\n
      !!   Flag_Convg(logical): If converged
      subroutine Convergence_Check(Input,Chisq_old,Chisq, &
                                   Flag_Convg)

      ! IO
      type(Input_class), intent(in):: Input
      logical, intent(out):: Flag_Convg
      double precision, intent(in):: Chisq_old, Chisq

      ! Local
      double precision:: Chisq_fraction


      ! Get fraction
      Chisq_fraction = (Chisq_old-Chisq)/Chisq

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,es15.4)') '   Chi2 RC: ',Chisq_fraction
        call verboseI(3)

      end if

      ! Converged if the chi2 is absolutely smaller than
      ! threshold or changed too little, below the threshold
      Flag_Convg = (Chisq.lt.Input%Threshold_chisq).or. &
                   (Chisq_fraction.lt.Input%Chisq_fraction)

      return

      end subroutine Convergence_Check

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the error in the parameters\n
      !!     LM_Stru(LMFIT_class): Structure with Jacobian and other
      !!                           LM quantities\n
      !!       Input(Input_class): Structure with settings data\n
      !! Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                           data\n
      !!      Sol(Solution_class): Class with the data of the RT
      !!                           solution\n
      !!   Inf_Nodes(Nodes_class): Structure with nodes data
      subroutine Err(LM_Stru,Input,Inf_Stokes,Sol,Inf_Nodes)

      ! IO
      type(LMFIT_class), intent(in):: LM_Stru
      type(Input_class), intent(in):: Input
      type(Stokes_class), intent(in):: Inf_Stokes
      type(Solution_class), intent(in):: Sol
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local
      integer:: tmp, i, j, k

      double precision:: Error1, Error2, Num, Den
      double precision, dimension(:), allocatable:: Errors


      ! Allocate errors
      allocate(Errors(LM_Stru%Num))

      ! If error is default or hessian
      if (Input%Err_Type.eq.0.or.Input%Err_Type.eq.1) then

        ! For each element in the Hessian (all variables and nodes)
        do i=1,LM_Stru%Num

          ! Get error from Hessian
          Errors(i) = sqrt(LM_Stru%Chisq_og/ &
                           LM_Stru%Hessian_og(i,i)/ &
                           dble(LM_Stru%Num))
        end do

      ! If error is RF
      else if (Input%Err_Type.eq.2) then

        ! If only intensity
        if (Inf_Nodes%Nodes_Type.eq.0) then

          ! If sigma values
          if (Inf_Stokes%Sigma_flag) then

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! Add contribution to numerator
                Num = Num + max(LM_Stru%ResidualI(j)* &
                                LM_Stru%ResidualI(j), &
                                Inf_Stokes%Sigma_W(0,j)* &
                                Inf_Stokes%Sigma_W(0,j))* &
                      LM_Stru%WeightI(j)*LM_STru%WeightI(j)

                ! Add contribution to denominator
                Den = Den + LM_Stru%JacobianI(j,i)* &
                            LM_Stru%JacobianI(j,i)* &
                            LM_Stru%WeightI(j)* &
                            LM_Stru%WeightI(j)

              end do ! Wavelengths

              ! Get error from RF
              Errors(i) = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

            end do

          ! No sigmas
          else

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! Add contribution to numerator
                Num = Num + LM_Stru%ResidualI(j)* &
                            LM_Stru%ResidualI(j)* &
                            LM_Stru%WeightI(j)* &
                            LM_STru%WeightI(j)

                ! Add contribution to denominator
                Den = Den + LM_Stru%JacobianI(j,i)* &
                            LM_Stru%JacobianI(j,i)* &
                            LM_Stru%WeightI(j)* &
                            LM_Stru%WeightI(j)

              end do ! Wavelengths

              ! Get error from RF
              Errors(i) = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

            end do

          end if ! Sigma

        ! If polarization
        else

          ! If sigma
          if (Inf_Stokes%Sigma_flag) then

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! For all Stokes parameters
                do k=0,3

                  ! Add contribution to numerator
                  Num = Num + max(LM_Stru%Residual(k,j)* &
                                  LM_Stru%Residual(k,j), &
                                  Inf_Stokes%Sigma_W(k,j)* &
                                  Inf_Stokes%Sigma_W(k,j))* &
                        LM_Stru%Weight(k,j)*LM_STru%Weight(k,j)

                  ! Add contribution to denominator
                  Den = Den + LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Weight(k,j)* &
                              LM_Stru%Weight(k,j)

                end do ! Stokes parameters
              end do ! Wavelengths

              ! Get error from RF
              Errors(i) = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

            end do

          ! No sigma
          else

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! For all Stokes parameters
                do k=0,3

                  ! Add contribution to numerator
                  Num = Num + LM_Stru%Residual(k,j)* &
                              LM_Stru%Residual(k,j)* &
                              LM_Stru%Weight(k,j)* &
                              LM_STru%Weight(k,j)

                  ! Add contribution to denominator
                  Den = Den + LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Weight(k,j)* &
                              LM_Stru%Weight(k,j)

                end do ! Stokes parameters
              end do ! Wavelengths

              ! Get error from RF
              Errors(i) = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

            end do

          end if ! Sigma
        end if ! Intensity/polarization

      ! If error is Worst
      else if (Input%Err_Type.eq.3) then

        ! If only intensity
        if (Inf_Nodes%Nodes_Type.eq.0) then

          ! If sigma
          if (Inf_Stokes%Sigma_flag) then

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Get error from Hessian
              Error1 = sqrt(LM_Stru%Chisq_og/ &
                            LM_Stru%Hessian_og(i,i)/ &
                            dble(LM_Stru%Num))

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! Add contribution to numerator
                Num = Num + max(LM_Stru%ResidualI(j)* &
                                LM_Stru%ResidualI(j), &
                                Inf_Stokes%Sigma_W(0,j)* &
                                Inf_Stokes%Sigma_W(0,j))* &
                      LM_Stru%WeightI(j)*LM_STru%WeightI(j)

                ! Add contribution to denominator
                Den = Den + LM_Stru%JacobianI(j,i)* &
                            LM_Stru%JacobianI(j,i)* &
                            LM_Stru%WeightI(j)* &
                            LM_Stru%WeightI(j)

              end do ! Wavelengths

              ! Get error from RF
              Error2 = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

              ! Keep worst
              Errors(i) = max(Error1,Error2)

            end do

          ! No sigma
          else

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Get error from Hessian
              Error1 = sqrt(LM_Stru%Chisq_og/ &
                            LM_Stru%Hessian_og(i,i)/ &
                            dble(LM_Stru%Num))

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! Add contribution to numerator
                Num = Num + LM_Stru%ResidualI(j)* &
                            LM_Stru%ResidualI(j)* &
                            LM_Stru%WeightI(j)* &
                            LM_STru%WeightI(j)

                ! Add contribution to denominator
                Den = Den + LM_Stru%JacobianI(j,i)* &
                            LM_Stru%JacobianI(j,i)* &
                            LM_Stru%WeightI(j)* &
                            LM_Stru%WeightI(j)

              end do ! Wavelengths

              ! Get error from RF
              Error2 = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

              ! Keep worst
              Errors(i) = max(Error1,Error2)

            end do

          end if ! Sigma

        ! If polarization
        else

          ! If sigma
          if (Inf_Stokes%Sigma_flag) then

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Get error from Hessian
              Error1 = sqrt(LM_Stru%Chisq_og/ &
                            LM_Stru%Hessian_og(i,i)/ &
                            dble(LM_Stru%Num))

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! For all Stokes parameters
                do k=0,3

                  ! Add contribution to numerator
                  Num = Num + max(LM_Stru%Residual(k,j)* &
                                  LM_Stru%Residual(k,j), &
                                  Inf_Stokes%Sigma_W(k,j)* &
                                  Inf_Stokes%Sigma_W(k,j))* &
                        LM_Stru%Weight(k,j)*LM_Stru%Weight(k,j)

                  ! Add contribution to denominator
                  Den = Den + LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Weight(k,j)* &
                              LM_Stru%Weight(k,j)

                end do ! Stokes parameters
              end do ! Wavelengths

              ! Get error from RF
              Error2 = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

              ! Keep worst
              Errors(i) = max(Error1,Error2)

            end do

          ! No sigma
          else

            ! For all variables and nodes
            do i=1,LM_Stru%Num

              ! Get error from Hessian
              Error1 = sqrt(LM_Stru%Chisq_og/ &
                            LM_Stru%Hessian_og(i,i)/ &
                            dble(LM_Stru%Num))

              ! Initialize
              Num = 0d0
              Den = 0d0

              ! For all wavelengths
              do j=1,Sol%Num_wavelength

                ! For all Stokes parameters
                do k=0,3

                  ! Add contribution to numerator
                  Num = Num + LM_Stru%Residual(k,j)* &
                              LM_Stru%Residual(k,j)* &
                              LM_Stru%Weight(k,j)* &
                              LM_STru%Weight(k,j)

                  ! Add contribution to denominator
                  Den = Den + LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Jacobian(k,j,i)* &
                              LM_Stru%Weight(k,j)* &
                              LM_Stru%Weight(k,j)

                end do ! Stokes parameters
              end do ! Wavelengths

              ! Get error from RF
              Error2 = sqrt(2d0*Num/Den/dble(LM_Stru%Num))

              ! Keep worst
              Errors(i) = max(Error1,Error2)

            end do

          end if ! Sigma
        end if ! Intensity/polarization
      end if ! Type of error

      ! Initialize index shift
      tmp = 0

      ! For each variable index to consider
      do i = Inf_Nodes%Indx_b, Inf_Nodes%Indx_e

        ! If inverting that variable
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! Initialize error
          Inf_Nodes%Node(i)%Errors = 0d0

          ! Get error in the varying nodes from the Hessian errors
          Inf_Nodes%Node(i)%Errors(Inf_Nodes%Node_Vary(1,i): &
                                   Inf_Nodes%Node_Vary(2,i)) = &
            Inf_Nodes%Scal(i)* &
            Errors(tmp+1:tmp+Inf_Nodes%Num_Vary(i))

          ! Advance shift
          tmp = tmp + Inf_Nodes%Num_Vary(i)

        end if ! If inverting

      end do ! Variables

      ! Free memory
      deallocate(Errors)

      return

      end subroutine Err

!#####################################################################
!#####################################################################
!#####################################################################

      !> Do a LM step and get new Stokes profiles\n
      !!       LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                             other LM quantities
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!             Lambda(double): LM lambda factor\n
      !!        Solution(double(:)): SVD solution\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!         Input(Input_class): Structure with settings data\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution
      !!     SolF(Solution_F_class): Full solution of the RT problem
      subroutine Trial_Synthesis(LM_Stru,Inf_Nodes,Lambda,Solution, &
                                 Atmo,Bfield,Atom,Atomb,Mol,Geom, &
                                 GeomI,Flgsg,Frec,fudge,kurucz,MPID, &
                                 Input,Sol,SolF)

      ! IO
      type(LMFIT_class), intent(inout):: LM_Stru
      type(Nodes_class), intent(in):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
      type(Input_class), intent(inout):: Input
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      double precision, intent(in):: Lambda
      double precision, dimension(:), allocatable:: Solution

      ! Local
      type(Atmo_class):: Tmp_Atmo
      type(Bfield_class):: Tmp_Bfield
      type(Nodes_class):: Inf_Nodes_tmp

      integer:: i

      double precision, dimension(:), allocatable:: Jacfvec_new
      double precision, dimension(:,:), allocatable:: Hessian_new


      ! Copy current Nodes
      Inf_Nodes_tmp = Inf_Nodes

      ! Allocate new Hessian and Jacobian vector
      allocate(Hessian_new(LM_Stru%Num,LM_Stru%Num))
      allocate(Jacfvec_new(LM_Stru%Num))

      ! Copy current Hessian and Jacobian vectors
      Hessian_new = LM_Stru%Hessian
      Jacfvec_new = LM_Stru%Jacfvec

      ! Add lambda contribution to diagonal in Hessian matrix
      do i = 1, LM_Stru%Num
        Hessian_new(i,i) = Hessian_new(i,i) + Lambda*LM_Stru%diag(i)
      end do

      ! SVD solution
      call SVD_Solve(Hessian_new, Jacfvec_new, Solution, &
                     LM_Stru%Num, Inf_Nodes_tmp, Input%SVD_Type)
      if (laborted) return

      ! Modify the nodes with the LM solution
      call Nodes_Modify(Solution, Inf_Nodes_tmp)

      ! Call verbose
      if (pid.eq.0) &
        call Verbose_Model(Inf_Nodes_tmp,.True.)

      ! If thermal inversion only
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Get copy of atmosphere
        call cAtmo(Atmo,Tmp_Atmo)

        ! Generate new stratification
        call Intpol_Atmo(Inf_Nodes_tmp, Tmp_Atmo, Atom, Atomb, &
                         Mol, Input, fudge)

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Tmp_Atmo,Bfield,Input, &
                        Sol,SolF,.False.)

        ! Wipe Tmp_Atmo
        call free_Atmo(Tmp_Atmo,.True.)

      ! If magnetic inversion only
      else if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Get a copy of the magnetic field
        Tmp_Bfield = Bfield

        ! Generate new stratification
        call Intpol_Bfield(Inf_Nodes_tmp, Atmo, Tmp_Bfield)
        if (laborted) return

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Tmp_Bfield,Input, &
                        Sol,SolF,.False.)

      ! If inverting all
      else if (Inf_Nodes%Nodes_Type.eq.2) then

        ! Get a copy of atmosphere and field
        call cAtmo(Atmo,Tmp_Atmo)
        Tmp_Bfield = Bfield

        ! Get new stratification
        call Intpol_Atmo_all(Inf_Nodes_tmp, Tmp_Atmo, Tmp_Bfield, &
                             Atom, Atomb, Mol, Input, fudge)

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Tmp_Atmo,Tmp_Bfield, &
                        Input,Sol,SolF,.False.)

        ! Wipe Tmp_Atmo
        call free_Atmo(Tmp_Atmo,.True.)

      ! Error
      else

        ! Aborting
        umsg = 'The index of the node type is not correct'
        urou = 'Trial_Synthesis'
        call aborted
        return

      end if ! Type of inversion

      ! Deallocate auxiliar Hessian and Jacobian
      deallocate(Hessian_new,Jacfvec_new)

      ! Error
      if (laborted) return

      ! If regularizing, get regulatization
      IF (Inf_Nodes%Regul_Flag) &
        call Get_Regl_all(Inf_Nodes_tmp, .False., LM_Stru%Rgl)

      ! Factorize
      LM_Stru%Rgl%Penalty = LM_Stru%Rgl%Penalty*Input%Regul_factor

      return

      end subroutine Trial_Synthesis

!#####################################################################
!#####################################################################
!#####################################################################

      !> Optimize lambda parameter with backtracking algorithm\n
      !!       LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                             other LM quantities
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!    Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!     SolF(Solution_F_class): Class with the full RT problem
      !!                             solution\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!         Input(Input_class): Structure with settings data\n
      !!   Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                             data\n
      !!        Solution(double(:)): Current SVD solution\n
      !!        Stokes(double(:,:)): Current minimum chi2 Stokes
      !!                             parameters\n
      subroutine Backtracking(LM_Stru,Atom,Atomb,Mol,Geom,GeomI, &
                              Flgsg,Frec,fudge,kurucz,MPID,Atmo, &
                              Bfield,Sol,SolF,Inf_Nodes,Input, &
                              Inf_Stokes,Solution_Min,Stokes_Min)

      ! IO
      type(LMFIT_class), intent(inout):: LM_Stru
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Fctsg_class):: Flgsg
      type(Geometry_class):: GeomI, Geom
      type(Frequency_class):: Frec
      type(fudge_class):: fudge
      type(kurucz_class):: kurucz
      type(MPI_class):: MPID
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      double precision, dimension(:), allocatable:: Solution_Min
      double precision, dimension(:,:), allocatable:: Stokes_Min

      ! Local
      character(len=9):: ID

      logical:: Bracketed, up_first, converged, changed
      logical:: second_lap, early_exit

      integer, parameter:: Length = 10
      integer, parameter:: Max_fail = 3
      integer:: indx,jndx,kndx,Chisq_indx,nfail

      double precision:: Chisq_old, Chisq_best, daux
      double precision:: Lambda_TMP, Penalty_TMP
      double precision, dimension(:), allocatable:: Lambda_array
      double precision, dimension(:), allocatable:: Chisq_array
      double precision, dimension(:), allocatable:: Penalty_array
      double precision, dimension(:), allocatable:: Solution


      ! Create ID
      ID = 'MINIMUMSL'

      ! Create auxiliar arrays
      allocate(Lambda_array(Length+1),Chisq_array(Length+1))
      allocate(Solution(LM_Stru%Num),Penalty_array(Length+1))

      ! Master
      if (pid.eq.0) then

          ! Verbose
          write(umsg,'(A)') ' - Apply backtracking method'
          call verboseI(3)

      end if ! Master

      ! Initialize
      Bracketed = .False.
      early_exit = .False.
      Penalty_TMP = LM_Stru%Rgl%Penalty
      Chisq_indx = -1
      if (Input%l_Lam_track) then
        Lambda_array(1) = LM_Stru%lambda/LM_Stru%factoraccept
      else
        Lambda_array(1) = 0.1
      end if
      ! Comply with limits
      if (Lambda_array(1).lt.LM_Stru%Lambda_bounds(1)) &
          Lambda_array(1) = LM_Stru%Lambda_bounds(1)
      if (Lambda_array(1).gt.LM_Stru%Lambda_bounds(2)) &
          Lambda_array(1) = LM_Stru%Lambda_bounds(2)

      chisq_old = LM_Stru%Chisq
      up_first = .True.
      Chisq_best = chisq_old
      jndx = 1
      second_lap = .False.

      ! Fake loop
      do while (.True.)

        ! For the hard-coded length
        do indx=jndx,Length

          ! Master
          if (pid.eq.0) then

            ! Verbose
            if (gpid.eq.0) then
              umsg = '------------------------------------'
              call verbose
              umsg = '| Trial synthesis for backtracking |'
              call verbose
              umsg = '------------------------------------'
              call verbose
            end if
            umsg = ' - Trial synthesis for backtracking'
            call verboseI(3)

          end if ! Master

          ! Initialize failures
          nfail = 0

          ! Fake loop
          do while (.True.)

            ! Try new solution
            call Trial_Synthesis(LM_Stru, Inf_Nodes, &
                                 Lambda_array(indx), Solution, Atmo, &
                                 Bfield, Atom, Atomb, Mol, Geom, &
                                 GeomI, Flgsg, Frec, fudge, kurucz, &
                                 MPID, Input, Sol, SolF)

            ! If there was an issue with the trial
            if (laborted) then

              ! Verbose
              if (gpid.eq.0) then
                umsg = ' # Trial synthesis failed'
                call verbose
              end if
              umsg = ' # Trial synthesis failed'
              call verboseI(3)

              ! Cancel the failure
              laborted = .False.

              ! Skip this lambda
              Lambda_array(indx) = Lambda_array(indx)* &
                                   LM_Stru%factorreject

              ! Check no other equal
              do while (.True.)

                ! Flag
                changed = .False.

                ! Check
                do kndx=1,indx-1

                  ! Repeated
                  if (abs(Lambda_array(kndx) - &
                          Lambda_array(indx)).lt.1d-3) then
                    ! Flag, edit, and leave
                    changed = .True.
                    Lambda_array(indx) = Lambda_array(indx)* &
                                         LM_Stru%factorreject
                    exit
                  end if ! Repeated
                end do ! Check older lambda

                ! Loop again?
                if (changed) continue

                ! Leave
                exit

              end do ! Check repetition

              ! If Lambda went beyond limits, leave
              if (Lambda_array(indx).lt. &
                  LM_Stru%Lambda_bounds(1).or. &
                  Lambda_array(indx).gt. &
                  LM_Stru%Lambda_bounds(2)) then
                  
                ! Verbose
                if (gpid.eq.0) then
                  umsg = ' # Could not find new lambda to try with'
                  call verbose
                end if
                umsg = ' # Could not find new lambda to try with'
                call verboseI(3)

                ! Exit trials
                early_exit = .True.
                exit

              end if

              ! Add to count
              nfail = nfail + 1

              ! If beyond saving
              if (nfail.ge.Max_fail) return

              ! Verbose
              if (gpid.eq.0) then
                umsg = ' # Try new lambda value'
                call verbose
              end if
              umsg = ' # Try new lambda value'
              call verboseI(3)

              ! Try again
              continue

            end if ! Failed the Trial

            ! Exit the fake loop
            exit

          end do ! Fake loop

          ! If ealy exit, leave
          if (early_exit) exit

          ! Get merit function
          call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                              Inf_Nodes%Nodes_Type, LM_Stru)

          ! If regularizing
          if (Inf_Nodes%Regul_Flag) then

            ! Get actual chi2
            LM_Stru%Chisq = LM_Stru%Chisq + &
                            LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio

            ! Store penalty
            Penalty_array(indx) = LM_Stru%Rgl%Penalty

          ! No regularizing
          else

            ! Store penalty (zero)
            Penalty_array(indx) = 0d0

          end if ! Regularizing

          ! Save this chi2
          Chisq_array(indx) = LM_Stru%Chisq

          ! Master
          if (pid.eq.0) then

            ! Verbose
            write(umsg,'(A,i2)') ' +'
            call verboseI(3)
            write(umsg,'(A,i2)') ' + Trial ',indx
            call verboseI(3)
            write(umsg,'(4(3x,A,es15.4))') &
              'chi2 = ',LM_Stru%Chisq, &
              'lambda = ',Lambda_array(indx), &
              'Penalty = ',LM_Stru%Rgl%Penalty, &
              'Ratio = ',LM_Stru%Rgl%Ratio
            call verboseI(3)

          end if ! Master

          ! If larger than best
          if (Chisq_array(indx).gt.Chisq_best) then

            ! If first worse
            if (up_first) then

              ! Flag found
              up_first = .False.

              ! If first element
              if (indx.eq.1) then

                ! Set best
                call set_best(SolF,.False.,.False.)

                ! Save first as best
                Chisq_indx = indx
                Chisq_best = Chisq_array(indx)
                Stokes_Min = Sol%Stokes_out
                Solution_Min = Solution

              ! If the second element
              else if (indx.eq.2) then

                ! Swap values with last
                daux = Chisq_array(indx)
                Chisq_array(indx) = Chisq_array(indx-1)
                Chisq_array(indx-1) = daux
                daux = Lambda_array(indx)
                Lambda_array(indx) = Lambda_array(indx-1)
                Lambda_array(indx-1) = daux
                daux = Penalty_array(indx)
                Penalty_array(indx) = Penalty_array(indx-1)
                Penalty_array(indx-1) = daux

                ! Move best to second position
                Chisq_indx = indx

              ! Not any of the first two
              else

                ! It is bracketed, we can leave
                Bracketed = .True.
                exit

              end if ! Index

              ! Get next lambda
              Lambda_array(indx+1) = Lambda_array(indx)* &
                                     LM_Stru%factorreject

            ! Is the worse result
            else

             !! If index is beyond the second
             !if (indx.ge.3) then

             !    ! It is bracketed
             !    Bracketed = .True.
             !    exit

             !! Otherwise
             !else

             !    ! It is not bracketed
             !    Bracketed = .False.
             !    exit

             !end if ! Index

              ! If index is beyond the second and it was
              ! improvement at some point
              if (indx.ge.3.and.Chisq_indx.ne.1) then

                  ! It is bracketed
                  Bracketed = .True.
                  exit

              end if ! Index

              ! Get next lambda
              Lambda_array(indx+1) = Lambda_array(indx)* &
                                     LM_Stru%factorreject

            end if ! First worst result

          ! If is an improvement
          else

            ! Set best
            call set_best(SolF,.False.,.False.)

            ! Store new best Stokes and solution
            Stokes_Min = Sol%Stokes_out
            Solution_Min = Solution

            ! If first not found
            if (up_first) then

              ! Tried two steps
              if (indx.ge.3) then

                ! Check maximum relative change below
                ! a small number
                if ((Chisq_best-Chisq_array(indx))/ &
                    Chisq_best.lt.1d-2) then

                  ! Update best index and chi2
                  Chisq_indx = indx
                  Chisq_best = Chisq_array(indx)

                  ! Not bracketed
                  Bracketed = .False.
                  exit

                end if ! Small maximum relative change
              end if ! Tried already two steps

              ! Reduce lambda
              Lambda_array(indx+1) = Lambda_array(indx)/ &
                                     LM_Stru%factoraccept

            ! First already found
            else

              ! Increment
              Lambda_array(indx+1) = Lambda_array(indx)* &
                                     LM_Stru%factorreject

            end if ! Found the first worse result

            ! Update best index and chi2
            Chisq_indx = indx
            Chisq_best = Chisq_array(indx)

          end if ! Improvement or not

          ! If Lambda went beyond limits, leave
          if (Lambda_array(indx+1).lt.LM_Stru%Lambda_bounds(1).or. &
              Lambda_array(indx+1).gt.LM_Stru%Lambda_bounds(2)) exit

        end do ! Try up to Length steps


        ! If variation is bracketed
        if (Bracketed) then

          ! Master
          if (pid.eq.0) then
            umsg = ' + Lambda bracketed'
            call verboseI(3)
          end if

          ! Parabolic interpolation
          call Parabolic(Lambda_array, Chisq_array, Chisq_indx, &
                         Lambda_TMP)

          ! Master
          if (pid.eq.0) then

            ! Verbose
            if (gpid.eq.0) then
              umsg = '---------------------------------------'// &
                     '----------'
              call verbose
              umsg = '| Trial synthesis for interpolated '// &
                     'backtracking |'
              call verbose
              umsg = '----------------------------------------'// &
                     '---------'
              call verbose
            end if
            umsg = ' - Trial synthesis for interpolated backtracking'
            call verboseI(3)

          end if ! Master

          ! Try new solution
          call Trial_Synthesis(LM_Stru, Inf_Nodes, Lambda_TMP, &
                               Solution, Atmo, Bfield, Atom, Atomb, &
                               Mol, Geom, GeomI, Flgsg, Frec, fudge, &
                               kurucz, MPID, Input, Sol, SolF)

          ! If failed
          if (laborted) then

            ! Cancel the failure
            laborted = .False.

            ! Verbose
            if (gpid.eq.0) then
              umsg = ' # Trial synthesis failed'
              call verbose
            end if
            umsg = ' # Trial synthesis failed'
            call verboseI(3)

            ! Exagerate new chi2
            LM_Stru%Chisq = Chisq_array(Chisq_indx)*1d4

          ! Did not fail
          else

            ! Get merit function
            call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                                Inf_Nodes%Nodes_Type, LM_Stru)

            ! Regularize
            if (Inf_Nodes%Regul_Flag) &
              LM_Stru%Chisq = LM_Stru%Chisq + LM_Stru%Rgl%Penalty* &
                                              LM_Stru%Rgl%Ratio

          end if

          ! If new chi2 is better than the best in backtracking
          if (LM_Stru%Chisq.lt.Chisq_array(Chisq_indx)) then

            ! Set best
            call set_best(SolF,.False.,.False.)

            ! New best
            Stokes_Min = Sol%Stokes_out
            Solution_Min = Solution
            LM_Stru%Lambda = Lambda_TMP

            ! Master
            if (pid.eq.0) then
              write(umsg,'(A)') &
                 ' - Improved chi2 with parabolic interpolation'
              call verboseI(3)
              write(umsg,'(3(A,es15.4))') &
                '   chi2: ',LM_Stru%chisq,'   <', &
                Chisq_array(Chisq_indx), &
                '   lambda:',LM_Stru%Lambda
              call verboseI(3)
            end if

          ! If it did not improve
          else

            ! Keep best from steps before
            LM_Stru%Lambda = Lambda_array(Chisq_indx)
            LM_Stru%Rgl%Penalty = Penalty_array(Chisq_indx)
            LM_Stru%Chisq = Chisq_array(Chisq_indx)

            ! Master
            if (pid.eq.0) then
              write(umsg,'(A)') &
                 ' - Chi2 did not improve with parabolic '// &
                 'interpolation'
              call verboseI(3)
              write(umsg,'(A,es15.4)') &
                '   lambda:',LM_Stru%Lambda
              call verboseI(3)
            end if

          end if ! The interpolation improved the result

        ! Variation not bracketed
        else

          ! Master
          if (pid.eq.0) then
            umsg = ' + Lambda not bracketed'
            call verboseI(3)
          end if

          ! If could not find
          if (Chisq_indx.lt.1) then

            ! Just take the last
            LM_Stru%Rgl%Penalty = Penalty_TMP
            LM_Stru%Chisq = chisq_old

          ! If could find
          else

            ! Just take the best
            LM_Stru%Lambda = Lambda_array(Chisq_indx)
            LM_Stru%Rgl%Penalty = Penalty_array(Chisq_indx)
            LM_Stru%Chisq = Chisq_array(Chisq_indx)

          end if ! Could find something

        end if ! Bracketed variation

        ! If desperate
        if (Input%LM_Back_Mode.eq.1.and..not.second_lap) then

          ! Check convergence here
          call Convergence_Check(Input, Chisq_old, &
                                 LM_Stru%Chisq, converged)

          ! If Chi2 did not improve or the improvement is
          ! too minor
          if (LM_Stru%Chisq.ge.Chisq_old.or. &
              (LM_Stru%Chisq.lt.Chisq_old.and.converged)) then

            ! Get last index
            if (indx.gt.Length) indx = Length

            ! If the minimum lambda was rather big or the
            ! maximum rather small, restart with the opposite
            if (minval(Lambda_array(1:indx)).gt. &
                Input%LM_lam_big_test.or. &
                maxval(Lambda_array(1:indx)).lt. &
                Input%LM_lam_small_test) then

              ! Get best
              Penalty_array(1) = Penalty_array(Chisq_indx)
              Chisq_array(1) = Chisq_array(Chisq_indx)
              Lambda_array(1) = Lambda_array(Chisq_indx)

              ! Restart
              Penalty_array(2:indx) = 0d0
              Chisq_array(2:indx) = 0d0
              Lambda_array(3:indx) = 0d0
              up_first = .True.
              jndx = 2

              ! If the minimum lambda was rather big, restart
              ! with something smaller
              if (minval(Lambda_array(1:indx)).gt. &
                  Input%LM_lam_big_test) then

                ! Next lambda
                Lambda_array(2) = Input%LM_lam_small_prove

                ! Master
                if (pid.eq.0) then

                  ! Verbose
                  write(umsg,'(A,i2)') ' +'
                  call verboseI(3)
                  umsg = ' + Failed to improve the fit. '// &
                         'The lambda values may be too large, '// &
                         'restarting backtracking with smaller '// &
                         'values'
                  call verboseI(3)
                  write(umsg,'(A,i2)') ' + Trial ',1
                  call verboseI(3)
                  write(umsg,'(4(3x,A,es15.4))') &
                    'chi2 = ',Chisq_array(1), &
                    'lambda = ',Lambda_array(1), &
                    'Penalty = ',Penalty_array(1)
                  call verboseI(3)

                end if ! Master

                ! Restart
                second_lap = .True.
                cycle

              end if ! Minimum lambda rather big

              ! If the maximum lambda was rather small, restart
              ! with something larger
              if (maxval(Lambda_array(1:indx)).lt. &
                  Input%LM_lam_small_test) then

                ! Next lambda
                Lambda_array(2) = Input%LM_lam_big_prove

                ! Master
                if (pid.eq.0) then

                  ! Verbose
                  write(umsg,'(A,i2)') ' +'
                  call verboseI(3)
                  umsg = ' + Failed to improve the fit. '// &
                         'The lambda values may be too small, '// &
                         'restarting backtracking with larger '// &
                         'values'
                  call verboseI(3)
                  write(umsg,'(A,i2)') ' + Trial ',1
                  call verboseI(3)
                  write(umsg,'(4(3x,A,es15.4))') &
                    'chi2 = ',Chisq_array(1), &
                    'lambda = ',Lambda_array(1), &
                    'Penalty = ',Penalty_array(1)
                  call verboseI(3)

                end if ! Master

                ! Restart
                second_lap = .True.
                cycle

              end if ! Maximum lambda rather small
            end if ! Single regime lambda
          end if ! Failed to improve significantly
        end if ! Desperate to reduce chi2

        ! Break the fake loop
        exit

      end do

      ! Free memory
      deallocate(Lambda_array,Chisq_array,Solution,Penalty_array)

      return

      end subroutine Backtracking

!#####################################################################
!#####################################################################
!#####################################################################

      !> Propose a lambda factor for the LM\n
      !! LM_Stru(LMFIT_class): Structure with Jacobian and other LM
      !!                       quantities\n
      !!       Lambda(double): The Lambda factor
      subroutine Lambda_propose_fix(LM_Stru,Lambda)

      ! IO
      type(LMFIT_class), intent(in):: LM_Stru
      double precision, intent(inout):: Lambda


      ! If it was accepted
      if(LM_Stru%accepted) then

        ! Reduce
        Lambda = Lambda/LM_Stru%factoraccept

      ! If it was rejected
      else

        ! Increase
        Lambda = Lambda*LM_Stru%factorreject

      endif

      return

      end subroutine Lambda_propose_fix

!#####################################################################
!#####################################################################
!#####################################################################

      !> Predict the next value for lambda\n
      !!        iter(integer): Current iteration\n
      !!       track(dble(:)): Lambda history\n
      !!      ntrack(integer): Size of information stored\n
      !! LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                       other LM quantities
      subroutine predict_lambda(iter,track,ntrack,LM_Stru)

      ! IO
      integer, intent(in):: iter,ntrack
      double precision, dimension(:), intent(in):: track
      type(LMFIT_class), intent(inout):: LM_Stru

      ! Local
      double precision:: d1,d2,d3
      double precision, dimension(ntrack):: x,a,b,c


      ! If first iteration, skip
      if (iter.eq.1) return

      ! Cases
      select case (ntrack)

        ! Constant
        case (1)

          ! Get new lambda
          LM_Stru%Lambda = track(1)

          ! Master
          if (pid.eq.0) then

            ! Verbose
            write(umsg,'(A,1x,es15.4)')  &
              ' - Tracking lambda:',track(1)

            if (gpid.eq.0) then
              call verboseI(0)
              call verboseI(4)
            else
              call verboseI(0)
              if (vlevel.eq.0) call verboseI(3)
            end if
          end if

          ! Sanity check
          if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1)) &
            LM_Stru%Lambda = LM_Stru%Lambda_bounds(1)
          if (LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) &
            LM_Stru%Lambda = LM_Stru%Lambda_bounds(2)

        ! Linear
        case (2)

          ! One one iteration
          if (iter.lt.3) then

            ! Get new lambda
            LM_Stru%Lambda = track(2)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,1x,es15.4)')  &
                ' - Tracking lambda:',track(2)

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if

          ! At least two iterations
          else

            ! Get new lambda
            LM_Stru%Lambda = 2*track(2) - track(1)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,2(1x,es15.4),A,es15.4)')  &
                ' - Tracking lambda:',track(1),track(2), &
                ' -> ',LM_Stru%Lambda

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if

          end if

          ! Sanity check
          if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1)) &
            LM_Stru%Lambda = LM_Stru%Lambda_bounds(1)
          if (LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) &
            LM_Stru%Lambda = LM_Stru%Lambda_bounds(2)

        ! Parabolic
        case (3)

          ! One one iteration
          if (iter.lt.3) then

            ! Get new lambda
            LM_Stru%Lambda = track(3)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,1x,es15.4)')  &
                ' - Tracking lambda:',track(3)

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if

            ! Sanity check
            if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1).or.&
                LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) &
                LM_Stru%Lambda = track(3)

          ! Only two iterations
          else if (iter.lt.4) then

            ! Get new lambda
            LM_Stru%Lambda = 2*track(3) - track(2)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,2(1x,es15.4),A,es15.4)')  &
                ' - Tracking lambda:',track(2),track(3), &
                ' -> ',LM_Stru%Lambda

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if

            ! Sanity check
            if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1).or.&
                LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) &
                LM_Stru%Lambda = track(3)

          ! At least three iterations
          else

            ! Prepare spline
            x = (/ 0d0,1d0,2d0 /)
            call spline(x,track,a,b,c,ntrack)

            ! Get new lambda
            LM_Stru%Lambda = track(3) + a(3) + b(3) + c(3)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,3(1x,es15.4),A,es15.4)')  &
                ' - Tracking lambda:',track(1),track(2), &
                track(3),' -> ',LM_Stru%Lambda

              if (gpid.eq.0) then
                call verboseI(0)
                call verboseI(4)
              else
                call verboseI(0)
                if (vlevel.eq.0) call verboseI(3)
              end if
            end if

            ! Sanity check
            if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1).or.&
                LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) then
                LM_Stru%Lambda = track(3)
            else

              ! Check derivatives
              d1 = track(2) - track(1)
              d2 = track(3) - track(2)
              d3 = LM_Stru%Lambda - track(2)

              ! If derivative changed sign in a monotonic case
              if (d3*d2.lt.0d0.and.d3*d1.lt.0d0) then

                ! If derivative became negative
                if (d3.lt.0d0) then

                  ! Check minimum
                  if (LM_Stru%Lambda.lt.track(1)) &
                    LM_Stru%Lambda = track(1)

                ! If  derivative became positive
                else

                  ! Check maximum
                  if (LM_Stru%Lambda.gt.track(1)) &
                    LM_Stru%Lambda = track(1)

                end if ! Derivative sign
              end if ! Inversion of monotonicity
            end if ! Sanity
          end if ! Number of iterations

      end select ! Type of prediction

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,1x,es15.4)')  &
          ' - New lambda:',LM_Stru%Lambda

        if (gpid.eq.0) then
          call verboseI(0)
          call verboseI(4)
        else
          call verboseI(0)
          if (vlevel.eq.0) call verboseI(3)
        end if
      end if

      return

      end subroutine predict_lambda

!#####################################################################
!#####################################################################
!#####################################################################

      !> Update the tracking of lambda\n
      !!       track(dble(:)): Lambda history\n
      !!      ntrack(integer): Size of information stored\n
      !!       lambda(double): The Lambda factor
      subroutine update_lambda(track,ntrack,lambda)

      ! IO
      integer, intent(in):: ntrack
      double precision, intent(in):: Lambda
      double precision, dimension(:), intent(inout):: track

      ! Local
      integer:: i

      ! Slide
      do i=1,ntrack-1
        track(i) = track(i+1)
      end do

      ! Add
      track(ntrack) = lambda

      return

      end subroutine update_lambda

!#####################################################################
!#####################################################################
!#####################################################################

      !> Modify the node values according to the SVD solution\n
      !!    Solution(double(:)): Last solution of the Hessian\n
      !! Inf_Nodes(Nodes_class): Structure with nodes data
      subroutine Nodes_Modify(Solution,Inf_Nodes)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      double precision, dimension(:), allocatable:: Solution

      ! Local
      integer:: tmp, i


      ! Initialize running index
      tmp = 0

      ! For each parameter to modify
      do i=Inf_Nodes%Indx_b,Inf_Nodes%Indx_e

        ! If inverting the variable
        if (Inf_Nodes%Nodes_Flags(i)) then

          ! If nodes are values
          if (Inf_Nodes%Node_Type(i).le.3) then

            ! Change the values in the varying nodes
            Inf_Nodes%Node(i)%Var(Inf_Nodes%Node_Vary(1, i): &
                                  Inf_Nodes%Node_Vary(2, i)) = &
               Inf_Nodes%Node(i)%Var(Inf_Nodes%Node_Vary(1, i): &
                                     Inf_Nodes%Node_Vary(2, i)) + &
               solution(tmp+1:tmp+Inf_Nodes%Num_Vary(i))* &
               Inf_Nodes%Scal(i)

            ! Check boundaries
            call CheckBounds(Inf_Nodes%Node(i), &
                             Inf_Nodes%Num_Nodes(i))

          ! If nodes are corrections
          else

            ! Just save solution
            Inf_Nodes%Node(i)%Var(Inf_Nodes%Node_Vary(1, i): &
                                  Inf_Nodes%Node_Vary(2, i)) = &
                             Inf_Nodes%Scal(i)* &
                             solution(tmp+1:tmp+Inf_Nodes%Num_Vary(i))

          end if ! Values or corrections

          ! Shift index counter
          tmp = tmp + Inf_Nodes%Num_Vary(i)

        end if ! Inverting variable

      end do ! Variables

      return

      end subroutine Nodes_Modify

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check limits of lambda factor\n
      !!           Lambda(double): Value of lambda factor\n
      !! Lambda_limits(double(2)): Boundary limits for lambda factor
      subroutine CheckLambda(Lambda,Lambda_limits)

      ! IO
      double precision, intent(inout):: Lambda
      double precision, dimension(2), intent(in):: Lambda_limits

      ! If above limits
      if (Lambda.ge.Lambda_limits(2)) then

        ! Force limit
        Lambda = Lambda_limits(2)

      ! If below limits
      else if (Lambda.le.Lambda_limits(1)) then

        ! Force limit
        Lambda = Lambda_limits(1)

      end if

      return

      end subroutine CheckLambda

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set best solution\n
      !!  Sol(Solution_F_class): Class with the data of the RT
      !!                         solution\n
      !!          best(logical): If really the best or just in
      !!                         backtracking\n
      !!          copy(logical): If getting from backtracking
      subroutine set_best(Sol,best,copy)

      ! I/O
      type(Solution_F_class), intent(inout):: Sol
      logical, intent(in):: best, copy

      ! Slaves, leave
      if (pid.gt.0) return

      ! Truly the best
      if (best) then

        ! If copying
        if (copy) then

          ! Copy what is allocated
          if (allocated(Sol%i_J00_t)) Sol%i_J00_b = Sol%i_J00_t
          if (allocated(Sol%i_J00C_t)) Sol%i_J00C_b = Sol%i_J00C_t
          if (allocated(Sol%e_tau1_t)) Sol%e_tau1_b = Sol%e_tau1_t
          if (allocated(Sol%i_J00P_t)) Sol%i_J00P_b = Sol%i_J00P_t
          if (allocated(Sol%e_Stk_t)) Sol%e_Stk_b = Sol%e_Stk_t
          if (allocated(Sol%e_Ctr_t)) Sol%e_Ctr_b = Sol%e_Ctr_t
          if (allocated(Sol%i_StkI_t)) Sol%i_StkI_b = Sol%i_StkI_t
          if (allocated(Sol%i_Stk_t)) Sol%i_Stk_b = Sol%i_Stk_t
          if (allocated(Sol%i_JKQ_t)) Sol%i_JKQ_b = Sol%i_JKQ_t
          if (allocated(Sol%i_JKQC_t)) Sol%i_JKQS_b = Sol%i_JKQS_t
          if (allocated(Sol%i_JKQC_t)) Sol%i_JKQC_b = Sol%i_JKQC_t
          if (allocated(Sol%i_rhoes_t)) Sol%i_rhoes_b = Sol%i_rhoes_t

        ! Not copying
        else

          ! Copy what is allocated
          if (allocated(Sol%i_J00)) Sol%i_J00_b = Sol%i_J00
          if (allocated(Sol%i_J00C)) Sol%i_J00C_b = Sol%i_J00C
          if (allocated(Sol%e_tau1)) Sol%e_tau1_b = Sol%e_tau1
          if (allocated(Sol%i_J00P)) Sol%i_J00P_b = Sol%i_J00P
          if (allocated(Sol%e_Stk)) Sol%e_Stk_b = Sol%e_Stk
          if (allocated(Sol%e_Ctr)) Sol%e_Ctr_b = Sol%e_Ctr
          if (allocated(Sol%i_StkI)) Sol%i_StkI_b = Sol%i_StkI
          if (allocated(Sol%i_Stk)) Sol%i_Stk_b = Sol%i_Stk
          if (allocated(Sol%i_JKQ)) Sol%i_JKQ_b = Sol%i_JKQ
          if (allocated(Sol%i_JKQC)) Sol%i_JKQS_b = Sol%i_JKQS
          if (allocated(Sol%i_JKQC)) Sol%i_JKQC_b = Sol%i_JKQC
          if (allocated(Sol%i_rhoes_t)) Sol%i_rhoes_b = Sol%i_rhoes

        end if ! Copying from backtrace?

      ! Provisional best
      else

        ! Copy what is allocated
        if (allocated(Sol%i_J00)) Sol%i_J00_t = Sol%i_J00
        if (allocated(Sol%i_J00C)) Sol%i_J00C_t = Sol%i_J00C
        if (allocated(Sol%e_tau1)) Sol%e_tau1_t = Sol%e_tau1
        if (allocated(Sol%i_J00P)) Sol%i_J00P_t = Sol%i_J00P
        if (allocated(Sol%e_Stk)) Sol%e_Stk_t = Sol%e_Stk
        if (allocated(Sol%e_Ctr)) Sol%e_Ctr_t = Sol%e_Ctr
        if (allocated(Sol%i_StkI)) Sol%i_StkI_t = Sol%i_StkI
        if (allocated(Sol%i_Stk)) Sol%i_Stk_t = Sol%i_Stk
        if (allocated(Sol%i_JKQ)) Sol%i_JKQ_t = Sol%i_JKQ
        if (allocated(Sol%i_JKQC)) Sol%i_JKQS_t = Sol%i_JKQS
        if (allocated(Sol%i_JKQC)) Sol%i_JKQC_t = Sol%i_JKQC
        if (allocated(Sol%i_rhoes)) Sol%i_rhoes_t = Sol%i_rhoes

      end if ! Truly the best

      end subroutine set_best

!#####################################################################
!#####################################################################
!#####################################################################

      end module lmfit_mod
