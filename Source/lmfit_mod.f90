      !> Levenberg-Marquardt fit
      module lmfit_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     22/03/2023
!  Last version:
!     18/12/2025 V4.0.9
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/12/2025:    V4.0.9 - Added reduced mode, a last chance
!                             iteration which modifies the weights
!                             to try compensating for a bad initial
!                             weighting. It is optional (TdPA)
!                           - NaN checks now use ieee (TdPA)
!                           - Added a tracker for the lambda
!                             prediction to make it not depend on the
!                             particular iteration index (TdPA)
!                           - Added verbosity to highlight the change
!                             of chi^2 after the update of the
!                             regularization weight (TdPA)
!                           - Modified calls of routines that have
!                             changed their argument list (TdPA)
!                           - Added logic to stop the Backtracking
!                             when in gradient descend regime with
!                             little hope to improve the fit (TdPA)
!                           - In backtracking, if the chi^2 is a
!                             disaster, the lambda parameter is
!                             enhanced more than once (TdPA)
!                           - Added verbosity when leaving the
!                             backtracking (TdPA)
!                           - Bugfix: there was a typo in which some
!                             cycle instances were written as
!                             continue instead (TdPA)
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
!  LMFIT
!    Levenberg-Marquardt fit of Stokes parameters
!
!  write_result_inv
!    Prepare the model atmosphere and Stokes profiles to write the
!  current result in the output and call such writing
!
!  Convergence_Check
!    Check the convergence of the merit function
!
!  Err
!    Calculate an estimation of the error in the inverted parameters
!
!  Trial_Synthesis
!    Give a Levenberg-Marquardt step and compute emergent profiles
!
!  Backtracking
!    Optimize the Levenberg-Marquardt lambda parameter with the
!  backtracking algorithm
!
!  Lambda_propose_fix
!    Propose a new Levenberg-Marquardt lambda parameter
!
!  predict_lambda
!    Propose the value of the Levenberg-Marquardt lambda parameter for
!  the next iteration
!
!  update_lambda
!    Update history of Levenberg-Marquardt lambda parameter
!
!  Nodes_Modify
!    Modify the parameter node values according to the SVD solution
!
!  CheckLambda
!    Enforce limits of Levenberg-Marquardt lambda parameter
!
!  set_best
!    Set a self-consistent solution as the best
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

      !> Levenberg-Marquardt fit of Stokes parameters\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!            imask(integer): Indicate if this pixel is masked
      !!                            when restarting the inversion\n
      !!           saving(logical): If the result is to be stored
      subroutine LMFIT(Atom,Atomb,Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                       kurucz,MPID,Atmo,Bfield,Input,Inf_Stokes, &
                       Inf_Nodes,Sol,LM_Stru,imask,saving)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
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

      type(Solution_F_class):: SolF

      logical:: Flag_Convg,Flag_Jac,Flag_RDmode,skip_Jac

      integer:: indx_iter,indx_rej,i,Num_Broyden,max_iters

      double precision:: Chisq_old,Ratio,Chisq_RD0,frac
      double precision, dimension(:), allocatable:: Lam_track
      double precision, dimension(:), allocatable:: Solution,Errors
      double precision, dimension(:,:), allocatable:: Stokes_Min
      double precision, dimension(:,:), allocatable:: Stokes_best


      ! Count innate memory of local structure
      SRAMc = SRAMc + 1d-6*sizeof(SolF)

      ! Check allocations (in case of previous failure)
      if (allocated(LM_Stru%Hessian)) then

        ! Free residual for intensity if allocated
        if (allocated(LM_Stru%ResidualI)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%ResidualI)
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%WeightI)
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%JacobianI)
          deallocate(LM_Stru%ResidualI)
          deallocate(LM_Stru%WeightI)
          deallocate(LM_Stru%JacobianI)
        end if
        if (allocated(LM_Stru%WeightI_mod)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%WeightI_mod)
          deallocate(LM_Stru%WeightI_mod)
        end if

        ! Free residual for polarization if allocated
        if (allocated(LM_Stru%Residual)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Residual)
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Weight)
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacobian)
          deallocate(LM_Stru%Residual)
          deallocate(LM_Stru%Weight)
          deallocate(LM_Stru%Jacobian)
        end if
        if (allocated(LM_Stru%Weight_mod)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Weight_mod)
          deallocate(LM_Stru%Weight_mod)
        end if

        ! Free memory
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacfvec)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%diag)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Hessian)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Hessian_og)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacfvec_og)
        deallocate(LM_Stru%Jacfvec,LM_Stru%diag)
        deallocate(LM_Stru%Hessian,LM_Stru%Hessian_og)
        deallocate(LM_Stru%Jacfvec_og)

        ! Free solutions
        call free_inv_solution(SolF)

      end if ! Previous failure

      !
      ! If Thermal inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Allocate arrays for just intensity
        allocate(LM_Stru%ResidualI(Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%WeightI(Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%JacobianI(Sol%Num_Wavelength,LM_Stru%Num))
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%ResidualI)
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%WeightI)
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%JacobianI)

      ! Magnetic or full inversion
      else

        ! Allocate full Stokes arrays
        allocate(LM_Stru%Residual(0:3,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%Weight(0:3,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%Jacobian(0:3,Sol%Num_Wavelength,LM_Stru%Num))
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Residual)
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Weight)
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Jacobian)

      end if ! Type of inversion

      !
      ! Determine if trials change thermodynamics
      !


      ! If just magnetic
      if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Fixed thermodynamics
        Sol%fix_th = .True.
        Sol%fix_tp = .True.

      ! Not just magnetic
      else

        ! Initialize as fixed thermodynamics
        Sol%fix_th = .True.

        ! Check all thermal variables
        do i=Inf_Nodes%index_T,Inf_Nodes%index_J21R-1

          ! Check if not fixed
          if (Inf_Nodes%Nodes_Flags(i)) then

            ! Flag not fixed and leave
            Sol%fix_th = .False.
            exit

          end if ! Non-fixed

        end do ! Thermal variables

        ! Determine if trials change temperature or gas pressure
        Sol%fix_tp = &
                (.not.Inf_Nodes%Nodes_Flags(Inf_Nodes%index_T)).and. &
                (.not.Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg))

      end if ! Type of inversion

      ! Initialize ratio for the first iteration
      Ratio = 0.1d0

      ! If ratio larger than regularitaion limit, reduce it now
      if (Ratio.gt.Input%Regul_Limit) Ratio = Input%Regul_Limit

      ! Initialize as false
      LM_Stru%Flag_weight = .False.

      ! Flag solution buffer not initialized
      SolF%no_initialized = .True.

      ! Flag not in reduced mode
      Flag_RDmode = .False.

      ! Do not skip computing Jacobian
      skip_Jac = .False.

      ! Allocate solution and Hessian quantities
      allocate(Solution(LM_Stru%Num),Errors(LM_Stru%Num))
      allocate(LM_Stru%Jacfvec(LM_Stru%Num),LM_Stru%diag(LM_Stru%Num))
      allocate(LM_Stru%Hessian(LM_Stru%Num,LM_Stru%Num))
      allocate(LM_Stru%Hessian_og(LM_Stru%Num,LM_Stru%Num))
      allocate(LM_Stru%Jacfvec_og(LM_Stru%Num))
      MRAMc = MRAMc + 1d-6*sizeof(Solution)
      MRAMc = MRAMc + 1d-6*sizeof(Errors)
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Jacfvec)
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%diag)
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Hessian)
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Hessian_og)
      MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Jacfvec_og)

      ! Allocate lambda array
      if (Input%l_Lam_track) then
        allocate(Lam_track(Input%Lam_track))
        MRAMc = MRAMc + 1d-6*sizeof(Lam_track)
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

      end if ! Master

      ! Get first synthesis
      call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                      kurucz,MPID,Atmo,Bfield,Input,Sol,SolF,0)

      ! Set the first solution as current best
      call set_best(SolF,.True.,.False.)

      ! Compute initial merit function
      call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                          Inf_Nodes%Nodes_Type,LM_Stru)

      ! Check NaN
      if (ieee_is_nan(LM_Stru%Chisq)) then

        ! Issue error
        umsg = 'The first value of the merit function (without '// &
               'penalties) is NaN. Check that there are no '// &
               'zeros in the "sigma" values in the data'
        urou = 'LMFIT'
        call aborted
        goto 2000

      end if ! NaN chi^2

      ! If regularizing
      if (Inf_Nodes%Regul_Flag) then

        ! Allocate space for regularization
        allocate(LM_Stru%Rgl%Regul_H(LM_Stru%Num,LM_Stru%Num))
        allocate(LM_Stru%Rgl%Regul_F(LM_Stru%Num))
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Rgl%Regul_H)
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Rgl%Regul_F)

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
        if (vlevel.gt.0) call verboseI(3)

        ! Verbose merit function
        write(umsg,'(A,i4,3x,A,es15.4)')  &
           ' * Iteration = ',0, &
              'Chi2 = ',LM_Stru%Chisq

        ! Verbose depending on the MPI regime
        if (gpid.eq.0) then
          call verboseI(0)
          if (vlevel.gt.0) call verboseI(3)
        else
          call verboseI(0)
          if (vlevel.gt.0) call verboseI(3)
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
      LM_Stru%nLambda = 0
      LM_Stru%Flag_weight = .False.
      LM_Stru%factorreject = Input%factorreject
      LM_Stru%factoraccept = Input%factoraccept
      LM_Stru%accepted = .True.
      LM_Stru%perc = 0.25d0
      LM_Stru%ff_max = 25d0
      LM_Stru%ff_contr = 0.35d0
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

        ! Check if writing
        if (Input%storeinv.and.indx_iter.gt.1) then

          ! Check module
          if (mod(indx_iter-1,Input%storeinv_step).eq.0) then

            ! Write to files
            call write_result_inv(Atom,Atomb,Mol,GeomI,Geom,Frec, &
                                  fudge,Atmo,Bfield,Input, &
                                  Inf_Stokes,Inf_Nodes,Sol,SolF, &
                                  LM_Stru,saving,.False.)

          end if ! Is time to produce a partial result file
        end if ! If producing not finished result files

        ! Predict lambda for LM
        if (Input%l_Lam_track) &
          call predict_lambda(Lam_track, \
                              Input%Lam_track,LM_Stru)

        ! If after the first iteration the accepted flag is on,
        ! we are doing Broyden and we have not done three already
        if (Input%Broyden.and.LM_Stru%accepted.and. &
            Num_Broyden.lt.3.and.indx_iter.gt.1) then

          ! If skipping
          if (skip_Jac) then

            ! Flag false
            skip_Jac = .False.

          ! Normal
          else

            ! Call Broyden
            call Broyden_Rank1(Stokes_best, Stokes_Min, &
                               Inf_Stokes%Num_Wavelength, &
                               Solution, Inf_Nodes, LM_Stru)
            ! Deflag Jacobian
            Flag_Jac = .False.

          end if

        ! First iteration, not accepted LM, not doing Broyden,
        ! or already three Broyden
        else

          ! If skipping Jacobian (last iteration failed and
          ! now reduced mode is active)
          if (skip_Jac) then

            ! Flag false
            skip_Jac = .False.

          ! Normal
          else

            ! Get Jabocian
            call Jacobian_Compute(Input,Atom,Atomb,Mol,Geom, &
                                  GeomI,Flgsg,Frec,fudge,kurucz, &
                                  MPID,Atmo,Bfield,Inf_Nodes,Sol, &
                                  SolF,LM_Stru)
            if (laborted) exit

          end if ! If Jacobian is actually computed or not

          ! Flag Jacobian
          Flag_Jac = .True.

        end if ! Broyden or Jacobian

        ! Compute Hessian
        call Hessian_Compute(Inf_Nodes%Nodes_Type, &
                             LM_Stru,Flag_RDmode)

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

        end if ! Regularizing

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
              call Trial_Synthesis(LM_Stru,Inf_Nodes, &
                                   LM_Stru%Lambda,Solution, &
                                   Atmo,Bfield,Atom,Atomb, &
                                   Mol,Geom,GeomI,Flgsg,Frec, &
                                   fudge,kurucz,MPID,Input, &
                                   Sol,SolF)
              if (laborted) exit

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
            call Backtracking(LM_Stru,Atom,Atomb,Mol,Geom, &
                              GeomI,Flgsg,Frec,fudge,kurucz, &
                              MPID,Atmo,Bfield,Sol,SolF, &
                              Inf_Nodes,Input,Inf_Stokes, &
                              Solution,Stokes_Min)

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

              ! If Jacobi not flagged, overflow Broyden index
              if (.not.Flag_Jac) Num_Broyden = 10

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
            if (vlevel.gt.0) call verboseI(3)

            ! Verbose merit function
            write(umsg,'(A,i4,2(3x,A,es15.4))')  &
              ' * Iteration = ',indx_iter, &
              'Chi2 = ',LM_Stru%Chisq, &
              'Chi2 (no regularization) = ',LM_Stru%Chisq_og

            ! Verbose depending on MPI regime
            if (gpid.eq.0) then
              call verboseI(0)
              if (vlevel.gt.0) call verboseI(3)
            else
              call verboseI(0)
              if (vlevel.gt.0) call verboseI(3)
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
          call Convergence_Check(Input,Chisq_old, &
                                 LM_Stru%Chisq,Flag_Convg)

          ! If converged
          if (Flag_Convg) then

            ! If Flagged Jacobian
            if (Flag_Jac) then

              ! Master
              if (pid.eq.0) then

                ! Verbose merit function if in output
                umsg = ' * '
                call verboseI(0)
                if (vlevel.gt.0) call verboseI(3)

                ! Verbose
                write(umsg,'(A,es15.4)') &
                  ' * LM converged. Chi2 = ',LM_Stru%Chisq

                ! Verbose depending on the MPI regime
                if (gpid.eq.0) then
                  call verboseI(0)
                  if (vlevel.gt.0) call verboseI(3)
                else
                  call verboseI(0)
                  if (vlevel.gt.0) call verboseI(3)
                end if

              end if ! Master

              ! If in reduced mode, or not allowed, leave iterations
              if (Flag_RDmode.or..not.Input%allow_RD_mode) exit

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

              ! Verbose depending on the MPI regime
              if (gpid.eq.0) then
                call verboseI(0)
                if (vlevel.gt.0) call verboseI(3)
              else
                call verboseI(0)
                if (vlevel.gt.0) call verboseI(3)
              end if ! MPI regime
            end if ! Master

            ! Exit loop
            exit

          end if ! Last iteration

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
                write(umsg,'(A,es15.4,2(3x,A,es15.4))') &
                  '   New scaled penalty = ', &
                  LM_Stru%Rgl%Penalty*LM_Stru%Rgl%Ratio,  &
                  'Ratio = ',LM_Stru%Rgl%Ratio, &
                  'Penalty = ',LM_Stru%Rgl%Penalty
                call verboseI(3)
                write(umsg,'(A,es15.4,3x,A,es15.4)') &
                  ' * New Chi2 (total) = ', &
                  LM_Stru%Chisq, &
                  'Chi2 (no regulatization) = ', &
                  LM_Stru%Chisq_og
                call verboseI(3)

              end if ! Master
            end if ! Lower ratio, but acceptable
          end if ! Regularizing and ratio larger than 1

          ! If not converged and in reduced mode
          if (.not.Flag_Convg.and.Flag_RDmode) then

            ! Get fraction from last
            frac = (Chisq_RD0 - LM_Stru%chisq)/ &
                   LM_Stru%chisq

            ! If at least 10% difference
            if (frac.ge.0.1d0) then

              ! Deactivate and recalculate weights
              Flag_RDmode = .False.

              ! Refresh Lambda track
              LM_Stru%nLambda = 0
              LM_Stru%Lambda = 1d-1*Input%factoraccept

              ! If thermal inversion
              if (Inf_Nodes%Nodes_Type.eq.0) then

                ! If modified weights
                if (allocated(LM_Stru%WeightI_mod)) then
                  MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%WeightI_mod)
                  deallocate(LM_Stru%WeightI_mod)
                end if

              ! If magnetic
              else

                ! If modified weights
                if (allocated(LM_Stru%Weight_mod)) then
                  MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Weight_mod)
                  deallocate(LM_Stru%Weight_mod)
                end if

              end if ! Type of inversion

              ! Master
              if (pid.eq.0) then

                ! Verbose
                write(umsg,'(A)') ' +'
                call verboseI(3)
                write(umsg,'(A,3(1x,f6.3))') &
                           ' + Deactivated Reduced Mode'
                call verboseI(3)

              end if ! Master
            end if ! At least 10% difference
          end if ! Not converged and reduced mode

          ! Update chi2
          Chisq_old = LM_Stru%Chisq

        end if ! Accepted step

        ! If converged because RC
        if ((Flag_Convg.and. &
             LM_Stru%chisq.ge.Input%Threshold_chisq).or. &
            .not.LM_Stru%accepted) then
             
          ! If not in RD mode
          if (.not.Flag_RDmode) then

            ! Activate
            Flag_RDmode = .True.

            ! Redefine the weights
            call adjust_weights(LM_Stru,Inf_Nodes, &
                                Inf_Stokes,Input)

            ! If not accepted, keep Jacobian
            if (.not.LM_Stru%accepted) skip_Jac = .True.

            ! Refresh Lambda track
            LM_Stru%nLambda = 0
            LM_Stru%Lambda = 1d-1*Input%factoraccept

            ! Save current chi^2
            Chisq_RD0 = LM_Stru%chisq

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A)') ' +'
              call verboseI(3)
              write(umsg,'(A,3(1x,f6.3))') &
                                ' + Activate Reduced Mode: ', &
                                LM_Stru%perc,LM_Stru%ff_max, &
                                LM_Stru%ff_contr
              call verboseI(3)

            end if ! Master

            ! And cycle
            cycle

          ! Already in RD mode
          else

            ! If not accepted and master
            if (.not.LM_Stru%accepted.and.pid.eq.0) then

              ! Verbose
              umsg = ' * '
              call verboseI(0)
              if (vlevel.gt.0) call verboseI(3)

              ! Verbose
              write(umsg,'(A)') ' * Could not improve more'

              ! Verbose depending on the MPI regime
              if (gpid.eq.0) then
                call verboseI(0)
                if (vlevel.gt.0) call verboseI(3)
              else
                call verboseI(0)
                if (vlevel.gt.0) call verboseI(3)
              end if

            end if ! Master and not accepted

            ! Just leave
            exit

          end if ! If in RD mode or not
        end if ! If failed

        ! Predict lambda if in accepted step
        if (Input%l_Lam_track.and.LM_Stru%accepted) &
          call update_lambda(Lam_track,Input%Lam_track, &
                             LM_Stru%nLambda,LM_Stru%Lambda)

      end do ! Iterations

      ! If error
      if (laborted) goto 2000

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
          call Hessian_Compute(Inf_Nodes%Nodes_Type,LM_Stru,.False.)

        end if

      ! Not accepted last step or error not from Hessian
      else

        ! If in reduced mode, recalculate Hessian
        if (Flag_RDmode) &
          call Hessian_Compute(Inf_Nodes%Nodes_Type,LM_Stru,.False.)

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
      call write_result_inv(Atom,Atomb,Mol,GeomI,Geom,Frec, &
                            fudge,Atmo,Bfield,Input, &
                            Inf_Stokes,Inf_Nodes,Sol,SolF, &
                            LM_Stru,saving,.True.)


      ! If thermal inversion
2000  if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Deallocate intensity quantities
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%ResidualI)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%WeightI)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%JacobianI)
        deallocate(LM_Stru%ResidualI)
        deallocate(LM_Stru%WeightI)
        deallocate(LM_Stru%JacobianI)

        ! If modified weights
        if (allocated(LM_Stru%WeightI_mod)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%WeightI_mod)
          deallocate(LM_Stru%WeightI_mod)
        end if

      ! If magnetic
      else

        ! Deallocate polarization quantities
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Residual)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Weight)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacobian)
        deallocate(LM_Stru%Residual)
        deallocate(LM_Stru%Weight)
        deallocate(LM_Stru%Jacobian)

        ! If modified weights
        if (allocated(LM_Stru%Weight_mod)) then
          MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Weight_mod)
          deallocate(LM_Stru%Weight_mod)
        end if

      end if ! Type of inversion

      ! Free memory
      MRAMc = MRAMc - 1d-6*sizeof(Solution)
      MRAMc = MRAMc - 1d-6*sizeof(Errors)
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacfvec)
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%diag)
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Hessian)
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Hessian_og)
      MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Jacfvec_og)
      deallocate(Solution,Errors,LM_Stru%Jacfvec,LM_Stru%diag)
      deallocate(LM_Stru%Hessian,LM_Stru%Hessian_og)
      deallocate(LM_Stru%Jacfvec_og)
      if (Input%l_Lam_track) then
        MRAMc = MRAMc - 1d-6*sizeof(Lam_track)
        deallocate(Lam_track)
      end if

      ! If regularizing, free memory
      if (Inf_Nodes%Regul_Flag.and. &
          allocated(LM_Stru%Rgl%Regul_H)) then
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Rgl%Regul_H)
        MRAMc = MRAMc - 1d-6*sizeof(LM_Stru%Rgl%Regul_F)
        deallocate(LM_Stru%Rgl%Regul_H)
        deallocate(LM_Stru%Rgl%Regul_F)
      end if

      ! Free solutions
      call free_inv_solution(SolF)

      ! Free innate memory of local structure
      SRAMc = SRAMc - 1d-6*sizeof(SolF)

      return

      end subroutine LMFIT

!#####################################################################
!#####################################################################
!#####################################################################

      !> Prepare the model atmosphere and Stokes profiles to write the
      !! current result in the output and call such writing\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!       Atmo_in(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!    SolF(Solution_F_class): Structure with the solution of
      !!                            the self-consistent problem and
      !!                            the corresponding emergent
      !!                            profiles, contribution function,
      !!                            and height for optical depth
      !!                            equal to one\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!           saving(logical): If the result is to be stored\n
      !!           ifinal(logical): If the inversion is finished
      subroutine write_result_inv(Atom,Atomb,Mol,GeomI,Geom,Frec, &
                                  fudge,Atmo_in,Bfield,Input, &
                                  Inf_Stokes,Inf_Nodes,Sol,SolF, &
                                  LM_Stru,saving,ifinal)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Geometry_class), intent(in):: GeomI, Geom
      type(Frequency_class), intent(in):: Frec
      type(fudge_class), intent(in):: fudge
      type(Atmo_class), intent(inout):: Atmo_in
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(LMFIT_class), intent(inout):: LM_Stru
      type(Solution_F_class), intent(in):: SolF
      logical, intent(in):: saving
      logical, intent(in):: ifinal

      ! Local
      type(Atmo_class):: Atmo

      ! Control
      call control

      ! If final
      if (ifinal) then

        ! Prepare output model atmosphere
        call setup_Atmo_ouinv(Atom,Atomb,Mol,Atmo_in,Input,fudge)

        ! Get a copy of model atmosphere
        call cAtmo(Atmo_in,Atmo)

      ! If not final
      else

        ! Get a copy of model atmosphere
        call cAtmo(Atmo_in,Atmo)

        ! Prepare output model atmosphere
        call setup_Atmo_ouinv(Atom,Atomb,Mol,Atmo,Input,fudge)

      end if

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
                Input%Type_inversion.ne.4.and. &
                Input%Type_inversion.ne.5) then

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
              call writestk(Input%folder,0,0,Frec%omega,Geom, &
                            SolF%e_Stk_b(:,:,1,1),Input%lim_stk)
              if (laborted) exit


              ! Output contribution
              if (Input%out_contr) then
                call writectr_inv(Input%folder, &
                                  SolF%e_Ctr_b(:,:,:,1,1), &
                                  Input%lim_ctr)
                if (laborted) exit
              end if

            end if ! Thermal or not

            ! Thermal not from a full inversion or non-thermal
            if ((Inf_Nodes%Nodes_type.eq.0.and. &
                 Input%Type_inversion.ne.3.and. &
                 Input%Type_inversion.ne.4.and. &
                 Input%Type_inversion.ne.5).or. &
                (Inf_Nodes%Nodes_type.ne.0)) then

              ! Output tau1
              if (Input%out_tau1) &
                call writetau_inv(Input%folder, &
                                  SolF%e_tau1_b(:,1,1), &
                                  Input%lim_tau)

            end if ! Thermal not from full or non-thermal

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

      ! Wipe the copy clean
      call free_Atmo(Atmo,.True.)

      return

      end subroutine write_result_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check the convergence of the merit function\n
      !!   Input(Input_class): Structure with configuration data\n
      !!    Chisq_old(double): Previous chi^2\n
      !!        Chisq(double): Current chi^2\n
      !!  Flag_Convg(logical): If chi^2 is converged
      subroutine Convergence_Check(Input,Chisq_old,Chisq,Flag_Convg)

      ! I/O

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

      !> Calculate an estimation of the error in the inverted
      !! parameters\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data
      subroutine Err(LM_Stru,Input,Inf_Stokes,Sol,Inf_Nodes)

      ! I/O

      type(LMFIT_class), intent(in):: LM_Stru
      type(Input_class), intent(in):: Input
      type(Stokes_class), intent(in):: Inf_Stokes
      type(Solution_class), intent(in):: Sol
      type(Nodes_class), intent(inout):: Inf_Nodes

      ! Local

      integer:: tmp,i,j,k

      double precision:: Error1,Error2,Num,Den
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

            end do ! Inverted variables

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

      !> Give a Levenberg-Marquardt step and compute emergent
      !! profiles\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!            Lambda(double): Levenberg-Marquardt lambda
      !!                            parameter\n
      !!       Solution(double(:)): SVD solution\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!    SolF(Solution_F_class): Structure with the solution of
      !!                            the self-consistent problem and
      !!                            the corresponding emergent
      !!                            profiles, contribution function,
      !!                            and height for optical depth
      !!                            equal to one
      subroutine Trial_Synthesis(LM_Stru,Inf_Nodes,Lambda,Solution, &
                                 Atmo,Bfield,Atom,Atomb,Mol,Geom, &
                                 GeomI,Flgsg,Frec,fudge,kurucz,MPID, &
                                 Input,Sol,SolF)

      ! I/O

      type(LMFIT_class), intent(inout):: LM_Stru
      type(Nodes_class), intent(in):: Inf_Nodes
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
      type(Input_class), intent(inout):: Input
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      double precision, intent(in):: Lambda
      double precision, dimension(:), &
                        allocatable, intent(inout):: Solution

      ! Local

      type(Atmo_class):: Tmp_Atmo
      type(Bfield_class):: Tmp_Bfield
      type(Nodes_class):: Inf_Nodes_tmp

      integer:: i

      double precision, dimension(:), allocatable:: Jacfvec_new
      double precision, dimension(:,:), allocatable:: Hessian_new

      ! Count memory
      ! Neglecting the data in Inf_Nodes_tmp%Node and
      ! Inf_Nodes%Inf_Inv
      MRAMc = MRAMc + 1d-6*sizeof(Tmp_Atmo)
      MRAMc = MRAMc + 1d-6*sizeof(Inf_Nodes_tmp)

      ! Copy current Nodes
      Inf_Nodes_tmp = Inf_Nodes

      ! Allocate new Hessian and Jacobian vector
      allocate(Hessian_new(LM_Stru%Num,LM_Stru%Num))
      allocate(Jacfvec_new(LM_Stru%Num))
      MRAMc = MRAMc + 1d-6*sizeof(Hessian_new)
      MRAMc = MRAMc + 1d-6*sizeof(Hessian_new)

      ! Copy current Hessian and Jacobian vectors
      Hessian_new = LM_Stru%Hessian
      Jacfvec_new = LM_Stru%Jacfvec

      ! Add lambda contribution to diagonal in Hessian matrix
      do i = 1, LM_Stru%Num
        Hessian_new(i,i) = Hessian_new(i,i) + Lambda*LM_Stru%diag(i)
      end do

      ! SVD solution
      call SVD_Solve(Hessian_new,Jacfvec_new,Solution, &
                     LM_Stru%Num,Inf_Nodes_tmp,Input%SVD_Type, &
                     .True.)
      if (laborted) goto 2000

      ! Predict improvement and check step
      call predict_improvement(LM_Stru%Hessian,LM_Stru%Jacfvec, &
                               Solution,LM_Stru%pred)

      ! Save norm (max)
      LM_Stru%step_norm = maxval(abs(Solution))

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
        call Intpol_Atmo(Inf_Nodes_tmp,Tmp_Atmo,Atom,Atomb, &
                         Mol,Input,fudge)

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Tmp_Atmo,Bfield,Input, &
                        Sol,SolF,1)

        ! Wipe Tmp_Atmo
        call free_Atmo(Tmp_Atmo,.True.)

      ! If magnetic inversion only
      else if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Get a copy of the magnetic field
        call cBfield(Bfield,Tmp_Bfield)

        ! Generate new stratification
        call Intpol_Bfield(Inf_Nodes_tmp, Atmo, Tmp_Bfield)
        if (laborted) goto 2000

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Tmp_Bfield,Input, &
                        Sol,SolF,1)

      ! If inverting all
      else if (Inf_Nodes%Nodes_Type.eq.2) then

        ! Get a copy of atmosphere and field
        call cAtmo(Atmo,Tmp_Atmo)
        call cBfield(Bfield,Tmp_Bfield)

        ! Get new stratification
        call Intpol_Atmo_all(Inf_Nodes_tmp,Tmp_Atmo,Tmp_Bfield, &
                             Atom,Atomb,Mol,Input,fudge)

        ! Get Stokes profiles
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Tmp_Atmo,Tmp_Bfield, &
                        Input,Sol,SolF,1)

        ! Wipe Tmp_Atmo
        call free_Atmo(Tmp_Atmo,.True.)

      ! Error
      else

        ! Aborting
        umsg = 'The index of the node type is not correct'
        urou = 'Trial_Synthesis'
        call aborted

      end if ! Type of inversion

      ! Deallocate auxiliar Hessian and Jacobian
2000  MRAMc = MRAMc - 1d-6*sizeof(Hessian_new)
      MRAMc = MRAMc - 1d-6*sizeof(Hessian_new)
      deallocate(Hessian_new,Jacfvec_new)

      ! Remove non-array memory of local structures
      MRAMc = MRAMc - 1d-6*sizeof(Tmp_Atmo)
      MRAMc = MRAMc - 1d-6*sizeof(Inf_Nodes_tmp)

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

      !> Optimize the Levenberg-Marquardt lambda parameter with the
      !! backtracking algorithm\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      Atomb(Atom_class(:)): Structures with atomic data for
      !!                            background atoms\n
      !!         Mol(Mol_class(:)): Structures with molecular data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!     GeomI(Geometry_class): Structure with geometric data for
      !!                            the intensity problem\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        fudge(fudge_class): Structure with fudge data\n
      !!      kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!       Sol(Solution_class): Structure with the frequency and
      !!                            synthetic Stokes parameters in the
      !!                            frequency range of the inverted
      !!                            data\n
      !!    SolF(Solution_F_class): Structure with the solution of
      !!                            the self-consistent problem and
      !!                            the corresponding emergent
      !!                            profiles, contribution function,
      !!                            and height for optical depth
      !!                            equal to one\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!        Input(Input_class): Structure with configuration
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!   Solution_Min(double(:)): SVD solution\n
      !!   Stokes_Min(double(:,:)): Currently best Stokes parameters
      subroutine Backtracking(LM_Stru,Atom,Atomb,Mol,Geom,GeomI, &
                              Flgsg,Frec,fudge,kurucz,MPID,Atmo, &
                              Bfield,Sol,SolF,Inf_Nodes,Input, &
                              Inf_Stokes,Solution_Min,Stokes_Min)

      ! I/O

      type(LMFIT_class), intent(inout):: LM_Stru
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(Geometry_class), intent(inout):: GeomI, Geom
      type(Frequency_class), intent(inout):: Frec
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(MPI_class), intent(inout):: MPID
      type(Atmo_class), intent(inout):: Atmo
      type(Bfield_class), intent(inout):: Bfield
      type(Input_class), intent(inout):: Input
      type(Stokes_class), intent(inout):: Inf_Stokes
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      double precision, dimension(:), &
                        allocatable, intent(inout):: Solution_Min
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: Stokes_Min

      ! Local

      logical:: Bracketed,up_first,converged,changed
      logical:: second_lap,early_exit,could_be_gd

      integer, parameter:: Length = 10
      integer, parameter:: Max_fail = 3
      integer:: indx,jndx,kndx,Chisq_indx,nfail

      double precision:: Chisq_old,Chisq_best,daux,rel,last_norm
      double precision:: last_best,relb
      double precision:: Lambda_TMP,Penalty_TMP
      double precision, dimension(:), allocatable:: Lambda_array
      double precision, dimension(:), allocatable:: Chisq_array
      double precision, dimension(:), allocatable:: Penalty_array
      double precision, dimension(:), allocatable:: Solution


      ! Create auxiliar arrays
      allocate(Lambda_array(Length+1),Chisq_array(Length+1))
      allocate(Solution(LM_Stru%Num),Penalty_array(Length+1))
      MRAMc = MRAMc + 1d-6*sizeof(Lambda_array)
      MRAMc = MRAMc + 1d-6*sizeof(Chisq_array)
      MRAMc = MRAMc + 1d-6*sizeof(Solution)
      MRAMc = MRAMc + 1d-6*sizeof(Penalty_array)

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
      Chisq_old = LM_Stru%Chisq
      up_first = .True.
      Chisq_best = chisq_old
      jndx = 1
      second_lap = .False.
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

      ! Fake loop
      do while (.True.)

        ! No GD candidate
        could_be_gd = .False.

        ! Reset last registered norm
        last_norm = 1e9

        ! For the hard-coded length
        do indx=jndx,Length

          ! Initialize failures
          nfail = 0

          ! Fake loop
          do while (.True.)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,i2)') ' +'
              call verboseI(3)
              write(umsg,'(A,i2)') ' + Trial ',indx
              call verboseI(3)
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

            ! Try new solution
            call Trial_Synthesis(LM_Stru,Inf_Nodes, &
                                 Lambda_array(indx),Solution,Atmo, &
                                 Bfield,Atom,Atomb,Mol,Geom, &
                                 GeomI,Flgsg,Frec,fudge,kurucz, &
                                 MPID,Input,Sol,SolF)

            ! If the Trial did not fail, calculate here the
            ! merit function
            if (.not.laborted) then

              ! Get merit function
              call Merit_function(Inf_Stokes, Sol%Stokes_out, &
                                  Inf_Nodes%Nodes_Type, LM_Stru)

              ! If the merit function is not a number, consider
              ! it a failure
              if (.not.ieee_is_finite(LM_Stru%Chisq)) then

                ! Flag
                laborted = .True.

                ! Verbose
                if (gpid.eq.0) then
                  umsg = ' # Trial synthesis lead to a '// &
                         'non-sensical merit function'
                  call verbose
                end if
                umsg = ' # Trial synthesis lead to a '// &
                       'non-sensical merit function'
                call verboseI(3)

              end if ! The merit function is wrong

            ! If the trial failed
            else

              ! Verbose
              if (gpid.eq.0) then
                umsg = ' # Trial synthesis failed'
                call verbose
              end if
              umsg = ' # Trial synthesis failed'
              call verboseI(3)

            end if ! The trial did not fail

            ! If there was an issue with the trial
            if (laborted) then

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
                if (changed) cycle

                ! Leave
                exit

              end do ! Check repetition

              ! If lambda below limits
              if (Lambda_array(indx).lt. &
                  LM_Stru%Lambda_bounds(1)) then

                ! If there is previous
                if (indx.gt.1) then

                  ! If already was below limit
                  if (Lambda_array(indx-1).le. &
                      LM_Stru%Lambda_bounds(1)) then

                    ! Verbose
                    if (gpid.eq.0) then
                      umsg = ' # Could not find new lambda to '// \
                             'try with'
                      call verbose
                    end if
                    umsg = ' # Could not find new lambda to try with'
                    call verboseI(3)

                    ! Exit trials
                    early_exit = .True.
                    exit

                  end if ! Already was below limits
                end if ! Not the first

                ! Force limit
                Lambda_array(indx) = LM_Stru%Lambda_bounds(1)

              ! If Lambda above limits
              else if (Lambda_array(indx).gt. &
                  LM_Stru%Lambda_bounds(2)) then

                ! If there is previous
                if (indx.gt.1) then

                  ! If already was above limit
                  if (Lambda_array(indx-1).ge. &
                      LM_Stru%Lambda_bounds(2)) then

                    ! Verbose
                    if (gpid.eq.0) then
                      umsg = ' # Could not find new lambda to '// &
                             'try with'
                      call verbose
                    end if
                    umsg = ' # Could not find new lambda to try with'
                    call verboseI(3)

                    ! Exit trials
                    early_exit = .True.
                    exit

                  end if ! Already was below limits
                end if ! Not the first

                ! Force limit
                Lambda_array(indx) = LM_Stru%Lambda_bounds(2)

              end if ! If Lambda beyond limits

              ! Add to count
              nfail = nfail + 1

              ! If beyond saving
              if (nfail.ge.Max_fail) goto 1000

              ! Verbose
              if (gpid.eq.0) then
                umsg = ' # Try new lambda value'
                call verbose
              end if
              umsg = ' # Try new lambda value'
              call verboseI(3)

              ! Try again
              cycle

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
            write(umsg,'(4(3x,A,es15.4))') &
              'chi2 = ',LM_Stru%Chisq, &
              'lambda = ',Lambda_array(indx), &
              'Penalty = ',LM_Stru%Rgl%Penalty, &
              'Ratio = ',LM_Stru%Rgl%Ratio
            call verboseI(3)
            write(umsg,'(A,es15.4)') '   Improvement ratio ', &
                                     (Chisq_old - LM_Stru%Chisq)/ &
                                     LM_Stru%pred
            call verboseI(3)

          end if ! Master

          ! If larger than best
          if (Chisq_array(indx).gt.Chisq_best) then

            ! Get relative change
            rel = (Chisq_array(indx) - Chisq_old) / &
                   Chisq_old

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

                ! Check how bad it was
                daux = (Chisq_array(indx-1) - Chisq_best) / &
                        Chisq_best

              ! Not any of the first two
              else

                ! It is bracketed, we can leave
                Bracketed = .True.
                exit

              end if ! Index

              ! If too wrong
              if (rel.gt.1e4) then

                  ! Get next lambda
                  Lambda_array(indx+1) = Lambda_array(indx)* &
                                         LM_Stru%factorreject* &
                                         LM_Stru%factorreject* &
                                         LM_Stru%factorreject

              ! If quite wrong
              else if (rel.gt.1e2) then

                  ! Get next lambda
                  Lambda_array(indx+1) = Lambda_array(indx)* &
                                         LM_Stru%factorreject* &
                                         LM_Stru%factorreject

              ! Normal
              else

                  ! Get next lambda
                  Lambda_array(indx+1) = Lambda_array(indx)* &
                                         LM_Stru%factorreject

              end if ! How wrong was it

            ! Is the worse result
            else

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

              ! Flag
              could_be_gd = .False.

            ! First already found
            else

              ! Update best index and chi2
              Chisq_indx = indx
              Chisq_best = Chisq_array(indx)

              ! Not the first
              if (indx.gt.1) then

                ! Larger lambda
                if (Lambda_array(indx).gt.Lambda_array(indx-1)) then

                  ! We had one potential GD
                  if (could_be_gd) then

                    ! GD regime
                    if (LM_Stru%step_norm/last_norm.gt..9d0) then

                      ! Check convergence
                      rel = (last_best - Chisq_best)/last_best

                      ! If too small an improvement
                      if (rel.lt.1d-2) then

                        ! Master
                        if (pid.eq.0) then

                          ! Verbose
                          write(umsg,'(A,es15.4)') &
                               ' # Seems like it is exploring '// &
                               'the GD regime without gain'
                          call verboseI(3)
                          write(umsg,'(A,es15.4,A,es15.4)') &
                              '   chi2 changed only from ',last_best, &
                              ' to ',Chisq_best
                          call verboseI(3)
                          write(umsg,'(A,es15.4,A,es15.4)') &
                              '   step norm changed only '// &
                              'from ',last_norm, &
                              ' to ',LM_Stru%step_norm
                          call verboseI(3)
                          write(umsg,'(A,es15.4,A,es15.4)') &
                              '   stop Backtracking'
                          call verboseI(3)

                        end if ! Master

                        ! Leave
                        exit

                      end if ! GD regime

                    ! There was improvement
                    else

                      ! Keep this new chi2
                      last_best = Chisq_best

                    end if

                  ! No potential GD yet
                  else

                    ! Flag
                    could_be_gd = .True.
                    last_best = Chisq_best

                  end if ! GD flagged
                end if ! Larger lambda
              end if ! Not the first

              ! Increment
              Lambda_array(indx+1) = Lambda_array(indx)* &
                                     LM_Stru%factorreject

            end if ! Found the first worse result

          end if ! Improvement or not

          ! If Lambda below limit
          if (Lambda_array(indx+1).lt.LM_Stru%Lambda_bounds(1)) then

            ! If previous was already at limit, exit
            if (Lambda_array(indx).le.LM_Stru%Lambda_bounds(1)) exit

            ! If not put to limit
            Lambda_array(indx+1) = LM_Stru%Lambda_bounds(1)

          ! If lambda above limit
          else if (Lambda_array(indx+1).gt. &
                   LM_Stru%Lambda_bounds(2)) then

            ! If previous was already at limit, exit
            if (Lambda_array(indx).ge.LM_Stru%Lambda_bounds(2)) exit

            ! If not put to limit
            Lambda_array(indx+1) = LM_Stru%Lambda_bounds(2)

          end if ! Lambda in bounds

          ! Save last step norm
          last_norm = LM_Stru%step_norm

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
          call Trial_Synthesis(LM_Stru,Inf_Nodes,Lambda_TMP, &
                               Solution,Atmo,Bfield,Atom,Atomb, &
                               Mol,Geom,GeomI,Flgsg,Frec,fudge, &
                               kurucz,MPID,Input, Sol,SolF)

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
              if (minval(Lambda_array(1:2)).gt. &
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
                  write(umsg,'(3(3x,A,es15.4))') &
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
              if (maxval(Lambda_array(1:2)).lt. &
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

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,i2)') ' +'
        call verboseI(3)
        write(umsg,'(A,2(3x,A,es15.4))') &
          ' + Leave backtracking with:', &
          'chi2 = ',LM_Stru%Chisq, &
          'lambda = ',LM_Stru%Lambda
        call verboseI(3)

      end if ! Master

      ! Free memory
1000  MRAMc = MRAMc - 1d-6*sizeof(Lambda_array)
      MRAMc = MRAMc - 1d-6*sizeof(Chisq_array)
      MRAMc = MRAMc - 1d-6*sizeof(Solution)
      MRAMc = MRAMc - 1d-6*sizeof(Penalty_array)
      deallocate(Lambda_array,Chisq_array,Solution,Penalty_array)

      return

      end subroutine Backtracking

!#####################################################################
!#####################################################################
!#####################################################################

      !> Propose a new Levenberg-Marquardt lambda parameter\n
      !!  LM_Stru(LMFIT_class): Structure with data for the
      !!                        Levenberg–Marquardt\n
      !!        Lambda(double): Levenberg–Marquardt lambda parameter
      subroutine Lambda_propose_fix(LM_Stru,Lambda)

      ! I/O
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

      !> Propose the value of the Levenberg-Marquardt lambda parameter
      !! for the next iteration\n
      !!      track(double(:)): Lambda parameter history\n
      !!       ntrack(integer): Size of history stored\n
      !!  LM_Stru(LMFIT_class): Structure with data for the
      !!                        Levenberg–Marquardt
      subroutine predict_lambda(track,ntrack,LM_Stru)

      ! I/O

      integer, intent(in):: ntrack
      double precision, dimension(:), intent(in):: track
      type(LMFIT_class), intent(inout):: LM_Stru

      ! Local

      double precision:: d1,d2,d3
      double precision, dimension(ntrack):: x,a,b,c


      ! If first iteration, skip
      if (LM_Stru%nLambda.lt.1) return

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

            ! Verbose depending on MPI regime
            if (gpid.eq.0) then
              call verboseI(3)
            else
              call verboseI(3)
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
          if (LM_Stru%nLambda.lt.2) then

            ! Get new lambda
            LM_Stru%Lambda = track(2)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,1x,es15.4)')  &
                ' - Tracking lambda:',track(2)

              if (gpid.eq.0) then
                call verboseI(3)
              else
                call verboseI(3)
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

              ! Verbose depending on MPI regime
              if (gpid.eq.0) then
                call verboseI(3)
              else
                call verboseI(3)
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
          if (LM_Stru%nLambda.lt.2) then

            ! Get new lambda
            LM_Stru%Lambda = track(3)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,1x,es15.4)')  &
                ' - Tracking lambda:',track(3)

              ! Verbose depending on MPI regime
              if (gpid.eq.0) then
                call verboseI(3)
              else
                call verboseI(3)
              end if
            end if

            ! Sanity check
            if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1).or.&
                LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) &
                LM_Stru%Lambda = track(3)

          ! Only two iterations
          else if (LM_Stru%nLambda.lt.3) then

            ! Get new lambda
            LM_Stru%Lambda = 2*track(3) - track(2)

            ! Master
            if (pid.eq.0) then

              ! Verbose
              write(umsg,'(A,2(1x,es15.4),A,es15.4)')  &
                ' - Tracking lambda:',track(2),track(3), &
                ' -> ',LM_Stru%Lambda

              ! Verbose depending on MPI regime
              if (gpid.eq.0) then
                call verboseI(3)
              else
                call verboseI(3)
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

              ! Verbose depending on MPI regime
              if (gpid.eq.0) then
                call verboseI(3)
              else
                call verboseI(3)
              end if
            end if

            ! Sanity check
            if (LM_Stru%Lambda.lt.LM_Stru%Lambda_bounds(1).or.&
                LM_Stru%Lambda.gt.LM_Stru%Lambda_bounds(2)) then

                ! Set
                LM_Stru%Lambda = track(3)

            ! Invalid value
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

        ! Verbose depending on MPI regime
        if (gpid.eq.0) then
          call verboseI(3)
        else
          call verboseI(3)
        end if
      end if

      return

      end subroutine predict_lambda

!#####################################################################
!#####################################################################
!#####################################################################

      !> Update history of Levenberg-Marquardt lambda parameter\n
      !!  track(double(:)): Lambda parameter history\n
      !!   ntrack(integer): Size of history stored\n
      !!  nlambda(integer): Number of updates\n
      !!    lambda(double): Levenberg-Marquardt lambda parameter
      subroutine update_lambda(track,ntrack,nlambda,lambda)

      ! I/O
      integer, intent(in):: ntrack
      integer, intent(inout):: nlambda
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
      nlambda = nlambda + 1

      return

      end subroutine update_lambda

!#####################################################################
!#####################################################################
!#####################################################################

      !> Modify the parameter node values according to the SVD
      !! solution\n
      !!     Solution(double(:)): SVD solution\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      subroutine Nodes_Modify(Solution,Inf_Nodes)

      ! I/O

      type(Nodes_class), intent(inout):: Inf_Nodes
      double precision, dimension(:), &
                        allocatable, intent(in):: Solution

      ! Local

      integer:: tmp,i


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

      !> Enforce limits of Levenberg-Marquardt lambda parameter\n
      !!            Lambda(double): Levenberg-Marquardt lambda
      !!                            parameter\n
      !!  Lambda_limits(double(:)): Boundary limits for
      !!                            Levenberg-Marquardt lambda
      !!                            parameter
      subroutine CheckLambda(Lambda,Lambda_limits)

      ! I/O

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

      !> Set a self-consistent solution as the best\n
      !!  Sol(Solution_class): Structure with the frequency and
      !!                       synthetic Stokes parameters in the
      !!                       frequency range of the inverted data\n
      !!        best(logical): If the proposed solution is really the
      !!                       best and not just the best of the
      !!                       backtracking\n
      !!        copy(logical): If the solution comes from the
      !!                       backtracking
      subroutine set_best(Sol,best,copy)

      ! I/O

      type(Solution_F_class), intent(inout):: Sol
      logical, intent(in):: best, copy

      ! Local

      logical:: ex1,ex2,ex3


      ! Slaves, leave
      if (pid.gt.0) return

      ! Truly the best
      if (best) then

        ! Existence flags
        ex1 = allocated(Sol%e_Stk_b)
        ex2 = allocated(Sol%e_tau1_b)
        ex3 = allocated(Sol%e_Ctr_b)

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
          if (allocated(Sol%i_JKQS_t)) Sol%i_JKQS_b = Sol%i_JKQS_t
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
          if (allocated(Sol%i_JKQS)) Sol%i_JKQS_b = Sol%i_JKQS
          if (allocated(Sol%i_JKQC)) Sol%i_JKQC_b = Sol%i_JKQC
          if (allocated(Sol%i_rhoes_t)) Sol%i_rhoes_b = Sol%i_rhoes

        end if ! Copying from backtrace?

        ! Memory
        if (.not.ex1) SRAMc = SRAMc + 1d-6*sizeof(Sol%e_Stk_b)
        if (.not.ex2.and.allocated(Sol%e_Ctr_b)) &
          SRAMc = SRAMc + 1d-6*sizeof(Sol%e_Ctr_b)
        if (.not.ex3.and.allocated(Sol%e_tau1_b)) &
          SRAMc = SRAMc + 1d-6*sizeof(Sol%e_tau1_b)

      ! Provisional best
      else

        ! Existence flags
        ex1 = allocated(Sol%e_Stk_t)
        ex2 = allocated(Sol%e_tau1_t)
        ex3 = allocated(Sol%e_Ctr_t)

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
        if (allocated(Sol%i_JKQS)) Sol%i_JKQS_t = Sol%i_JKQS
        if (allocated(Sol%i_JKQC)) Sol%i_JKQC_t = Sol%i_JKQC
        if (allocated(Sol%i_rhoes)) Sol%i_rhoes_t = Sol%i_rhoes

        ! Memory
        if (.not.ex1) SRAMc = SRAMc + 1d-6*sizeof(Sol%e_Stk_t)
        if (.not.ex2.and.allocated(Sol%e_Ctr_t)) &
          SRAMc = SRAMc + 1d-6*sizeof(Sol%e_Ctr_t)
        if (.not.ex3.and.allocated(Sol%e_tau1_t)) &
          SRAMc = SRAMc + 1d-6*sizeof(Sol%e_tau1_t)

      end if ! Truly the best

      end subroutine set_best

!#####################################################################
!#####################################################################
!#####################################################################

      end module lmfit_mod
