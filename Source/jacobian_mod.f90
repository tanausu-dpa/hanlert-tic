      !> Manages the jacobian for the LM
      module jacobian_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC/NSSCC)
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     24/02/2023
!  Last version:
!     18/12/2025 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     18/12/2025:    V4.0.2 - Added routines predict_improvement,
!                             angle_and_norm,
!                             get_percentile_and_maximum,
!                             adjust_weights_hessian,
!                             adjust_weights_trial, and
!                             adjust_weights (TdPA)
!                           - Changed Hessian_compute routine to
!                             allow for the use of two different
!                             weight variables (TdPA)
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
!  Merit_function
!    Compute the L2 merit function of Stokes profiles
!
!  predict_improvement
!    Predict the improvement for the calculated step
!
!  angle_and_norm
!    Get angle between the two vectors and their norm
!
!  get_percentile_and_maximum
!    Get percentile and maximum of a given vector ignoring zeros
!
!  adjust_weights_hessian
!    Compute the Hessian from the Jacobian with given weights
!
!  adjust_weights_trial
!    Modify weights to enhance larger differences in the fit for a
!  given set of parameters
!
!  adjust_weights
!    Modify weights to enhance larger differences in the fit
!
!  Jacobian_Compute
!    Compute the Jacobian for every inversion parameter
!
!  Hessian_Compute
!    Compute the Hessian from the Jacobian
!
!  Broyden_Rank1
!    Rank 1 update of the Jacobian following Broyden's method
!
!#####################################################################
!#####################################################################
!#####################################################################

      use commons_mod
      use qsort_mod
      use rf_mod
      use svd_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the L2 merit function of Stokes profiles\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!   Stokes_out(double(:,:)): Stokes parameters\n
      !!       Nodes_type(integer): Type of inversion\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt
      subroutine Merit_function(Inf_Stokes,Stokes_out,Nodes_Type, &
                                LM_Stru)

      ! I/O

      type(Stokes_class), intent(in):: Inf_Stokes
      type(LMFIT_class), intent(inout):: LM_Stru
      integer, intent(in):: Nodes_Type
      double precision, dimension(:,:), &
                        allocatable, intent(in):: Stokes_out


      ! Initialize chi
      LM_Stru%Chisq_og = 0d0

      ! If thermal
      if (Nodes_Type.eq.0) then

        ! Compute difference
        LM_Stru%ResidualI = Stokes_out(0,:) - &
                            Inf_Stokes%Stokes_Ob(0,:)

        ! If weight not flagged in LM structure
        if (.not.LM_Stru%Flag_weight) then

          ! If wavelength dependent sigma
          if (Inf_Stokes%Sigma_Flag) then

            ! Compute weight
            LM_Stru%WeightI = Inf_Stokes%weight(0,:)* &
                              Inf_Stokes%weight(0,:)/ &
                              (Inf_Stokes%Sigma_W(0,:)* &
                               Inf_Stokes%Sigma_W(0,:))

          ! No wavelength dependent sigma
          else

            ! Weight
            LM_Stru%WeightI = Inf_Stokes%weight(0,:)* &
                              Inf_Stokes%weight(0,:)

          end if ! Frequency dependent sigma

          ! Flag weight as computed
          LM_Stru%Flag_weight = .True.

        end if ! Flagged weight in LM

        ! Get chi^2
        LM_Stru%Chisq_og = sum(LM_Stru%ResidualI* &
                               LM_Stru%ResidualI* &
                               LM_Stru%WeightI)
      ! Non-thermal
      else

        ! Compute difference
        LM_Stru%Residual = Stokes_out - &
                           Inf_Stokes%Stokes_Ob

        ! If weight not flagged in LM structure
        if (.not.LM_Stru%Flag_weight) then

          ! Weight contribution
          LM_Stru%weight = Inf_Stokes%Weight*Inf_Stokes%Weight

          ! If wavelength dependent sigma
          if (Inf_Stokes%Sigma_Flag) then

            ! Compute weight
            LM_Stru%Weight = LM_Stru%weight/ &
                             (Inf_Stokes%Sigma_W*Inf_Stokes%Sigma_W)

          end if

          ! Flag weight as computed
          LM_Stru%Flag_weight = .True.

        end if ! Weight not flagged in LM

        ! Compute chi2
        LM_Stru%Chisq_og = sum(LM_Stru%Residual* &
                               LM_Stru%Residual* &
                               LM_Stru%Weight)

      end if ! Type of inversion

      ! Update chi2 in LM structure
      LM_Stru%Chisq = LM_Stru%Chisq_og

      return

      end subroutine Merit_function

!#####################################################################
!#####################################################################
!#####################################################################

      !> Predict the improvement for the calculated step\n
      !!         AA(double(:,:)): System matrix\n
      !!           gg(double(:)): Independent term\n
      !!           xx(double(:)): Solution\n
      !!            pred(double): Predicted reduction
      subroutine predict_improvement(AA,gg,xx,pred)

      ! I/O

      double precision, intent(out):: pred
      double precision, dimension(:), intent(in):: gg
      double precision, dimension(:), intent(in):: xx
      double precision, dimension(:,:), intent(in):: AA

      ! Local

      integer:: NN,ii,jj


      ! Size of vector
      NN = size(xx)

      ! Initialize prediction
      pred = 0d0

      ! For each row
      do ii=1,NN

        ! For each column
        do jj=1,NN

          ! Add to prediction
          pred = pred - 0.5d0*xx(ii)*AA(ii,jj)*xx(jj)

        end do ! Columns

        ! Add to prediction (gg is defined as negative)
        pred = pred + gg(ii)*xx(ii)

      end do ! Rows

      end subroutine predict_improvement

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get angle between the two vectors and their norm\n
      !!    aa(double(:)): One of the vectors\n
      !!    bb(double(:)): The other vector\n
      !!      NN(integer): Size of the vectors\n
      !!   cosine(double): Cosine of the angle between the vectors\n
      !!    anorm(double): Norm of the first vector\n
      !!    bnorm(double): Norm of the second vector
      subroutine angle_and_norm(aa,bb,NN,cosine,anorm,bnorm)

      ! I/O

      integer, intent(in):: NN
      double precision, intent(out):: cosine,anorm,bnorm
      double precision, dimension(:), intent(in):: aa,bb

      ! Local

      integer:: ii


      ! Initialize
      cosine = 0d0
      anorm = 0d0
      bnorm = 0d0

      ! For each element
      do ii=1,NN

        ! Numerator
        cosine = cosine + aa(ii)*bb(ii)

        ! Norms
        anorm = anorm + aa(ii)*aa(ii)
        bnorm = bnorm + bb(ii)*bb(ii)

      end do

      ! Norm
      anorm = sqrt(anorm)
      bnorm = sqrt(bnorm)

      ! Sanity
      if (anorm.le.0d0) anorm = 1d0
      if (bnorm.le.0d0) bnorm = 1d0

      ! Cosine
      cosine = cosine/anorm/bnorm

      end subroutine angle_and_norm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get percentile and maximum of a given vector ignoring zeros\n
      !!      arr(double(:)): Given vector\n
      !!         nn(integer): Original number of elements\n
      !!         nc(integer): Number of non-zero elements\n
      !!  percentile(double): Percentile to take\n
      !!        perc(double): Percentile to return\n
      !!     maximum(double): Maximum of vector
      subroutine get_percentile_and_maximum(arr,nn,nc,percentile, &
                                            perc,maximum)

      ! I/O

      integer, intent(in):: nn,nc
      double precision, intent(in):: percentile
      double precision, intent(out):: perc,maximum
      double precision, dimension(nn), intent(in):: arr

      ! Local

      integer:: n0,ii

      double precision, dimension(nn):: arr2


      ! Copy array
      arr2 = arr

      ! Order array
      call QsortC(arr2)

      ! Number of zeros
      n0 = nn - nc

      ! Move to remove zeros
      if (n0.gt.0) arr2(1:nc) = arr2(n0+1:nn)

      ! Maximum is the last one
      maximum = arr2(nc)

      ! Get percentile
      ii = nint(dble(nc)*percentile)
      perc = arr2(ii)

      end subroutine get_percentile_and_maximum

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the Hessian from the Jacobian with given weights\n
      !!      thermal(logical): If thermal inversion\n
      !!        Regul(logical): If there are regularizations\n
      !!  LM_Stru(LMFIT_class): Structure with data for the
      !!                        Levenberg–Marquardt\n
      !!       WW(double(:,:)): Weights to use to compute Hessian and
      !!                        Jacobian\n
      !!         gg(double(:)): Computed Jacobian\n
      !!         DD(double(:)): Diagonal of the Hessian\n
      !!       AA(double(:,:)): Computed Hessian
      subroutine adjust_weights_hessian(thermal,Regul,LM_Stru,WW, &
                                        gg,DD,AA)

      ! I/O

      type(LMFIT_class), intent(inout):: LM_Stru
      logical, intent(in):: thermal,Regul
      double precision, dimension(:,:), &
                        allocatable, intent(in):: WW
      double precision, dimension(:), intent(out):: gg,DD
      double precision, dimension(:,:), intent(out):: AA

      ! Local

      integer:: i, j


      ! Initialize
      gg = 0d0
      AA = 0d0

      ! If thermal
      if (thermal) then

        ! Row (Hessian column)
        do i=1,LM_Stru%Num

          ! Jacobian vector
          gg(i) = -1d0*sum(LM_Stru%JacobianI(:,i)* &
                           LM_Stru%ResidualI*WW(0,:))

          ! Hessian matrix diagonal
          AA(i,i) = sum(LM_Stru%JacobianI(:,i)* &
                        LM_Stru%JacobianI(:,i)*WW(0,:))

          ! Hessian row
          do j=i+1,LM_Stru%Num

            ! Hessian matrix
            AA(j,i) = sum(LM_Stru%JacobianI(:,i)* &
                          LM_Stru%JacobianI(:,j)*WW(0,:))

            ! Symmetric matrix
            AA(i,j) = AA(j,i)

          end do ! Hessian row
        end do ! Row (Hessian column)

      ! Polarization inversion
      else

        ! Row (Hessian column)
        do i=1,LM_Stru%Num

          ! Jacobian vector
          gg(i) = -1d0*sum(LM_Stru%Jacobian(:,:,i)* &
                           LM_Stru%Residual*WW)

          ! Hessian matrix diagonal
          AA(i,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                        LM_Stru%Jacobian(:,:,i)*WW)

          ! Hessian row
          do j=i+1,LM_Stru%Num

            ! Hessian matrix off-diagonal
            AA(j,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                          LM_Stru%Jacobian(:,:,j)*WW)

            ! Symmetric matrix
            AA(i,j) = AA(j,i)

          end do ! Hessian row
        end do ! Row (Hessian column)

      end if ! Type of inversion

      ! Regularization
      if (Regul) then

        ! Add to Hessian
        AA = AA + LM_Stru%Rgl%Regul_H*LM_Stru%Rgl%Ratio

        ! Add to vector
        gg = gg + LM_Stru%Rgl%Regul_F*LM_Stru%Rgl%Ratio

      end if ! Regularizing

      ! Diagonal
      do i=1,LM_Stru%Num
        DD(i) = AA(i,i)
      end do

      return

      end subroutine adjust_weights_hessian

!#####################################################################
!#####################################################################
!#####################################################################

      !> Modify weights to enhance larger differences in the fit for
      !! a given set of parameters\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!          thermal(logical): If thermal inversion\n
      !!        param_perc(double): Percentile to calculate
      !!                            percentile value\n
      !!      param_ff_max(double): Maximum enhancement of a weight\n
      !!    param_ff_contr(double): Contribution fraction to which
      !!                            apply enhancement\n
      !!            Weight(double): Modified weights to return
      subroutine adjust_weights_trial(LM_Stru,Inf_Stokes,thermal, &
                                      param_perc,param_ff_max, &
                                      param_ff_contr,Weight)

      ! I/O

      type(LMFIT_class), intent(in):: LM_Stru
      type(Stokes_class), intent(in):: Inf_Stokes
      logical, intent(in):: thermal
      double precision, intent(in):: param_perc,param_ff_max
      double precision, intent(in):: param_ff_contr
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: Weight

      ! Local

      integer:: ii,jj,kk,ntotal
      integer, dimension(1):: tmp

      double precision:: S0,perc,maximum,ff,cResi,tResi
      double precision, dimension(:), allocatable:: WW
      double precision, dimension(:), allocatable:: Resi

      ! Thermal
      if (thermal) then

        ! Allocate
        allocate(Resi(Inf_Stokes%Num_Wavelength))
        allocate(WW(Inf_Stokes%Num_Wavelength))

        ! Square residual
        Resi = LM_Stru%ResidualI*LM_Stru%ResidualI

        ! If using sigma
        if (Inf_Stokes%Sigma_Flag) then

          ! Get residual
          Resi = Resi/(Inf_Stokes%Sigma_W(0,:)* &
                       Inf_Stokes%Sigma_W(0,:))

        end if ! Using sigma

        ! Remove neglected
        ntotal = 0
        do ii=1,Inf_Stokes%Num_Wavelength
          if (Inf_Stokes%weight(0,ii).le.0d0) then
            Resi(ii) = 0d0
          else
            ntotal = ntotal + 1
          end if
        end do

      ! Non-thermal
      else

        ! Allocate
        allocate(Resi(Inf_Stokes%Num_Wavelength*4))
        allocate(WW(Inf_Stokes%Num_Wavelength*4))

        ! For each Stokes
        do ii=0,3

          ! Square residual
          Resi(ii*Inf_Stokes%Num_Wavelength+1: &
               (ii+1)*Inf_Stokes%Num_Wavelength) = &
                                             LM_Stru%Residual(ii,:)* &
                                             LM_Stru%Residual(ii,:)

        end do ! For each Stokes

        ! If using sigma
        if (Inf_Stokes%Sigma_Flag) then

          ! For each Stokes
          do ii=0,3

            ! Get residual
            Resi(ii*Inf_Stokes%Num_Wavelength+1: &
                 (ii+1)*Inf_Stokes%Num_Wavelength) = &
                             Resi(ii*Inf_Stokes%Num_Wavelength+1: &
                                  (ii+1)*Inf_Stokes%Num_Wavelength)/ &
                                          (Inf_Stokes%Sigma_W(ii,:)* &
                                           Inf_Stokes%Sigma_W(ii,:))
          end do ! Stokes parameters

        end if ! Using sigma

        ! Remove neglected
        ntotal = 0
        kk = 0
        do ii=0,3
          do jj=1,Inf_Stokes%Num_Wavelength
            kk = kk + 1
            if (Inf_Stokes%weight(ii,jj).le.0d0) then
              Resi(kk) = 0d0
            else
              ntotal = ntotal + 1
            end if
          end do
        end do

      end if ! Thermal inversion

      ! Get median and maximum
      call get_percentile_and_maximum(Resi, &
                                      Inf_Stokes%Num_Wavelength, &
                                      ntotal,param_perc, &
                                      perc,maximum)

      ! Enhancement factor
      ff = maximum/perc
      ff = min(ff,param_ff_max)

      ! Thermal
      if (thermal) then

        ! Get actual contributions to chi^2
        Resi = LM_Stru%ResidualI*LM_Stru%ResidualI*LM_Stru%WeightI

        ! L2 norm of old weights
        S0 = sum(LM_Stru%WeightI)

        ! Target and initial residual
        tResi = sum(Resi)*param_ff_contr
        cResi = 0d0

        ! Till enough contribution enhanced
        do while (cResi.lt.tResi)

          ! Location of maximum
          tmp = maxloc(Resi)
          ii = tmp(1)

          ! Process
          Weight(0,ii) = Weight(0,ii)*ff
          cResi = cResi + Resi(ii)
          Resi(ii) = 0d0

        end do

        ! New weights
        Weight = Weight*S0/sum(LM_Stru%WeightI)

      ! Non-thermal
      else

        ! For each Stokes
        do ii=0,3

          ! Get actual contributions to chi^2
          Resi(ii*Inf_Stokes%Num_Wavelength+1: &
               (ii+1)*Inf_Stokes%Num_Wavelength) = &
                                           LM_Stru%Residual(ii,:)* &
                                           LM_Stru%Residual(ii,:)* &
                                           LM_Stru%Weight(ii,:)

        end do ! For each Stokes

        ! L2 norm of old weights
        S0 = sum(LM_Stru%Weight)

        ! Target and initial residual
        tResi = sum(Resi)*param_ff_contr
        cResi = 0d0

        ! Till enough contribution enhanced
        do while (cResi.lt.tResi)

          ! Location of maximum
          tmp = maxloc(Resi)
          kk = tmp(1)

          ! Get other indexes
          ii = (kk-1)/Inf_Stokes%Num_Wavelength
          jj = kk - ii*Inf_Stokes%Num_Wavelength

          ! Process
          Weight(ii,jj) = Weight(ii,jj)*ff
          cResi = cResi + Resi(kk)
          Resi(kk) = 0d0

        end do

        ! New weights
        Weight = Weight*S0/sum(LM_Stru%Weight)

      end if ! Thermal

      return

      end subroutine adjust_weights_trial

!#####################################################################
!#####################################################################
!#####################################################################

      !> Modify weights to enhance larger differences in the fit\n
      !!      LM_Stru(LMFIT_class): Structure with data for the
      !!                            Levenberg–Marquardt\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!        Input(Input_class): Structure with configuration data
      subroutine adjust_weights(LM_Stru,Inf_Nodes,Inf_Stokes,Input)

      ! I/O

      type(Input_class), intent(in):: Input
      type(Nodes_class), intent(in):: Inf_Nodes
      type(Stokes_class), intent(in):: Inf_Stokes
      type(LMFIT_class), intent(inout):: LM_Stru

      ! Local

      type(Nodes_class):: Inf_Nodes_tmp

      logical:: thermal

      integer:: i1,i2,i3,i4,ii

      double precision:: step,pred,cost,cosm,normg,normx
      double precision, dimension(:), allocatable:: gg,DD,xx
      double precision, dimension(:,:), allocatable:: HH,AA,Weight

      ! Parameters
      integer, parameter:: n1=4
      double precision, dimension(n1), &
             parameter:: c_perc = (/ 0.15d0,0.25d0,0.3d0,0.5d0 /)
      integer, parameter:: n2=4
      double precision, dimension(n2), &
             parameter:: c_ff_max = (/ 5d0,10d0,20d0,30d0 /)
      integer, parameter:: n3=4
      double precision, dimension(n3), &
             parameter:: c_ff_contr = (/ 0.15d0,0.25d0,0.4d0,0.6d0 /)
      integer, parameter:: n4=4
      double precision, dimension(n4), &
             parameter:: c_lambda = (/ 5d0,50d0,500d0,2500d0 /)


      ! If thermal
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Allocate weight
        allocate(Weight(0:0,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%WeightI_mod(Inf_Stokes%Num_Wavelength))
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%WeightI_mod)

        ! Flag thermal
        thermal = .True.

      ! Non-thermal
      else

        ! Allocate weight
        allocate(Weight(0:3,Inf_Stokes%Num_Wavelength))
        allocate(LM_Stru%Weight_mod(0:3,Inf_Stokes%Num_Wavelength))
        MRAMc = MRAMc + 1d-6*sizeof(LM_Stru%Weight_mod)

        ! Flag thermal
        thermal = .False.

      end if

      ! Allocate space for new Hessian
      allocate(gg(LM_Stru%Num),DD(LM_Stru%Num),xx(LM_Stru%Num))
      allocate(HH(LM_Stru%Num,LM_Stru%Num))
      allocate(AA(LM_Stru%Num,LM_Stru%Num))

      ! Initialize
      LM_Stru%perc = 0.25d0
      LM_Stru%ff_max = 20d0
      LM_Stru%ff_contr = 0.35d0
      step = -1d0

      ! For each percentile
      do i1=1,n1

        ! For each maximum factor
        do i2=1,n2

          ! For each chi2 contribution
          do i3=1,n3

            ! Copy originals
            if (thermal) then
              Weight(0,:) = LM_Stru%WeightI
            else
              Weight = LM_Stru%Weight
            end if

            ! Trial
            call adjust_weights_trial(LM_Stru,Inf_Stokes,thermal, &
                                      c_perc(i1),c_ff_max(i2), &
                                      c_ff_contr(i3),Weight)

            ! Compute Hessian
            call adjust_weights_hessian(thermal, &
                                        Inf_Nodes%Regul_flag, &
                                        LM_Stru,Weight,gg,DD,HH)

            ! For each lambda
            do i4=1,n4

              ! Get A matrix
              AA = HH
              do ii=1,LM_Stru%Num
                AA(ii,ii) = AA(ii,ii) + c_lambda(i4)*DD(ii)
              end do

              ! Copy current Nodes
              Inf_Nodes_tmp = Inf_Nodes

              ! SVD solution
              call SVD_Solve(AA,gg,xx,LM_Stru%Num, &
                             Inf_Nodes_tmp,Input%SVD_Type, &
                             .False.)

              ! Failed
              if (laborted) then
                laborted = .False.
                cycle
              end if

              ! Angle model
              call angle_and_norm(gg,LM_Stru%Jacfvec, &
                                  LM_Stru%Num, &
                                  cosm,normg,normx)

              ! Sanity
              if (cosm.lt.0.0d0) cycle

              ! Angle consistency
              call angle_and_norm(gg,xx,LM_Stru%Num, &
                                  cost,normg,normx)

              ! Sanity
              if (cost.lt.0.3d0) cycle
              if (normx.gt.0.2d0) cycle

              ! Predict improvement and check step
              call predict_improvement(HH,gg,xx,pred)

              ! Score
              pred = pred*(cosm*cosm)/(1d0 + normx)

              ! If better score
              if (pred.gt.step) then

                ! Initialize
                LM_Stru%perc = c_perc(i1)
                LM_Stru%ff_max = c_ff_max(i2)
                LM_Stru%ff_contr = c_ff_contr(i3)
                step = pred

              end if ! Better score

            end do ! Lambdas
          end do ! ff_contr
        end do ! ff_max
      end do ! perc

      ! Copy originals
      if (thermal) then
        Weight(0,:) = LM_Stru%WeightI
      else
        Weight = LM_Stru%Weight
      end if

      ! Final
      call adjust_weights_trial(LM_Stru,Inf_Stokes,thermal, &
                                LM_Stru%perc,LM_Stru%ff_max, &
                                LM_Stru%ff_contr,Weight)

      ! Copy
      if (thermal) then
        LM_Stru%WeightI_mod = Weight(0,:)
      else
        LM_Stru%Weight_mod = Weight
      end if

      ! Deallocate
      deallocate(Weight)

      return

      end subroutine adjust_weights

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the Jacobian for every inversion parameter\n
      !!      Input(Input_class): Structure with configuration data\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!    Atomb(Atom_class(:)): Structures with atomic data for
      !!                          background atoms\n
      !!       Mol(Mol_class(:)): Structures with molecular data\n
      !!    Geom(Geometry_class): Structure with geometric data\n
      !!   GeomI(Geometry_class): Structure with geometric data for
      !!                          the intensity problem\n
      !!      Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                          and J-symbols\n
      !!   Frec(Frequency_class): Structure with frequency data\n
      !!      fudge(fudge_class): Structure with fudge data\n
      !!    kurucz(kurucz_class): Structure with Kurucz line data\n
      !!         MPID(MPI_class): Structure with MPI data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Bfield(Bfield_class): Structure with magnetic field data\n
      !!  Inf_Nodes(Nodes_class): Structure with inversion node data\n
      !!     Sol(Solution_class): Structure with the frequency and
      !!                          synthetic Stokes parameters in the
      !!                          frequency range of the inverted
      !!                          data\n
      !!  SolF(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one\n
      !!    LM_Stru(LMFIT_class): Structure with data for the
      !!                          Levenberg–Marquardt
      subroutine Jacobian_Compute(Input,Atom,Atomb,Mol,Geom,GeomI, &
                                  Flgsg,Frec,fudge,kurucz,MPID,Atmo, &
                                  Bfield,Inf_Nodes,Sol,SolF,LM_Stru)

      ! I/O

      type(Input_class), intent(inout):: Input
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
      type(Nodes_class), intent(inout)::Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      type(LMFIT_class), intent(inout):: LM_Stru

      ! Local

      integer:: i,j,Indx_Para


      ! Master
      if (pid.eq.0) then

        ! Verbose
        umsg = ' - Compute the Jacobian'
        call verboseI(3)

      end if ! Master


      ! If inversion is only thermal
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Initialize j index
        j = 0

        ! For each thermal node
        do i=1,Inf_Nodes%Num_fit

          ! Get parameter
          Indx_Para = Inf_Nodes%Inf_Inv(1,i)

          ! If not thermal, skip
          if (Indx_Para.lt.Inf_Nodes%index_f.or. &
              Indx_Para.gt.Inf_Nodes%index_Pg) cycle

          ! Advance index
          j = j + 1

          ! Get response function
          call RF_Thermo(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                         fudge,kurucz,MPID,Atmo,Bfield, &
                         Inf_Nodes,i,LM_Stru%JacobianI(:,j), &
                         Sol,SolF,Input)

        end do ! Thermal nodes

      ! If inversion is magnetic
      else if (Inf_Nodes%Nodes_Type.eq.1) then

        ! Initialize j index
        j = 0

        ! For each magnetic node
        do i=1,Inf_Nodes%Num_fit

          ! Get parameter
          Indx_Para = Inf_Nodes%Inf_Inv(1,i)

          ! If not magnetic, skip
          if (Indx_Para.gt.Inf_Nodes%index_f) cycle

          ! Advance index
          j = j + 1

          ! Get response function
          call RF_Mag(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                      kurucz,MPID,Atmo,Bfield,Inf_Nodes,i, &
                      LM_Stru%Jacobian(:,:,i),Sol,SolF,Input)

        end do ! Magnetic nodes

      ! If inverting all
      else if (Inf_Nodes%Nodes_Type.eq.2) then

        ! For all nodes
        do i=1,Inf_Nodes%Num_Fit

          ! Get response function
          call RF_ALL(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                      kurucz,MPID,Atmo,Bfield,Inf_Nodes,i, &
                      LM_Stru%Jacobian(:,:,i),Sol,SolF,Input)

        end do ! All nodes

      end if ! Type of inversion

      ! Flag as true Jacobian
      LM_Stru%Flag_Jac = .True.

      return

      end subroutine Jacobian_Compute

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the Hessian from the Jacobian\n
      !!   Nodes_Type(integer): Type of inversion\n
      !!  LM_Stru(LMFIT_class): Structure with data for the
      !!                        Levenberg–Marquardt\n
      !!       RDmode(logical): If we need to use modified weights
      subroutine Hessian_Compute(Nodes_Type,LM_Stru,RDmode)

      ! I/O

      type(LMFIT_class), intent(inout):: LM_Stru
      logical, intent(in):: RDmode
      integer, intent(in):: Nodes_Type

      ! Local

      integer:: i, j


      ! Initialize
      LM_Stru%Jacfvec_og = 0d0
      LM_Stru%Hessian_og = 0d0

      ! Reduced Mode
      if (RDmode) then

        ! If thermal inversion (only intensity)
        if (Nodes_Type.eq.0) then

          ! Row (Hessian column)
          do i=1,LM_Stru%Num

            ! Jacobian vector
            LM_Stru%Jacfvec_og(i) = -1d0*sum(LM_Stru%JacobianI(:,i)* &
                                             LM_Stru%ResidualI* &
                                             LM_Stru%WeightI_mod)

            ! Hessian matrix diagonal
            LM_Stru%Hessian_og(i,i) = sum(LM_Stru%JacobianI(:,i)* &
                                          LM_Stru%JacobianI(:,i)* &
                                          LM_Stru%WeightI_mod)

            ! Hessian row
            do j=i+1,LM_Stru%Num

              ! Hessian matrix
              LM_Stru%Hessian_og(j,i) = sum(LM_Stru%JacobianI(:,i)* &
                                            LM_Stru%JacobianI(:,j)* &
                                            LM_Stru%WeightI_mod)

              ! Symmetric matrix
              LM_Stru%Hessian_og(i,j) = LM_Stru%Hessian_og(j,i)

            end do ! Hessian row
          end do ! Row (Hessian column)

        ! Polarization inversion
        else

          ! Row (Hessian column)
          do i=1,LM_Stru%Num

            ! Jacobian vector
            LM_Stru%Jacfvec_og(i) = -1d0*sum(LM_Stru%Jacobian(:,:,i)*&
                                             LM_Stru%Residual* &
                                             LM_Stru%Weight_mod)

            ! Hessian matrix diagonal
            LM_Stru%Hessian_og(i,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                          LM_Stru%Jacobian(:,:,i)* &
                                          LM_Stru%Weight_mod)

            ! Hessian row
            do j=i+1,LM_Stru%Num

              ! Hessian matrix off-diagonal
              LM_Stru%Hessian_og(j,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                            LM_Stru%Jacobian(:,:,j)* &
                                            LM_Stru%Weight_mod)

              ! Symmetric matrix
              LM_Stru%Hessian_og(i,j) = LM_Stru%Hessian_og(j,i)

            end do ! Hessian row
          end do ! Row (Hessian column)

        end if ! Type of inversion

      ! Normal mode
      else

        ! If thermal inversion (only intensity)
        if (Nodes_Type.eq.0) then

          ! Row (Hessian column)
          do i=1,LM_Stru%Num

            ! Jacobian vector
            LM_Stru%Jacfvec_og(i) = -1d0*sum(LM_Stru%JacobianI(:,i)* &
                                             LM_Stru%ResidualI* &
                                             LM_Stru%WeightI)

            ! Hessian matrix diagonal
            LM_Stru%Hessian_og(i,i) = sum(LM_Stru%JacobianI(:,i)* &
                                          LM_Stru%JacobianI(:,i)* &
                                          LM_Stru%WeightI)

            ! Hessian row
            do j=i+1,LM_Stru%Num

              ! Hessian matrix
              LM_Stru%Hessian_og(j,i) = sum(LM_Stru%JacobianI(:,i)* &
                                            LM_Stru%JacobianI(:,j)* &
                                            LM_Stru%WeightI)

              ! Symmetric matrix
              LM_Stru%Hessian_og(i,j) = LM_Stru%Hessian_og(j,i)

            end do ! Hessian row
          end do ! Row (Hessian column)

        ! Polarization inversion
        else

          ! Row (Hessian column)
          do i=1,LM_Stru%Num

            ! Jacobian vector
            LM_Stru%Jacfvec_og(i) = -1d0*sum(LM_Stru%Jacobian(:,:,i)*&
                                             LM_Stru%Residual* &
                                             LM_Stru%Weight)

            ! Hessian matrix diagonal
            LM_Stru%Hessian_og(i,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                          LM_Stru%Jacobian(:,:,i)* &
                                          LM_Stru%Weight)

            ! Hessian row
            do j=i+1,LM_Stru%Num

              ! Hessian matrix off-diagonal
              LM_Stru%Hessian_og(j,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                            LM_Stru%Jacobian(:,:,j)* &
                                            LM_Stru%Weight)

              ! Symmetric matrix
              LM_Stru%Hessian_og(i,j) = LM_Stru%Hessian_og(j,i)

            end do ! Hessian row
          end do ! Row (Hessian column)

        end if ! Type of inversion
      end if ! Reduced or normal mode

      ! Copy results
      LM_Stru%Jacfvec = LM_Stru%Jacfvec_og
      LM_Stru%Hessian = LM_Stru%Hessian_og

      return

      end subroutine Hessian_Compute

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rank 1 update of the Jacobian following Broyden's method\n
      !!  Stokes_Old(double(:,:)): Best fit Stokes parameters\n
      !!  Stokes_New(double(:,:)): Backtracking proposed Stokes
      !!                           parameters\n
      !!        Num_Wave(integer): Number of wavelengths\n
      !!      Solution(double(:)): Last solution from the Hessian\n
      !!   Inf_Nodes(Nodes_class): Structure with inversion node
      !!                           data\n
      !!     LM_Stru(LMFIT_class): Structure with data for the
      !!                           Levenberg–Marquardt
      subroutine Broyden_Rank1(Stokes_Old,Stokes_New,Num_wave, &
                               Solution,Inf_Nodes,LM_Stru)

      ! I/O

      type(Nodes_class), intent(in):: Inf_Nodes
      type(LMFIT_class), intent(inout):: LM_Stru
      integer, intent(in):: Num_wave
      double precision, dimension(:), intent(in):: Solution
      double precision, dimension(:,:), intent(in):: Stokes_Old
      double precision, dimension(:,:), intent(in):: Stokes_New

      ! Local

      integer:: i,j,k

      double precision:: Suma
      double precision, dimension(:), allocatable:: StokesI
      double precision, dimension(:,:), allocatable:: Stokes

      ! Get quadratic sum of solution
      Suma = 1d0/sum(Solution(1:LM_Stru%Num)* &
                     Solution(1:LM_Stru%Num))

      ! If intensity only inversion
      if (Inf_Nodes%Nodes_Type.eq.0) then

        ! Allocate space for Stokes parameters
        allocate(StokesI(Num_wave))

        ! For each wavelength
        do j=1,Num_wave

          ! Initialize to difference between new and old Stokes
          StokesI(j) = Stokes_New(1,j) - Stokes_Old(1,j)

          ! Add contribution from Jacobian
          StokesI(j) = StokesI(j) + &
                       sum(LM_Stru%JacobianI(j,:)*Solution)

        end do  ! Wavelenghts

        ! For each variable/node
        do k=1,LM_Stru%Num

          ! Modify Jacobian
          LM_Stru%JacobianI(:,k) = StokesI*Solution(k)*Suma

        end do ! Variables/nodes

        ! Free memory
        deallocate(StokesI)

      ! Polarization
      else

        ! Allocate space for Stokes parameters
        allocate(Stokes(0:3,Num_wave))

        ! Initialize Stokes to the difference
        Stokes = Stokes_New - Stokes_Old

        ! For each wavelength
        do j=1,Num_wave
          ! For each Stokes parameter
          do i=0,3

            ! Add contribution from Jacobian
            Stokes(i,j) = Stokes(i,j) + &
                          sum(LM_Stru%Jacobian(i,j,:)*Solution)

          end do ! Stokes parameters
        end do ! Wavelength

        ! For each variable/node
        do k=1,LM_Stru%Num

          ! Modify Jacobian
          LM_Stru%Jacobian(:,:,k) = Stokes*Solution(k)*Suma

        end do ! Variables/nodes

        ! Free memory
        deallocate(Stokes)

      end if ! Type of inversion

      ! Flag Jacobian as false
      LM_Stru%Flag_Jac = .False.

      return

      end subroutine Broyden_Rank1

!#####################################################################
!#####################################################################
!#####################################################################

      end module jacobian_mod
