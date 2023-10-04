      !> Manages the jacobian for the LM
      module jacobian_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/24/2023
!  Last version:
!     09/08/2023 V3.0.7
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/08/2023:    V3.0.7 - Verbosity update (TdPA)
!
!     08/11/2023:    V3.0.6 - The inversion weights are now fully
!                             wavelength dependent (TdPA)
!
!     07/31/2023:    V3.0.5 - Change the verbosity level in the
!                             inversion (HL)
!
!     07/03/2023:    V3.0.4 - Made changes to run over the relevant
!                             variable indexes, as the addition of
!                             diffuse light made it non trivial (TdPA)
!
!     06/12/2023:    V3.0.3 - Rename the variable Inf_File (HL)
!
!     04/11/2023:    V3.0.2 - Update the weiths for multi-wavelength
!                             ranges (HL)
!
!     03/15/2023:    V3.0.1 - Removed some commented lines (TdPA)
!                           - The Blos variables are in the same
!                             structure than the polar ones (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/24/2023:    V0.0.0 - Started from 12/05/2020
!                             TIC@jacobian_mod.f90 revision (TdPA)
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
!    Merit_function:
!      Compute the merit function from the Stokes profiles
!
!    Jacobian_Compute:
!      Compute the Jacobian of with respect to the parameters
!
!    Hessian_Compute:
!      Compute the Hessian from the Jacobian
!
!    Broyden_Rank1:
!      Rank 1 update of the Jacobian following Broyden's method
!
!#####################################################################
!#####################################################################
!#####################################################################

      use commons_mod
      use rf_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the merit function from the Stokes profiles
      !!  Inf_Stokes(Stokes_class): Structure with Stokes parameters
      !!                            data\n
      !!       Stokes(double(:,:)): Emerging\n
      !!    Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!      LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                            other LM quantities
      subroutine Merit_function(Inf_Stokes,Stokes_out,Nodes_Type, &
                                LM_Stru)

      ! IO
      type(Stokes_class), intent(in):: Inf_Stokes
      type(LMFIT_class), intent(inout):: LM_Stru
      integer, intent(in):: Nodes_Type
      double precision, dimension(:,:), allocatable, intent(in):: &
                                                            Stokes_out


      ! Initialize chi
      LM_Stru%Chisq_og = 0d0

      ! If thermal
      if (Nodes_Type.eq.0) then

        ! Compute difference
        LM_Stru%ResidualI = Stokes_out(0,:) - &
                            Inf_Stokes%Stokes_Ob(0,:)

        ! If weight not flagged in LM
        if (.not.LM_Stru%Flag_weight) then

          ! If wavelength dependent sigma
          if (Inf_Stokes%Sigma_Flag) then

            ! Compute weight
            LM_Stru%WeightI = Inf_Stokes%weight(0,:)* &
                              Inf_Stokes%weight(0,:)/ &
                              (Inf_Stokes%Sigma_W(0,:)* &
                               Inf_Stokes%sigma_W(0,:))/ &
                              dble(Inf_Stokes%Num_freedomI)

          ! No wavelength dependent sigma
          else

            ! Weight
            LM_Stru%WeightI = Inf_Stokes%weight(0,:)* &
                              Inf_Stokes%weight(0,:)/ &
                              dble(Inf_Stokes%Num_freedomI)

          end if ! Frequency dependent sigma

          ! Flag weight as computed
          LM_Stru%Flag_weight = .True.

        end if ! Flagged weight in LM

        ! Get chi^2
        LM_Stru%Chisq_og = sum(LM_Stru%ResidualI* &
                               LM_Stru%ResidualI* &
                               LM_STru%WeightI)
      ! Non-thermal
      else

        ! Compute difference
        LM_Stru%Residual = Stokes_out - &
                           Inf_Stokes%Stokes_Ob

        ! If weight not flagged in LM
        if (.not.LM_Stru%Flag_weight) then

          ! Weight contribution
          LM_Stru%weight = Inf_Stokes%Weight*Inf_Stokes%Weight/ &
                           dble(Inf_Stokes%Num_freedom)

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

      ! Update chi2
      LM_Stru%Chisq = LM_Stru%Chisq_og

      return

      end subroutine Merit_function

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the Jacobian of with respect to the parameters\n
      !!         Input(Input_class): Structure with settings data\n
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
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!     SolF(Solution_F_class): Class with the full RT
      !!                             solution\n
      !!       LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                             other LM quantities
      subroutine Jacobian_Compute(Input,Atom,Atomb,Mol,Geom, &
                                  GeomI,Flgsg,Frec,fudge,kurucz, &
                                  MPID,Atmo,Bfield,Inf_Nodes,Sol, &
                                  SolF,LM_Stru)

      ! I/O
      type(Input_class), intent(inout):: Input
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
      type(Nodes_class), intent(inout)::Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      type(LMFIT_class), intent(inout):: LM_Stru

      ! Local
      integer:: i, j, Indx_Para


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
      !!  Nodes_Type(integer): Type of inversion\n
      !! LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                       other LM quantities
      subroutine Hessian_Compute(Nodes_Type,LM_Stru)

      ! IO
      type(LMFIT_class), intent(inout):: LM_Stru
      integer, intent(in):: Nodes_Type

      ! Local
      integer:: i, j


      ! Initialize
      LM_Stru%Jacfvec_og = 0d0
      LM_Stru%Hessian_og = 0d0

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
          LM_Stru%Jacfvec_og(i) = -1d0*sum(LM_Stru%Jacobian(:,:,i)* &
                                           LM_Stru%Residual* &
                                           LM_Stru%Weight)

          ! Hessian matrix diagonal
          LM_Stru%Hessian_og(i,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                        LM_Stru%Jacobian(:,:,i)* &
                                        LM_Stru%Weight)

          ! Hessian row
          do j=i+1,LM_Stru%Num

            ! Hessian matrix diagonal
            LM_Stru%Hessian_og(j,i) = sum(LM_Stru%Jacobian(:,:,i)* &
                                          LM_Stru%Jacobian(:,:,j)* &
                                          LM_Stru%Weight)

            ! Symmetric matrix
            LM_Stru%Hessian_og(i,j) = LM_Stru%Hessian_og(j,i)

          end do ! Hessian row
        end do ! Row (Hessian column)

      end if ! Type of inversion

      ! Copy results
      LM_Stru%Jacfvec = LM_Stru%Jacfvec_og
      LM_Stru%Hessian = LM_Stru%Hessian_og

      return

      end subroutine Hessian_Compute

!#####################################################################
!#####################################################################
!#####################################################################

      !> Rank 1 update of the Jacobian following Broyden's method\n
      !! Stokes_Old(double(:,:)): Stokes parameters for best fit\n
      !! Stokes_New(double(:,:)): Stokes parameters proposed by
      !!                          the backtracking\n
      !!       Num_Wave(integer): Number of wavelengths\n
      !!     Solution(double(:)): Last solution of the Hessian\n
      !!  Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!      LM_Stru(LMFIT_class): Structure with Jacobian and
      !!                            other LM quantities
      subroutine Broyden_Rank1(Stokes_Old,Stokes_New,Num_wave, &
                               Solution,Inf_Nodes,LM_Stru)

      ! IO
      type(Nodes_class), intent(in):: Inf_Nodes
      type(LMFIT_class), intent(inout):: LM_Stru
      integer:: Num_wave
      double precision, dimension(:), intent(in):: Solution
      double precision, dimension(:,:), intent(in):: Stokes_Old
      double precision, dimension(:,:), intent(in):: Stokes_New

      ! Local
      integer:: i, j, k

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
