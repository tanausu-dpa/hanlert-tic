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
!  Merit_function
!    Compute the L2 merit function of Stokes profiles
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
      use rf_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the L2 merit function of Stokes profiles\n
      !!  Inf_Stokes(Stokes_class): Structure with inversion Stokes
      !!                            parameters data\n
      !!   Stokes_out(double(:,:)): Stokes parameters\n
      !!    Inf_Nodes(Nodes_class): Structure with inversion node
      !!                            data\n
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
                               Inf_Stokes%sigma_W(0,:))

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
                               LM_STru%WeightI)
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
      type(Atom_class), dimension(:), intent(inout):: Atom
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
      !!                        Levenberg–Marquardt
      subroutine Hessian_Compute(Nodes_Type,LM_Stru)

      ! I/O

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
