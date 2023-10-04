      !> Compute response functions
      module rf_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Hao Li (IAC)
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!  Start:
!     02/27/2023
!  Last version:
!     09/08/2023 V3.0.5
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     09/08/2023:    V3.0.5 - Verbosity update (TdPA)
!
!     07/31/2023:    V3.0.4 - Change the verbosity level in the
!                             inversion (HL)
!
!     07/03/2023:    V3.0.3 - Renamed Inf_File to Input (TdPA)
!                           - Removed checks in the parameter index
!                             in the RF functions, as the calling
!                             function takes care of calling the
!                             correct ones (TdPA)
!                           - RF routines no longer have a saying
!                             in the synthesis configuration, they
!                             just call HanleRTTIC indicating that
!                             a RF is being calculated (TdPA)
!                           - The variable factor is no longer
!                             hard-coded. Instead, the user can
!                             specify a minimum relative perturbation
!                             and factor will be adjusted to this
!                             request (TdPA)
!                           - Added RF_fInt and RF_fPol to
!                             compute the diffuse light factor RF
!                             without calling a synthesis (TdPA)
!                           - Added Atmo_Modify_var to get the
!                             stratification for a perturbation of
!                             a given variable, what allows to
!                             significantly simplify the source (TdPA)
!                           - Hydrostatic equilibrium is not
!                             mandatory (TdPA)
!                           - Copies of model atmosphere are made
!                             using cAtmo, ensuring no undesired
!                             memory space sharing (TdPA)
!                           - Copies of model atmosphere are fred
!                             with free_Atmo, ensuring that the
!                             memory is correctly released (TdPA)
!                           - Added get_value_corr to significantly
!                             simplify the source code when acquiring
!                             the value of any variable when it is
!                             set as correction (TdPA)
!                           - Fixed typo in verbosity (TdPA)
!
!     04/11/2023:    V3.0.2 - Remove the keyword Hanle_Effect (HL)
!
!     03/15/2023:    V3.0.1 - CheckPerturb now accounts for the
!                             additional limits possible in the
!                             input (TdPA)
!                           - Atmo_Modify does not need the Flgsg
!                             argument (TdPA)
!                           - The Blos variables are in the same
!                             structure than the polar ones (TdPA)
!                           - Removed some commented lines (TdPA)
!                           - Removed some superflous lines (TdPA)
!
!     03/08/2023:    V3.0.0 - First working version (TdPA)
!
!     02/27/2023:    V0.0.0 - Started from 12/05/2020
!                             TIC@rf_mod.f90 revision (TdPA)
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
!    RF_Thermo:
!      Compute the response function of the intensity to the
!    perturbation in a thermal parameter node
!
!    RF_Mag:
!      Compute the response function of the polarization to the
!    perturbation in a magnetic parameter node
!
!    RF_All:
!      Compute the response function of the polarization to the
!    perturbation in a node
!
!    RF_fInt:
!      Compute the response function of the intensity to a
!    perturbation of the diffuse light factor
!
!    RF_fPol:
!      Compute the response function of the polarization to a
!    perturbation of the diffuse light factor
!
!    CheckPerturb:
!      Check that a perturbation is suitable, both up and down,
!    keeping the variable within the valid limits
!
!    Nodes_Perturb:
!      Perturbe a given node for a given variable in a given sign
!    direction
!
!    get_value_corr:
!      Get the value of a given variable from the model atmosphere,
!    for nodes set to specify corrections
!
!    Atmo_Modify_var:
!      Get the stratification for a given variable aster node
!    perturbation
!
!    Atmo_Modify:
!      Modify the atmospheric model to account for the perturbation
!    to one node
!
!#####################################################################
!#####################################################################
!#####################################################################

      use aborted_mod
      use free_mod
      use hanlert_mod
      use hydrostatic_mod
      use initmodel_mod
      use model_mod
      use ratmo_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the response function of the intensity to the
      !! perturbation in a thermal parameter node\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!              Indx(integer): Index of the node in the
      !!                             Jacobian\n
      !!              RF(double(:)): Response function\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!     SolF(Solution_F_class): Class with the full RT
      !!                             solution\n
      !!         Input(Input_class): Structure with settings data
      subroutine RF_Thermo(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                           fudge,kurucz,MPID,Atmo,Bfield,Inf_Nodes, &
                           Indx,RF,Sol,SolF,Input)

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
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      type(Input_class), intent(inout):: Input
      integer, intent(in):: Indx
      double precision, dimension(:), intent(out):: RF

      ! Local
      type(Atmo_class):: Tmp_Atmo

      logical, dimension(2):: Flag

      integer:: Indx_Nodes, Indx_Para

      double precision:: factor

      double precision:: tmp_H, tmp_value, Pg, daux
      double precision, dimension(:), allocatable:: Tmp_Var
      double precision, dimension(:), allocatable:: I_up, I_down
      double precision, dimension(:,:), allocatable:: Stokes_out


      ! Initialize
      factor = 1d0

      ! Get index of parameter and node
      Indx_Para = Inf_Nodes%Inf_Inv(1,indx)
      Indx_Nodes = Inf_Nodes%Inf_Inv(2,indx)

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,i5,3x,A,i5)') &
          ' - Get RF for model parameter = ', Indx_Para, &
          'and node = ', Indx_Nodes
        call verboseI(3)

      end if ! Master

      !
      ! If diffuse light
      !
      if (Indx_Para.eq.Inf_Nodes%index_f) then

        ! RF for diffuse light, only intensity
        call RF_fInt(Inf_Nodes,RF,Sol,Input%centered)

        ! We are done here
        return

      end if

      ! If hydrostatic equilibrium
      if (Inf_Nodes%hydroeq) then

        ! If inverting gas pressure
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg)) then

          ! Get value from node
          Pg = Inf_Nodes%Node(Inf_Nodes%index_Pg)%Var(1)

        ! If not inverting
        else

          ! Get from input
          Pg = Inf_Nodes%Pg_Bound

        end if ! Inverting gas pressure

      ! Otherwise, no value
      else

        Pg = 0d0

      end if

      ! Allocate space for Stokes parameters and temporal
      ! node values
      allocate(I_up(size(Sol%omega_input)))
      if (Input%centered) allocate(I_down(size(Sol%omega_input)))
      allocate(Stokes_out(0:3,size(Sol%omega_input)))
      allocate(Tmp_Var(Inf_Nodes%Num_Nodes(Indx_Para)))

      ! Initialize with current values to recover later
      Tmp_Var = Inf_Nodes%Node(Indx_Para)%Var
      Stokes_out = Sol%Stokes_out

      ! Copy atmosphere
      call cAtmo(Atmo,Tmp_Atmo)

      ! Save node height
      Tmp_H = Inf_Nodes%Node(Indx_Para)%H(Indx_Nodes)

      ! If by value
      if (Inf_Nodes%Node_Type(Indx_Para).le.3) then

        ! Get value from model
        tmp_value = Inf_Nodes%Node(Indx_Para)%Var(Indx_Nodes)

      ! If by correction
      else

        ! Get value from model
        call get_value_corr(Inf_Nodes,Atmo,Bfield, &
                            Indx_Para,Indx_Nodes,tmp_value)

      end if ! Type of node

      ! Correct factor?
      if (Inf_Nodes%min_rel_Pert(Indx_Para).gt.0d0) then

        ! Get absolute
        daux = abs(tmp_value)

        ! Only if larger than 0
        if (daux.gt.0d0) then

          daux = 1d0/daux

          do while (Inf_Nodes%Perturb(Indx_Para)*factor*daux.lt. &
                    Inf_Nodes%min_rel_Pert(Indx_Para))

            factor = factor*2d0

          end do

        end if ! Value larger than zero
      end if ! Correct perturbation

      ! Check if the perturbation is fine
      call CheckPerturb(tmp_H, tmp_value, &
                        Inf_Nodes%Perturb(Indx_Para)*factor, &
                        Inf_Nodes%Node(Indx_Para), Flag)
      if (laborted) return

      ! If not doing centered or if one of the two perturbations is
      ! not valid
      if ((.not.Input%Centered).or. &
          (.not.(Flag(1).and.Flag(2)))) then

        ! Perturb the node
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, Flag(1))

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Bfield,Tmp_Atmo,Atom,Atomb,Mol,Input, &
                         fudge,Inf_Nodes,Indx_Para,Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '-----------------------------------------'
          call verbose
          umsg = '| Synthesis to compute non-centered RF |'
          call verbose
          umsg = '-----------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes in intensity "up"
        I_up(:) = Sol%Stokes_out(0,:)

        ! Get response function
        RF = (I_up - Stokes_out(0,:))* &
             Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

        ! If valid down, and not up
        if (.not.Flag(1)) then

          ! Negative numerator
          RF = -1d0*RF

        end if ! Not valid up

      ! If centered and both perturbations are valid
      else

        ! Perturb the node up
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .True.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Bfield,Tmp_Atmo,Atom,Atomb,Mol, &
                         Input,fudge,Inf_Nodes,Indx_Para, &
                         Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '------------------------------------------'
          call verbose
          umsg = '| Synthesis "up" to compute centered RF |'
          call verbose
          umsg = '------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes in intensity "up"
        I_up(:) = Sol%Stokes_out(0,:)

        ! Perturb the node down
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .False.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Bfield,Tmp_Atmo,Atom,Atomb,Mol,Input,&
                         fudge,Inf_Nodes,Indx_Para,Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '--------------------------------------------'
          call verbose
          umsg = '| Synthesis "down" to compute centered RF |'
          call verbose
          umsg = '--------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes in intensity "down"
        I_down(:) = Sol%Stokes_out(0,:)

        ! Get response
        RF = (I_up - I_down)* &
             0.5d0*Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

      end if ! Centered and valid perturbations in both directions

      ! Recover original output
      Sol%Stokes_out = Stokes_out

      ! Free memory
      deallocate(I_up, Stokes_out, Tmp_Var)
      if (Input%centered) deallocate(I_down)
      call free_Atmo(Tmp_Atmo,.True.)

      return

      end subroutine RF_Thermo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the response function of the polarization to the
      !! perturbation in a magnetic parameter node\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!              Indx(integer): Index of the node in the
      !!                             Jacobian\n
      !!              RF(double(:)): Response function\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!     SolF(Solution_F_class): Class with the full RT
      !!                             solution\n
      !!         Input(Input_class): Structure with settings data
      subroutine RF_Mag(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Bfield,Inf_Nodes,Indx,RF, &
                        Sol,SolF,Input)

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
      type(Nodes_class), intent(inout)::Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      type(Input_class), intent(inout):: Input
      integer, intent(in):: Indx
      double precision, dimension(:,:):: RF

      ! Local
      type(Bfield_class):: Tmp_Bfield

      logical, dimension(2):: Flag

      integer:: Indx_Nodes, Indx_Para

      double precision:: factor

      double precision:: tmp_H, tmp_value, Pg, daux
      double precision, dimension(:), allocatable:: Tmp_Var
      double precision, dimension(:,:), allocatable:: Stokes_up
      double precision, dimension(:,:), allocatable:: Stokes_down
      double precision, dimension(:,:), allocatable:: Stokes_out


      ! Initialize
      factor = 1d0

      ! Get index of parameter and node
      Indx_Para = Inf_Nodes%Inf_Inv(1,indx)
      Indx_Nodes = Inf_Nodes%Inf_Inv(2,indx)

      ! Master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,i5,3x,A,i5)') &
          ' - Get RF for model parameter = ', Indx_Para, &
          'and node = ', Indx_Nodes
        call verboseI(3)

      end if ! Master

      !
      ! If diffuse light
      !
      if (Indx_Para.eq.Inf_Nodes%index_f) then

        ! RF for diffuse light, only intensity
        call RF_fPol(Inf_Nodes,RF,Sol,Input%centered)

        ! We are done here
        return

      end if

      ! Allocate space for Stokes parameters and temporal
      ! node values
      allocate(Stokes_up(0:3,size(Sol%omega_input)))
      allocate(Stokes_out(0:3,size(Sol%omega_input)))
      if (Input%centered) &
        allocate(Stokes_down(0:3,size(Sol%omega_input)))
      allocate(Tmp_Var(Inf_Nodes%Num_Nodes(Indx_Para)))

      ! Store current values to recover later
      Tmp_Var = Inf_Nodes%Node(Indx_Para)%Var
      Tmp_Bfield = Bfield
      Stokes_out = Sol%Stokes_out

      ! Save node height
      Tmp_H = Inf_Nodes%Node(Indx_Para)%H(Indx_Nodes)

      ! If by value
      if (Inf_Nodes%Node_Type(Indx_Para).le.3) then

        ! Get value from model
        tmp_value = Inf_Nodes%Node(Indx_Para)%Var(Indx_Nodes)

      ! By correction
      else

        ! Get value from model
        call get_value_corr(Inf_Nodes,Atmo,Bfield, &
                            Indx_Para,Indx_Nodes,tmp_value)

      end if ! Value or correction


      ! Correct factor?
      if (Inf_Nodes%min_rel_Pert(Indx_Para).gt.0d0) then

        ! Get absolute
        daux = abs(tmp_value)

        ! Only if larger than 0
        if (daux.gt.0d0) then

          daux = 1d0/daux

          do while (Inf_Nodes%Perturb(Indx_Para)*factor*daux.lt. &
                    Inf_Nodes%min_rel_Pert(Indx_Para))

            factor = factor*2d0

          end do

        end if
      end if ! Correct perturbation

      ! Check if the perturbation is fine
      call CheckPerturb(tmp_H, tmp_value, &
                        Inf_Nodes%Perturb(Indx_Para)*factor, &
                        Inf_Nodes%Node(Indx_Para), Flag)
      if (laborted) return

      ! If not doing centered or if one of the two perturbations is
      ! not valid
      if ((.not.Input%Centered).or. &
          (.not.(Flag(1).and.Flag(2)))) then

        ! Perturb the node
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, Flag(1))

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Atmo,Atom,Atomb,Mol,Input, &
                         fudge,Inf_Nodes,Indx_Para,Tmp_Var,0d0)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '-----------------------------------------'
          call verbose
          umsg = '| Synthesis to compute non-centered RF |'
          call verbose
          umsg = '-----------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Tmp_Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes "up"
        Stokes_up = Sol%Stokes_out

        ! Get response function
        RF = (Stokes_up - Stokes_out)* &
             Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

        ! If valid down, and not up
        if (.not.Flag(1)) then

          ! Negative numerator
          RF = -1d0*RF

        end if ! Not valid up

      ! If centered and both perturbations are valid
      else

        ! Perturb the node up
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .True.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Atmo,Atom,Atomb,Mol,Input, &
                         fudge,Inf_Nodes,Indx_Para,Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '------------------------------------------'
          call verbose
          umsg = '| Synthesis "up" to compute centered RF |'
          call verbose
          umsg = '------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Tmp_Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes "up"
        Stokes_up = Sol%Stokes_out

        ! Perturb the node down
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .False.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Atmo,Atom,Atomb,Mol,Input, &
                         fudge,Inf_Nodes,Indx_Para,Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '--------------------------------------------'
          call verbose
          umsg = '| Synthesis "down" to compute centered RF |'
          call verbose
          umsg = '--------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Tmp_Bfield,Input, &
                        Sol,SolF,.True.)
        if (laborted) return


        ! Save Stokes "down"
        Stokes_down = Sol%Stokes_out

        ! Get response function
        RF = (Stokes_up - Stokes_down)* &
             0.5d0*Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

      end if ! Centered derivative

      ! Recover solution
      Sol%Stokes_out = Stokes_out

      ! Free memory
      if (allocated(Stokes_up)) deallocate(Stokes_up)
      if (allocated(Stokes_down)) deallocate(Stokes_down)
      if (allocated(Stokes_out)) deallocate(Stokes_out)
      if (allocated(Tmp_Var)) deallocate(Tmp_Var)

      return

      end subroutine RF_Mag

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the response function of the polarization to the
      !! perturbation in a node\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with the geometry data
      !!                             for the intensity problem\n
      !!       Geom(Geometry_class): Structure with the geometry
      !!                             data\n
      !!         Flgsg(Fctsg_class): Structure with factorials and
      !!                             signs\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!            MPID(MPI_class): Structure with MPI data
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!              Indx(integer): Index of the node in the
      !!                             Jacobian\n
      !!              RF(double(:)): Response function\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!     SolF(Solution_F_class): Class with the full RT
      !!                             solution\n
      !!         Input(Input_class): Structure with settings data
      subroutine RF_All(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec,fudge, &
                        kurucz,MPID,Atmo,Bfield,Inf_Nodes,Indx,RF, &
                        Sol,SolF,Input)

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
      type(Nodes_class), intent(inout)::Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      type(Solution_F_class), intent(inout):: SolF
      type(Input_class), intent(inout):: Input
      integer, intent(in):: Indx
      double precision, dimension(:,:):: RF

      ! Local
      type(Bfield_class):: Tmp_Bfield
      type(Atmo_class):: Tmp_Atmo

      logical, dimension(2):: Flag

      integer:: Indx_Nodes, Indx_Para

      double precision:: factor

      double precision:: tmp_H, tmp_value, Pg, daux
      double precision, dimension(:), allocatable:: Tmp_Var
      double precision, dimension(:,:), allocatable:: Stokes_up
      double precision, dimension(:,:), allocatable:: Stokes_down
      double precision, dimension(:,:), allocatable:: Stokes_out


      ! Initialize
      factor = 1d0

      ! Get index of parameter and node
      Indx_Para = Inf_Nodes%Inf_Inv(1,indx)
      Indx_Nodes = Inf_Nodes%Inf_Inv(2,indx)

      ! If master
      if (pid.eq.0) then

        ! Verbose
        write(umsg,'(A,i5,3x,A,i5)') &
          ' - Get RF for model parameter = ', Indx_Para, &
          'and node = ', Indx_Nodes
        call verboseI(3)

      end if

      !
      ! If diffuse light
      !
      if (Indx_Para.eq.Inf_Nodes%index_f) then

        ! RF for diffuse light, only intensity
        call RF_fPol(Inf_Nodes,RF,Sol,Input%centered)

        ! We are done here
        return

      end if

      ! If hydrostatic equilibrium
      if (Inf_Nodes%hydroeq) then

        ! If inverting gas pressure
        if (Inf_Nodes%Nodes_Flags(Inf_Nodes%index_Pg)) then

          ! Get value from node
          Pg = Inf_Nodes%Node(Inf_Nodes%index_Pg)%Var(1)

        ! If not inverting
        else

          ! Get from input
          Pg = Inf_Nodes%Pg_Bound

        end if ! Inverting gas pressure

      ! Otherwise, no value
      else

        Pg = 0d0

      end if

      allocate(Stokes_up(0:3,size(Sol%omega_input)))
      allocate(Stokes_out(0:3,size(Sol%omega_input)))
      allocate(Stokes_down(0:3,size(Sol%omega_input)))
      allocate(Tmp_Var(Inf_Nodes%Num_Nodes(Indx_Para)))

      ! Store current values to recover later
      Tmp_Var = Inf_Nodes%Node(Indx_Para)%Var
      Tmp_Bfield = Bfield
      Stokes_out = Sol%Stokes_out

      ! Copy atmosphere
      call cAtmo(Atmo,Tmp_Atmo)

      ! Save node height
      Tmp_H = Inf_Nodes%Node(Indx_Para)%H(Indx_Nodes)

      ! If by value
      if (Inf_Nodes%Node_Type(Indx_Para).le.3) then

        ! Get value from model
        tmp_value = Inf_Nodes%Node(Indx_Para)%Var(Indx_Nodes)

      ! If by correction
      else

        ! Get value from model
        call get_value_corr(Inf_Nodes,Atmo,Bfield, &
                            Indx_Para,Indx_Nodes,tmp_value)
      end if ! Type of node

      ! Correct factor?
      if (Inf_Nodes%min_rel_Pert(Indx_Para).gt.0d0) then

        ! Get absolute
        daux = abs(tmp_value)

        ! Only if larger than 0
        if (daux.gt.0d0) then

          daux = 1d0/daux

          do while (Inf_Nodes%Perturb(Indx_Para)*factor*daux.lt. &
                    Inf_Nodes%min_rel_Pert(Indx_Para))

            factor = factor*2d0

          end do

        end if
      end if ! Correct perturbation

      ! If inverting asymmetry, signal input
      if (Inf_Nodes%Num_Asymmetry.gt.0) Input%nasym = 1

      ! Check if perturbation is fine
      call CheckPerturb(tmp_H, tmp_value, &
                        Inf_Nodes%Perturb(Indx_Para)*factor, &
                        Inf_Nodes%Node(Indx_Para), Flag)
      if (laborted) return

      ! If not doing centered or if one of the two perturbations is
      ! not valid
      if ((.not.Input%Centered).or. &
          (.not.(Flag(1).and.Flag(2)))) then

        ! Perturb the node
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, Flag(1))

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Tmp_Atmo,Atom,Atomb,Mol, &
                         Input,fudge,Inf_Nodes,Indx_Para, &
                         Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '-----------------------------------------'
          call verbose
          umsg = '| Synthesis to compute non-centered RF |'
          call verbose
          umsg = '-----------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Tmp_Bfield, &
                        Input,Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes "up"
        Stokes_up = Sol%Stokes_out

        ! Get response function
        RF = (Stokes_up - Stokes_out)* &
             Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

        ! If valid down, and not up
        if (.not.Flag(1)) then

          ! Negative numerator
          RF = -1d0*RF

        end if ! Not valid up

      ! If centered and both perturbations are valid
      else

        ! Perturb the node up
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .True.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Tmp_Atmo,Atom,Atomb,Mol, &
                         Input,fudge,Inf_Nodes,Indx_Para, &
                         Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '------------------------------------------'
          call verbose
          umsg = '| Synthesis "up" to compute centered RF |'
          call verbose
          umsg = '------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Tmp_Bfield, &
                        Input,Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes "up"
        Stokes_up = Sol%Stokes_out

        ! Perturb the node down
        call Nodes_Perturb(Indx_Para, Indx_Nodes, Inf_Nodes, &
                           Tmp_Var, factor, .False.)

        ! Modify the model atmosphere accordingly
        call Atmo_Modify(Tmp_Bfield,Tmp_Atmo,Atom,Atomb,Mol, &
                         Input,fudge,Inf_Nodes,Indx_Para, &
                         Tmp_Var,Pg)

        ! If Master
        if (gpid.eq.0) then

          ! Verbose
          umsg = '--------------------------------------------'
          call verbose
          umsg = '| Synthesis "down" to compute centered RF |'
          call verbose
          umsg = '--------------------------------------------'
          call verbose

        end if ! Master

        ! Call the synthesis
        call HanleRTTIC(Atom,Atomb,Mol,GeomI,Geom,Flgsg,Frec, &
                        fudge,kurucz,MPID,Tmp_Atmo,Tmp_Bfield, &
                        Input,Sol,SolF,.True.)
        if (laborted) return

        ! Save Stokes "down"
        Stokes_down = Sol%Stokes_out

        ! Get response function
        RF = (Stokes_up - Stokes_down)* &
             0.5d0*Inf_Nodes%Scal(Indx_Para)/ &
             Inf_Nodes%Perturb(Indx_Para)/ &
             factor

      end if ! Centered derivative

      ! Recover solution
      Sol%Stokes_out = Stokes_out

      ! Free memory
      if (allocated(Stokes_up)) deallocate(Stokes_up)
      if (allocated(Stokes_down)) deallocate(Stokes_down)
      if (allocated(Stokes_out)) deallocate(Stokes_out)
      if (allocated(Tmp_Var)) deallocate(Tmp_Var)
      call free_Atmo(Tmp_Atmo,.True.)

      return

      end subroutine RF_All

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the response function of the intensity to the
      !! perturbation in the diffuse light factor\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!              RF(double(:)): Response function\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!          centered(logical): If centered derivative
      subroutine RF_fInt(Inf_Nodes,RF,Sol,centered)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: centered
      double precision, dimension(:), intent(out):: RF

      ! Local
      logical, dimension(2):: Flag

      integer:: indx,ir

      double precision:: tmp_H, tmp_value
      double precision, dimension(:), allocatable:: Tmp_Var
      double precision, dimension(:), allocatable:: I_up, I_down
      double precision, dimension(:), allocatable:: I_clean
      double precision, dimension(:), allocatable:: I_diff


      ! Variable index
      indx = Inf_Nodes%index_f

      ! Allocate space for Stokes parameters and temporal
      ! node values
      allocate(I_up(size(Sol%omega_input)))
      if (centered) allocate(I_down(size(Sol%omega_input)))
      allocate(I_clean(size(Sol%omega_input)))
      allocate(I_diff(size(Sol%omega_input)))
      allocate(Tmp_Var(1))

      ! Get scaled diffuse light
      do ir=1,Sol%Num_Range

        ! Scale to the range
        I_diff(Sol%Range(ir,1):Sol%Range(ir,2)) = &
               Sol%Stokes_diff(0,Sol%Range(ir,1):Sol%Range(ir,2))/ &
               Sol%Scal_Stokes(ir)

      end do

      ! Initialize with current values to recover later
      Tmp_Var = Inf_Nodes%Node(indx)%Var(1)
      Tmp_H = Inf_Nodes%Node(indx)%H(1)
      tmp_value = Tmp_var(1)

      ! Get clean profile
      I_clean = Sol%Stokes_out(0,:)
      I_clean = (I_clean - Tmp_var(1)*I_diff)/(1d0 - Tmp_var(1))

      ! Check if the perturbation is fine
      call CheckPerturb(tmp_H, tmp_value, &
                        Inf_Nodes%Perturb(indx), &
                        Inf_Nodes%Node(indx), Flag)
      if (laborted) return

      ! If not doing centered or if one of the two perturbations is
      ! not valid
      if ((.not.centered).or.(.not.(Flag(1).and.Flag(2)))) then

        ! Perturb the node
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, Flag(1))

        ! Save Stokes in intensity "up"
        I_up = I_clean*(1d0 - Tmp_var(1)) + Tmp_var(1)*I_diff

        ! Get response function
        RF = (I_up - Sol%Stokes_out(0,:))* &
             Inf_Nodes%Scal(indx)/ &
             Inf_Nodes%Perturb(indx)

        ! If valid down, and not up
        if (.not.Flag(1)) then

          ! Negative numerator
          RF = -1d0*RF

        end if ! Not valid up

      ! If centered and both perturbations are valid
      else

        ! Perturb the node up
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, .True.)

        ! Save Stokes in intensity "up"
        I_up = I_clean*(1d0 - Tmp_var(1)) + Tmp_var(1)*I_diff

        ! Perturb the node down
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, .False.)

        ! Save Stokes in intensity "down"
        I_down = I_clean*(1d0 - Tmp_var(1)) + Tmp_var(1)*I_diff

        ! Get response
        RF = (I_up - I_down)* &
             0.5d0*Inf_Nodes%Scal(indx)/ &
             Inf_Nodes%Perturb(indx)

      end if ! Centered and valid perturbations in both directions

      ! Free memory
      deallocate(I_up,I_diff,Tmp_Var)
      if (centered) deallocate(I_down)

      return

      end subroutine RF_fint

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the response function of the Stokes parameters to the
      !! perturbation in the diffuse light factor\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!            RF(double(:,:)): Response function\n
      !!        Sol(Solution_class): Class with the data of the RT
      !!                             solution\n
      !!          centered(logical): If centered derivative
      subroutine RF_fPol(Inf_Nodes,RF,Sol,centered)

      ! IO
      type(Nodes_class), intent(inout):: Inf_Nodes
      type(Solution_class), intent(inout):: Sol
      logical, intent(in):: centered
      double precision, dimension(:,:), intent(out):: RF

      ! Local
      logical, dimension(2):: Flag

      integer:: indx,ir

      double precision:: tmp_H, tmp_value
      double precision, dimension(:), allocatable:: Tmp_Var
      double precision, dimension(:,:), allocatable:: Stokes_up
      double precision, dimension(:,:), allocatable:: Stokes_down
      double precision, dimension(:,:), allocatable:: Stokes_clean
      double precision, dimension(:,:), allocatable:: Stokes_diff


      ! Variable index
      indx = Inf_Nodes%index_f

      ! Allocate space for Stokes parameters and temporal
      ! node values
      allocate(Stokes_up(0:3,size(Sol%omega_input)))
      if (centered) allocate(Stokes_down(0:3,size(Sol%omega_input)))
      allocate(Stokes_clean(0:3,size(Sol%omega_input)))
      allocate(Stokes_diff(0:3,size(Sol%omega_input)))
      allocate(Tmp_Var(1))

      ! Store current values to recover later
      Tmp_Var = Inf_Nodes%Node(indx)%Var(1)
      Tmp_H = Inf_Nodes%Node(indx)%H(1)
      tmp_value = Tmp_var(1)

      ! Get scaled diffuse light
      do ir=1,Sol%Num_Range

        ! Scale to the range
        Stokes_diff(:,Sol%Range(ir,1):Sol%Range(ir,2)) = &
               Sol%Stokes_diff(:,Sol%Range(ir,1):Sol%Range(ir,2))/ &
               Sol%Scal_Stokes(ir)

      end do

      !
      ! Get clean value
      !

      ! If fractional
      if (Sol%Fractional) then

        ! For each range in wavelengths, descale intensity
        do ir=1,Sol%Num_Range

          ! De-scale to the range just the intensity
          Stokes_clean(0,Sol%Range(ir,1):Sol%Range(ir,2)) = &
                  Sol%Stokes_out(0,Sol%Range(ir,1):Sol%Range(ir,2))* &
                  Sol%Scal_Stokes(ir)

        end do

        ! De-fractionize
        Stokes_clean(1,:) = 1d-2*Stokes_clean(0,:)
        Stokes_clean(2,:) = 1d-2*Stokes_clean(0,:)
        Stokes_clean(3,:) = 1d-2*Stokes_clean(0,:)

        ! For each range in wavelengths, scale
        do ir=1,Sol%Num_Range

          ! Scale to the range
          Stokes_clean(:,Sol%Range(ir,1):Sol%Range(ir,2)) = &
                  Stokes_clean(:,Sol%Range(ir,1):Sol%Range(ir,2))/ &
                  Sol%Scal_Stokes(ir)

        end do

      end if

      ! Get clean profile
      Stokes_clean = (Stokes_clean - Tmp_var(1)*Stokes_diff)/ &
                     (1d0 - Tmp_var(1))

      ! Check if the perturbation is fine
      call CheckPerturb(tmp_H, tmp_value, &
                        Inf_Nodes%Perturb(indx), &
                        Inf_Nodes%Node(indx), Flag)
      if (laborted) return

      ! If not doing centered or if one of the two perturbations is
      ! not valid
      if ((.not.centered).or.(.not.(Flag(1).and.Flag(2)))) then

        ! Perturb the node
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, Flag(1))

        ! Save Stokes "up"
        Stokes_up = Stokes_clean*(1d0 - Tmp_var(1)) + &
                    Tmp_var(1)*Stokes_diff

        ! Fractional?
        if (Sol%Fractional) then

          Stokes_up(1,:) = 1d2*Stokes_up(1,:)/Stokes_up(0,:)
          Stokes_up(2,:) = 1d2*Stokes_up(2,:)/Stokes_up(0,:)
          Stokes_up(3,:) = 1d2*Stokes_up(3,:)/Stokes_up(0,:)

        end if

        ! Get response function
        RF = (Stokes_up - Sol%Stokes_out)* &
             Inf_Nodes%Scal(indx)/ &
             Inf_Nodes%Perturb(indx)

        ! If valid down, and not up
        if (.not.Flag(1)) then

          ! Negative numerator
          RF = -1d0*RF

        end if ! Not valid up

      ! If centered and both perturbations are valid
      else

        ! Perturb the node up
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, .True.)

        ! Save Stokes "up"
        Stokes_up = Stokes_clean*(1d0 - Tmp_var(1)) + &
                    Tmp_var(1)*Stokes_diff

        ! Fractional?
        if (Sol%Fractional) then

          Stokes_up(1,:) = 1d2*Stokes_up(1,:)/Stokes_up(0,:)
          Stokes_up(2,:) = 1d2*Stokes_up(2,:)/Stokes_up(0,:)
          Stokes_up(3,:) = 1d2*Stokes_up(3,:)/Stokes_up(0,:)

        end if

        ! Perturb the node down
        call Nodes_Perturb(indx, 1, Inf_Nodes, Tmp_Var, 1d0, .False.)

        ! Save Stokes "down"
        Stokes_down = Stokes_clean*(1d0 - Tmp_var(1)) + &
                      Tmp_var(1)*Stokes_diff

        ! Fractional?
        if (Sol%Fractional) then

          Stokes_down(1,:) = 1d2*Stokes_down(1,:)/Stokes_down(0,:)
          Stokes_down(2,:) = 1d2*Stokes_down(2,:)/Stokes_down(0,:)
          Stokes_down(3,:) = 1d2*Stokes_down(3,:)/Stokes_down(0,:)

        end if

        ! Get response function
        RF = (Stokes_up - Stokes_down)* &
             0.5d0*Inf_Nodes%Scal(indx)/ &
             Inf_Nodes%Perturb(indx)

      end if ! Centered derivative

      ! Free memory
      if (allocated(Stokes_up)) deallocate(Stokes_up)
      if (allocated(Stokes_down)) deallocate(Stokes_down)
      if (allocated(Stokes_clean)) deallocate(Stokes_clean)
      if (allocated(Stokes_diff)) deallocate(Stokes_diff)
      if (allocated(Tmp_Var)) deallocate(Tmp_Var)

      return

      end subroutine RF_fPol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Check that a perturbation is suitable, keeping the variable
      !! within valid limits\n
      !!          H(double): Height\n
      !!        Par(double): Node value\n
      !!    Perturb(double): Perturbation\n
      !!   Node(Node_class): Structure with the nodes data\n
      !!   Flag(logical(2)): Validity flags for perturbations up and
      !!                     down
      subroutine CheckPerturb(H,Par,Perturb,Node,Flag)

      ! IO
      type(Node_class), intent(in):: Node
      logical, dimension(2), intent(inout):: Flag
      double precision, intent(in):: H,Par, Perturb

      ! Local
      integer:: j

      ! For each special limit
      do j=1,Node%nebound

        ! If height within this special limit
        if (H.ge.Node%ebound(1,j).and.H.le.Node%ebound(2,j)) then

          ! If perturbation up is below upper limit
          Flag(1) = (Par + Perturb).le.Node%ebound(4,j)

          ! If perturbation down is above upper limit
          Flag(2) = (Par - Perturb).ge.Node%ebound(3,j)

          ! If both invalid
          if ((.not.Flag(1)).and.(.not.Flag(2))) then
            umsg = 'Purturbation is too large'
            urou = 'CheckPerturb'
            call aborted
            return
          end if

          ! Return
          return

        end if ! Height within special limits

      end do ! For each special limit

      ! If the perturbation up is below the upper limit
      Flag(1) = (Par + Perturb).le.Node%Bounds(2)

      ! If the perturbation down is above the lower limit
      Flag(2) = (Par - Perturb).ge.Node%Bounds(1)

      ! If both invalid
      if ((.not.Flag(1)).and.(.not.Flag(2))) then
        umsg = 'Purturbation is too large'
        urou = 'CheckPerturb'
        call aborted
        return
      end if

      return

      end subroutine CheckPerturb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Perturbe a given node for a given variable in a given sign
      !! direction\n
      !!     Indx_Para(integer): Index of the parameter to perturb\n
      !!    Indx_Nodes(integer): Index of the node to perturb\n
      !! Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!     Tmp_var(double(:)): Array with original values for the
      !!                         node variable\n
      !!         factor(double): Factor to multiply the perturbation\n
      !!            Up(logical): If the perturbation is positive
      subroutine Nodes_Perturb(Indx_Para,Indx_Nodes,Inf_Nodes, &
                               Tmp_Var,factor,Up)

      ! IO
      type(Nodes_class), intent(in):: Inf_Nodes
      logical, intent(in):: Up
      integer, intent(in):: Indx_Para, Indx_Nodes
      double precision, intent(in):: factor
      double precision, dimension(:), intent(inout):: Tmp_Var

      ! Local
      double precision:: sig


      ! If up
      if (Up) then
        sig = 1d0
      else
        sig = -1d0
      end if

      ! If by value
      if (Inf_Nodes%Node_Type(Indx_Para).le.3) then

        ! Add the perturbation
        Tmp_Var(Indx_Nodes) = &
                         Inf_Nodes%Node(Indx_Para)%Var(Indx_Nodes) + &
                         Inf_Nodes%Perturb(Indx_Para)*factor*sig

      ! If by correction
      else

        ! Reset perturbations
        Tmp_Var = 0d0

        ! Set for node
        Tmp_Var(Indx_Nodes) = Inf_Nodes%Perturb(Indx_Para)*factor*sig

      end if ! Value or correction

      return

      end subroutine Nodes_Perturb

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get the value of the correction for the appropriate parameter
      !! and node\n
      !!  Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!        Atmo(Atmo_class): Structure with atmospheric data\n
      !!    Bfield(Bfield_class): Structure with the vertical
      !!                          magnetic field data\n
      !!      Indx_Para(integer): Parameter index\n
      !!     Indx_Nodes(integer): Node index\n
      !!       tmp_value(double): Output value
      subroutine get_value_corr(Inf_Nodes,Atmo,Bfield, &
                                Indx_Para,Indx_Nodes,tmp_value)
      ! I/O
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Nodes_class), intent(in)::Inf_Nodes
      integer, intent(in):: Indx_Para, Indx_Nodes
      double precision, intent(out):: tmp_value

      ! Magnetic field strengths or BLOS
      if (Indx_Para.eq.Inf_Nodes%index_B) then

        ! If vertical reference
        if (Inf_Nodes%Btype.eq.0) then

          ! Get value from stratification
          tmp_value = Bfield%Bstrength(Inf_Nodes% &
                                        Node(Indx_Para)% &
                                        Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Bfield%Blos(Inf_Nodes% &
                                      Node(Indx_Para)% &
                                      Tau_Indx(Indx_Nodes))
        end if ! Reference frame

      ! Magnetic field inclination or BPOS
      else if (Indx_Para.eq.Inf_Nodes%index_Bt) then

        ! If vertical reference
        if (Inf_Nodes%Btype.eq.0) then

          ! Get value from model
          tmp_value = Bfield%Btheta(Inf_Nodes% &
                                     Node(Indx_Para)% &
                                     Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Bfield%Bpos(Inf_Nodes% &
                                      Node(Indx_Para)% &
                                      Tau_Indx(Indx_Nodes))

        end if ! Reference frame

      ! Magnetic field azimuth
      else if (Indx_Para.eq.Inf_Nodes%index_Bp) then

        ! If vertical reference
        if (Inf_Nodes%Btype.eq.0) then

          ! Get value from model
          tmp_value = Bfield%Bphi(Inf_Nodes%Node(Indx_Para)% &
                                            Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Bfield%Azimuth(Inf_Nodes% &
                                         Node(Indx_Para)% &
                                         Tau_Indx(Indx_Nodes))

        end if ! Reference frame

      ! Temperature
      else if (Indx_Para.eq.Inf_Nodes%index_T) then

        ! Get value from atmosphere
        tmp_value = Atmo%T(Inf_Nodes%Node(Indx_Para)% &
                                       Tau_Indx(Indx_Nodes))
      ! X velocity
      else if (Indx_Para.eq.Inf_Nodes%index_vx) then

        ! If vertical reference
        if (Inf_Nodes%vtype.eq.0) then

          ! Get value from atmosphere
          tmp_value = Atmo%vx(Inf_Nodes%Node(Indx_Para)% &
                              Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Atmo%vpos(Inf_Nodes%Node(Indx_Para)% &
                                Tau_Indx(Indx_Nodes))

        end if ! Reference frame

      ! Y velocity
      else if (Indx_Para.eq.Inf_Nodes%index_vy) then

        ! If vertical reference
        if (Inf_Nodes%vtype.eq.0) then

          ! Get value from atmosphere
          tmp_value = Atmo%vy(Inf_Nodes%Node(Indx_Para)% &
                              Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Atmo%vphi(Inf_Nodes%Node(Indx_Para)% &
                                Tau_Indx(Indx_Nodes))

        end if ! Reference frame

      ! Vertical velocity
      else if (Indx_Para.eq.Inf_Nodes%index_vz) then

        ! If vertical reference
        if (Inf_Nodes%vtype.eq.0) then

          ! Get value from atmosphere
          tmp_value = Atmo%vz(Inf_Nodes%Node(Indx_Para)% &
                              Tau_Indx(Indx_Nodes))
        ! If LOS/POS
        else

          ! Get value from stratification
          tmp_value = Atmo%vlos(Inf_Nodes%Node(Indx_Para)% &
                                Tau_Indx(Indx_Nodes))

        end if ! Reference frame

      ! Microturbulent velocity
      else if (Indx_Para.eq.Inf_Nodes%index_vm) then

        ! Get value from atmosphere
        tmp_value = Atmo%vmi(Inf_Nodes%Node(Indx_Para)% &
                                       Tau_Indx(Indx_Nodes))

      ! Gas pressure
      else if (Indx_Para.eq.Inf_Nodes%index_Pg) then

        ! Hydrostatic
        if (Inf_Nodes%hydroeq) then

          ! Get value from atmosphere
          tmp_value = Atmo%Pg(1)

        ! No hydrostatic
        else

          ! Get value from atmosphere
          tmp_value = Atmo%Pg(Inf_Nodes%Node(Indx_Para)% &
                                        Tau_Indx(Indx_Nodes))
        end if

      ! Error
      else

        umsg = 'The index of the parameter is not correct'
        urou = 'get_value_corr'
        call aborted
        return

      end if ! Model parameter

      end subroutine get_value_corr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get stratification for a given variable after node
      !! perturbation\n
      !!  Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!          Var(double(:)): Node values\n
      !!        o_var(double(:)): Output stratification\n
      !!            z(double(:)): Optical depth axis\n
      !!             nn(integer): Dimension of optical depth axis\n
      !!           indx(integer): Parameter index
      subroutine Atmo_Modify_var(Inf_Nodes,Var,o_var,z,nn,indx)

      ! I/O
      type(Nodes_class), intent(in):: Inf_Nodes
      integer, intent(in):: nn, indx
      double precision, dimension(:), intent(in):: z, Var
      double precision, dimension(:), intent(inout):: o_var

      ! Local
      double precision, dimension(nn):: i_var

      ! If by value
      if (Inf_Nodes%Node_Type(indx).le.3) then

          ! Interpolate
          call Intpol(Inf_Nodes%Node(indx)%H, Var, &
                      Inf_Nodes%Num_Nodes(indx), z, o_var, nn, &
                      Inf_Nodes%Interpolation, 3)

      ! If by correction
      else

          ! Interpolate
          call Intpol(Inf_Nodes%Node(indx)%H, Var, &
                      Inf_Nodes%Num_Nodes(indx), z, i_var, nn, &
                      Inf_Nodes%Interpolation, 3)

          ! Get new value
          o_var = o_var + i_var

      end if

      end subroutine Atmo_Modify_Var

!#####################################################################
!#####################################################################
!#####################################################################

      !> Modify the atmospheric model to account for the perturbation
      !! to one node\n
      !!       Bfield(Bfield_class): Structure with the vertical
      !!                             magnetic field data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!           Atmo(Atmo_class): Structure with the model\n
      !!           Atom(Atom_class): Structure with the atomic data\n
      !!          Atomb(Atom_class): Structure with the atomic data
      !!                             for background opacities\n
      !!             Mol(Mol_class): Structure with the molecule
      !!                             data\n
      !!         Input(Input_class): Structure with settings data\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!     Inf_Nodes(Nodes_class): Structure with nodes data\n
      !!         Indx_Para(integer): Index of the parameter to
      !!                             perturb\n
      !!         Tmp_var(double(:)): Array with original values for
      !!                             the node variable\n
      !!                 Pg(double): Gas pressure boundary value
      subroutine Atmo_Modify(Bfield,Atmo,Atom,Atomb,Mol,Input, &
                             fudge,Inf_Nodes,Indx_Para, &
                             Tmp_Var,Pg)

      ! IO
      type(Bfield_class), intent(inout):: Bfield
      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:):: Atom
      type(Atom_class), dimension(:), allocatable:: Atomb
      type(Mol_class), dimension(:), allocatable:: Mol
      type(Input_class):: Input
      type(fudge_class):: fudge
      type(Nodes_class), intent(inout)::Inf_Nodes
      integer, intent(in):: Indx_Para
      double precision:: Pg
      double precision, dimension(:), intent(in):: Tmp_Var

      ! Local
      double precision:: Pg_new
      double precision, dimension(:), allocatable:: Z


      ! Allocate height array
      allocate(Z(Atmo%nZ))

      ! Get log tau
      Z = log10(Atmo%z)

      ! Initialize gas pressure
      Pg_new = Pg

      ! Magnetic field strength or BLOS
      if (Indx_para.eq.Inf_Nodes%index_B) then

        ! Vertical
        if (Inf_Nodes%Btype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Bstrength, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else if (Inf_Nodes%Btype.eq.1) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Blos, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Magnetic field inclination or BPOS
      else if (Indx_para.eq.Inf_Nodes%index_Bt) then

        ! Vertical
        if (Inf_Nodes%Btype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Btheta, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else if (Inf_Nodes%Btype.eq.1) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Bpos, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Magnetic field azimuth
      else if (Indx_para.eq.Inf_Nodes%index_Bp) then

        ! Vertical
        if (Inf_Nodes%Btype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Bphi, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else if (Inf_Nodes%Btype.eq.1) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Bfield%Azimuth, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Temperature
      else if (Indx_para.eq.Inf_Nodes%index_T) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%T, &
                             Z,Atmo%nz,Indx_Para)

      ! X velocity
      else if (Indx_para.eq.Inf_Nodes%index_vx) then

        ! Vertical
        if (Inf_Nodes%vtype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vx, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vpos, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Vertical velocity
      else if (Indx_para.eq.Inf_Nodes%index_vy) then

        ! Vertical
        if (Inf_Nodes%vtype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vy, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vphi, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Vertical velocity
      else if (Indx_para.eq.Inf_Nodes%index_vz) then

        ! Vertical
        if (Inf_Nodes%vtype.eq.0) then

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vz, &
                               Z,Atmo%nz,Indx_Para)

        ! LOS
        else

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%vlos, &
                               Z,Atmo%nz,Indx_Para)

        end if ! Vertical/LOS

      ! Microturbulent velocity
      else if (Indx_para.eq.Inf_Nodes%index_vm) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%vmi, &
                             Z,Atmo%nz,Indx_Para)

      ! Gas pressure
      else if (Indx_para.eq.Inf_Nodes%index_Pg) then

        ! If hydrostatic
        if (Inf_Nodes%hydroeq) then

          ! Get new value
          Pg_new = Tmp_Var(1)

        ! No Hydrostatic
        else

          ! Interpolate
          call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                               Atmo%Pg, &
                               Z,Atmo%nz,Indx_Para)

        end if

      ! J21R
      else if (Indx_para.eq.Inf_Nodes%index_J21R) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%JKQin(4*Atmo%nz+1:5*Atmo%nz), &
                             Z,Atmo%nz,Indx_Para)

      ! J21I
      else if (Indx_para.eq.Inf_Nodes%index_J21I) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%JKQin(5*Atmo%nz+1:6*Atmo%nz), &
                             Z,Atmo%nz,Indx_Para)

      ! J22R
      else if (Indx_para.eq.Inf_Nodes%index_J22R) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%JKQin(6*Atmo%nz+1:7*Atmo%nz), &
                             Z,Atmo%nz,Indx_Para)

      ! J22I
      else if (Indx_para.eq.Inf_Nodes%index_J22I) then

        ! Interpolate
        call Atmo_Modify_Var(Inf_Nodes,Tmp_Var, &
                             Atmo%JKQin(7*Atmo%nz+1:8*Atmo%nz), &
                             Z,Atmo%nz,Indx_Para)


      ! Fail
      else

        ! Abort
        umsg = 'The index of the parameter is not correct'
        urou = 'Atmo_Modify'
        call aborted
        return

      end if ! Variable

      ! Free memory
      deallocate(Z)

      ! If magnetic field is in LOS and parameter is magnetic
      if (Inf_Nodes%Btype.eq.1.and. &
          (Indx_Para.eq.Inf_Nodes%index_B.or. &
           Indx_Para.eq.Inf_Nodes%index_Bt.or. &
           Indx_Para.eq.Inf_Nodes%index_Bp)) then

        ! Transform into vertical reference frame
        call Bconversion(Atmo%nZ, Inf_Nodes%mu, Inf_Nodes%azimuth, &
                         Bfield%Blos, Bfield%Bpos, Bfield%Azimuth, &
                         Bfield%Bstrength, Bfield%Btheta, &
                         Bfield%Bphi)

      end if ! LOS magnetic field parameter

      ! If velocity is in LOS and parameter is magnetic
      if (Inf_Nodes%vtype.eq.1.and. &
          (Indx_Para.eq.Inf_Nodes%index_vx.or. &
           Indx_Para.eq.Inf_Nodes%index_vy.or. &
           Indx_Para.eq.Inf_Nodes%index_vz)) then

        ! Transform into vertical reference frame
        call vconversion(Atmo%nZ, Inf_Nodes%mu, Inf_Nodes%azimuth, &
                         Atmo%vlos, Atmo%vpos, Atmo%vphi, &
                         Atmo%vx, Atmo%vy, Atmo%vz)

      end if ! LOS magnetic field parameter

      ! If hydrostatic eq.
      if (Inf_Nodes%hydroeq) then

        ! If perturbing temperature or gas pressure, redo
        ! hydrostatic equilibrium
        if ((Indx_Para.eq.Inf_Nodes%index_T).or. &
            (Indx_Para.eq.Inf_Nodes%index_Pg)) then

          ! Get pressures
          call Compute_Pressure_all(Atmo,Atom,Atomb,Mol, &
                                    Input,fudge,Pg_new)

        end if ! Perturbing T of Pg

      end if ! Hydrostatic equilibrium

      return

      end subroutine Atmo_Modify

!#####################################################################
!#####################################################################
!#####################################################################

      end module rf_mod
