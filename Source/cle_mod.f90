      !> CLE RT
      module cle_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     24/11/2022
!  Last version:
!     15/05/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.1 - Generalized declarations of Atom, Atomb,
!                             and Mol to allow for empty arrays for
!                             any of them (TdPA)
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
!  CLE
!    Solve the radiative transfer equation along the LOS
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use getrtcle_mod
      use iosolution_mod
      use parameters_mod , only: TINYF , vacuum
      use ratmo_mod
      use rtstep_mod
      use rtstepi_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Solve the radiative transfer equation along the LOS\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atomb(Atom_class(:)): Structures with atomic data for
      !!                             background atoms\n
      !!          Mol(Mol_class(:)): Structures with molecular data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!             Red(Red_class): Structure with redistribution
      !!                             input frequency data,
      !!                             redistribution function data, and
      !!                             profile or normalization data\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Flgsg(Fctsg_class): Structure with factorials,
      !!                             signs, and J-symbols\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           batmo(double(:)): Atmospheric model data\n
      !!               x(double(:)): LOS axis if cartesian mode\n
      !!                  y(double): Y coordinate in PoS\n
      !!                  z(double): Y coordinate in PoS\n
      !!           dims(integer(:)): Axes dimensions\n
      !!            bion(double(:)): Read ionization fraction data\n
      !! ion_column_ind(integer(:)): Index of column in buffer for
      !!                             ionization data\n
      !!  ion_value_ind(integer(:)): Index of value in value array for
      !!                             ionization data\n
      !!       ion_value(double(:)): Numeric constant ionization
      !!                             fraction values\n
      !!         spect(spect_class): Structure with the input spectra
      !!                             data\n
      !!     chianti(chianti_class): Structure with the CHIANTI data
      subroutine CLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red,Geom, &
                     Flgsg,fudge,kurucz,batmo,x,y,z,dims,bion, &
                     ion_column_ind,ion_value_ind,ion_value,spect, &
                     chianti)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atomb
      type(Mol_class), dimension(:), &
                       allocatable, intent(inout):: Mol
      type(Fctsg_class), intent(inout):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(Geometry_class), intent(inout):: Geom
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(in):: Red
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(in):: MPID
      type(spect_class), intent(inout):: spect
      type(chianti_class), intent(in):: chianti
      integer, dimension(:), intent(in):: ion_value_ind,ion_column_ind
      integer, dimension(:), intent(in):: dims
      double precision, intent(in):: y,z
      double precision, dimension(:), intent(in):: ion_value,bion,x
      double precision, dimension(:), intent(inout):: batmo

      ! Local

      type(Atmo_class):: Atmo

      logical:: skip,indisk

      integer:: if0,if1,m,o,p,ix0,ix,nx,ios

      double precision:: r,xl,dsm,dsp,maxab,zl
      double precision, dimension(:), allocatable, target:: tau
      double precision, dimension(:,:), allocatable, target:: tau1

      ! Pointers
      double precision, dimension(:,:,:), pointer:: data1M,data1O, &
                                                    data1P
      double precision, dimension(:,:), pointer:: p_K0M, p_K1M, &
                                                  p_K2M, &
                                                  p_SM, p_StkM
      double precision, dimension(:,:), pointer:: p_K0O, p_K1O, &
                                                  p_K2O, &
                                                  p_SO, p_StkO
      double precision, dimension(:,:), pointer:: p_K0P, p_SP
      double precision, dimension(:), pointer:: p_etaIM, p_etaIO
      double precision, dimension(:), pointer:: p_tauM, p_x


      !
      ! Initialize pointers
      !
      nullify(data1M,data1O,data1P)
      nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM,p_etaIM,p_tauM)
      nullify(p_K0O,p_K1O,p_K2O,p_SO,p_StkO,p_etaIO)
      nullify(p_K0P,p_SP,p_x)

      !
      ! Initialize indexes in Atmo and get LOS axis
      !
      call rAtmo_cle_init(batmo,Input,x,y,z,Atmo,dims)


      ! write geometry in output
      call write_CLEgeom(Input%folder,Atmo,Input%lim_stk, &
                         Input%out_tau1)

      !
      ! We are going to cheat the code to think that there is
      ! only one ''height'' node, so the internal routines
      ! can be more easily dealt with. We store the real
      ! number of nodes in nx
      !
      nx = nz
      nz = 1
      Rz0 = 1
      Rz1 = 1

      ! Limits in frequency for this CPU
      if0 = MPID%if0(pid)
      if1 = MPID%if1(pid)

      ! Allocate O pointer for RT coefficients
      allocate(data1O(0:3,if0:if1,0:5))

      ! If tau in output
      if (Input%out_tau1) then

        ! Allocate and initialize
        allocate(tau(MPID%nf(pid)))
        tau = 0d0

      ! Not computing tau
      else

        ! Allcoate dummy and initialize it
        allocate(tau(1))
        tau = 0d0

      end if ! Tau is needed

      !
      ! 3D modes
      !

      ! Cartesian or non-cartesian grids
      if (Atmo%mode.eq.0.or.Atmo%mode.eq.2) then

        ! Allocate M pointer for RT coefficients
        allocate(data1M(0:3,if0:if1,0:5))

        ! Allocate tau1 if in output
        if (Input%out_tau1) allocate(tau1(2,MPID%nf(pid)))

        !
        ! Solve RT transfer problem
        !

        ! Initialize initial ix
        ix0 = 1

        ! Initialize skip column
        skip = .False.

        ! Initialize in front of the disk
        indisk = .False.

        ! Check if the LOS intersects with the disk
        if (real(sqrt(Atmo%ypos*Atmo%ypos + &
                      Atmo%zpos*Atmo%zpos)).le. &
            1d0-TINYF) then

          ! If skipping disk
          if (Input%skip_disk) then

            ! Set to true
            skip = .True.

          ! If doing LOS even if the disk is behind
          else

            ! Start from 0 as a negative flag
            ix0 = 0

            ! For every point along the LOS
            do ix=1,Atmo%nz

              ! Get coordinate
              xl = Atmo%z(ix)

              ! If behind, continue searching
              if (xl.lt.0.) cycle

              ! If in front, compute distance from center
              r = sqrt(xl*xl + Atmo%ypos*Atmo%ypos + &
                               Atmo%zpos*Atmo%zpos)

              ! If above surface, you found the point
              if (r.gt.1d0) then
                ix0 = ix
                exit
              end if

            end do ! LOS positions

            ! If there are no points in front of the Sun, skip
            if (ix0.lt.1) then

              ! Set skip flag
              skip = .True.

            ! Found a point in front of the disk (and outside)
            else

              ! Idenfity first point as in front of the disk
              indisk = .True.

            end if ! Found an index in front of the disk
          end if ! Skipping disk
        end if ! LOS cuts disks somewhere


        !
        ! If skipping, just put zero in Stokes and point to it
        !
        if (skip) then

          data1M(:,:,5) = 0d0
          p_StkO => data1M(:,:,5)
          if (Input%out_tau1) tau = 0d0

        !
        ! Not skipping, do RT
        !
        else

          ! If there are more than one point, or is disk
          if (nx.gt.1.or.indisk) then

            !
            ! Back point
            !

            ! Identify index
            o = ix0

            ! Get RT coefficients for o
            call getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red, &
                          Geom,Flgsg,fudge,kurucz,o,if0,if1,batmo, &
                          bion,ion_column_ind,ion_value_ind, &
                          ion_value,spect,chianti,data1M, &
                          indisk,.True.)
            ! Error
            if (laborted) goto 1000

            ! If output tau1
            if (Input%out_tau1) then

              ! Initialize and point to it
              tau1 = 0d0
              p_tauM => tau1(1,:)

            end if

            ! More than one point
            if (nx.gt.1) then

              ! Identidy next point along LOS
              p = ix0 + 1

              ! Get RT coefficients for p
              call getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red, &
                            Geom,Flgsg,fudge,kurucz,p,if0,if1,batmo, &
                            bion,ion_column_ind,ion_value_ind, &
                            ion_value,spect,chianti,data1O, &
                            indisk,.False.)
              ! Error
              if (laborted) goto 1000

            end if ! More than one point
          end if ! More than one point or in disk

          !
          ! Intermediate points
          !
          do ix=ix0+1,nx-1

            ! Allocate P pointers
            allocate(data1P(0:3,if0:if1,0:5))

            ! Identify point in the LOS
            m = ix - 1
            o = ix
            p = ix + 1

            ! Calculate distance to previous point (the z variable
            ! is the x axis, which is the LOS)
            dsm = Atmo%z(o) - Atmo%z(m)
            dsm = dsm*Input%R_star

            ! Caculate quantities of the next point (the z variable
            ! is the x axis, which is the LOS)
            dsp = Atmo%z(p) - Atmo%z(o)
            dsp = dsp*Input%R_star

            ! Get RT coefficients for p
            call getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red, &
                          Geom,Flgsg,fudge,kurucz,p,if0,if1,batmo, &
                          bion,ion_column_ind,ion_value_ind, &
                          ion_value,spect,chianti,data1P, &
                          indisk,.False.)
            ! Error
            if (laborted) goto 1001

            ! Point to the data
            p_K0M  => data1M(:,:,0)
            p_K1M  => data1M(:,:,1)
            p_K2M  => data1M(:,:,2)
            p_SM   => data1M(:,:,4)
            p_StkM => data1M(:,:,5)
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)
            p_K0P  => data1P(:,:,0)
            p_SP   => data1P(:,:,4)

            ! Apply short characteristics BESSER
            call RTStep(o,1,1,MPID%nf(pid), &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.True.)

            ! Compute optical path
            if (Input%out_tau1) then

              ! Point
              p_etaIM => data1M(0,:,0)
              p_etaIO => data1O(0,:,0)

              ! Compute tau
              call RTtauI(dsm,MPID%nf(pid),Atmo%z(m), &
                          Atmo%z(o),p_etaIM,p_etaIO, &
                          p_tauM,tau,tau1)

              ! Shift tau data
              p_tauM => tau

            end if

            ! Shift data (O->M, P->O)
1001        deallocate(data1M)
            data1M => data1O
            data1O => data1P
            nullify(data1P)

            ! If error, leave now
            if (laborted) exit

          end do ! Intermediate heights

          ! Error
          if (laborted) goto 1000

          !
          ! Last height
          !

          ! If there was more than one point
          if (nx.gt.1) then

            ! Identify point in the LOS
            m = nx - 1
            o = nx

            ! Calculate distance to previous point (the z variable
            ! is the x axis, which is the LOS)
            dsm = Atmo%z(o) - Atmo%z(m)
            dsm = dsm*Input%R_star

            ! Point to the data
            p_K0M  => data1M(:,:,0)
            p_K1M  => data1M(:,:,1)
            p_K2M  => data1M(:,:,2)
            p_SM   => data1M(:,:,4)
            p_StkM => data1M(:,:,5)
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)
            ! These are not used, but just in case there was just
            ! one point
            p_K0P  => data1O(:,:,1)
            p_SP   => data1O(:,:,4)

            ! Apply short characteristics BESSER
            call RTStep(o,1,1,MPID%nf(pid), &
                        dsm,dsp,p_K0M,p_K1M,p_K2M, &
                        p_SM,p_K0O,p_K1O,p_K2O, &
                        p_SO,p_K0P,p_SP,p_StkM, &
                        p_StkO,.False.)

            ! Compute optical path
            if (Input%out_tau1) then

              ! Point
              p_etaIM => data1M(0,:,0)
              p_etaIO => data1O(0,:,0)

              ! Compute tau
              call RTtauI(dsm,MPID%nf(pid),Atmo%z(m), &
                          Atmo%z(o),p_etaIM,p_etaIO, &
                          p_tauM,tau,tau1)

              ! Shift tau data
              p_tauM => tau

            end if ! Tau
          end if ! Point to do something
        end if ! Skip LOS or not

        !
        ! Clean local pointers
        !
1000    nullify(p_K0M,p_K1M,p_K2M,p_SM,p_StkM)
        nullify(p_K0O,p_K1O,p_K2O,p_SO)
        nullify(p_K0P,p_SP)
        nullify(data1P)
        deallocate(data1M)
        nullify(data1M,data1P)
        if (Input%out_tau1) then
          nullify(p_etaIM,p_etaIO,p_tauM)
        end if

      !
      ! Slab model
      !
      else if (Atmo%mode.eq.1) then

        ! Height on the plane of the sky
        zl = (Atmo%z(1)+1d0)*sin(Atmo%ypos)

        ! Initialize skip column
        skip = .False.

        ! Initialize in front of the disk
        indisk = .False.

        ! If the LOS intersects with the disk
        if (real(zl).le.1d0-TINYF) then

          ! If skipping disk
          if (Input%skip_disk) then

            ! Set to true
            skip = .True.

          ! If doing LOS even if the disk is behind
          else

            ! Set in disk to true
            indisk = .True.

          end if ! Skipping disk
        end if ! LOS cuts disk


        !
        ! If skipping
        !
        if (skip) then

          ! Just make it zero and point to it
          data1O(:,:,5) = 0d0
          p_StkO => data1M(:,:,5)
          if (Input%out_tau1) tau = 0d0

        !
        ! Not skipping, RT
        !
        else

          ! Get RT coefficients
          call getRTCLE(Atom,Atomb,Mol,Atmo,MPID,Input,Frec,Red, &
                        Geom,Flgsg,fudge,kurucz,1,if0,if1,batmo, &
                        bion,ion_column_ind,ion_value_ind, &
                        ion_value,spect,chianti,data1O, &
                        indisk,.True.)
          ! Error
          if (laborted) goto 2000

          ! Tau (input)
          dsp = Atmo%zpos

          ! If no optical depth, just emissivity
          if (dsp.le.0d0) then

            ! For each stokes parameter
            do m=0,3
              data1O(m,:,5) = data1O(m,:,5) + &
                              data1O(m,:,4)*(data1O(0,:,0) + vacuum)
            end do

            ! Pointer
            p_StkO => data1O(:,:,5)

          ! There is optical depth
          else

            ! Point to the data
            p_K0O  => data1O(:,:,0)
            p_K1O  => data1O(:,:,1)
            p_K2O  => data1O(:,:,2)
            p_SO   => data1O(:,:,4)
            p_StkO => data1O(:,:,5)

            ! The initial maximum is the local
            maxab = maxval(p_K0O(1,:))

            ! If splitting in frequencies
            if (nproc.gt.1) then

              ! Get maximum absorption
              call MPI_ALLREDUCE(MPI_IN_PLACE,maxab,1, &
                                 MPI_DOUBLE_PRECISION, &
                                 MPI_MAX,MPI_COMM_RT,ios)

            end if ! Splitting in frequencies

            ! Reduction factor
            maxab = dsp/maxab

            ! Optical depth
            if (Input%out_tau1) tau = p_K0O(1,:)*maxab

            ! Compute Stokes
            call RTCStep(MPID%nf(pid),maxab,p_K0O, &
                         p_K1O,p_K2O,p_SO,p_StkO)

          end if ! Optical depth
        end if ! Skipping calculation

        !
        ! Clean local pointers
        !
2000    nullify(p_K0O,p_K1O,p_K2O,p_SO)

      end if ! Type of atmospheric model

      ! Error
      if (laborted) goto 3000

      !
      ! Write in output
      !
      call write_CLE(Input%folder,if0,if1,MPID%nf(pid),p_StkO, &
                     tau,Input%lim_stk,Input%out_tau1)


      !
      ! Clean local pointers
      !
3000  nullify(p_StkO)
      deallocate(data1O)
      nullify(data1O)
      deallocate(tau)

      ! Give back the dimensionality
      nz = nx

      return

      end subroutine CLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module cle_mod
