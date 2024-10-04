      !> Computes radiation transfer coefficients for the NLTE problem
      !! of the second kind
      module rtcoeff_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/27/2017
!  Last version:
!     10/04/2024 V3.0.15
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/04/2024:   V3.0.15 - Allocate a dummy p_prof in CLE (TdPA)
!
!     09/23/2024:   V3.0.14 - Added add_cont_cle argument to
!                             RTCoeff_CLE, that allows skipping the
!                             continuum contribution (TdPA)
!
!     04/01/2024:   V3.0.13 - Added optional (through asserts) call to
!                             routine to compute the magnetic first
!                             order radiative transfer coefficients
!                             together (TdPA)
!                           - Added call to routine to compute the non
!                             magnetic first order radiative transfer
!                             coefficients together (TdPA)
!
!     11/16/2023:   V3.0.12 - When doing LOS, AD second order
!                             emissivity should be called with the
!                             direction index always 1 (TdPA)
!
!     11/14/2023:   V3.0.11 - Split the calls to emiss2ord into AA
!                             and AD (TdPA)
!                           - The index for Red%indx should be for
!                             direction index jddir, and not the
!                             jcdir one (TdPA)
!
!     10/16/2023:   V3.0.10 - Made LTElines allocatable to satisfy
!                             memory warnings (TdPA)
!
!     10/04/2023:    V3.0.9 - Bugfix: When not keeping the LTE lines
!                             profiles, point to a dummy (TdPA)
!
!     09/29/2023:    V3.0.8 - Apply Krad limit to continuum (TdPA)
!
!     09/21/2023:    V3.0.7 - Outputs RT coefficients are always
!                             scaled to eta_I (TdPA)
!
!     08/17/2023:    V3.0.6 - Added pointer cleaning when there are
!                             file errors (TdPA)
!
!     08/07/2023:    V3.0.5 - Added the contribution of LTE lines to
!                             the RT coefficients (TdPA)
!
!     11/24/2022:    V3.0.4 - Added RTCoeff_CLE. The second order
!                             emissivity has not been implemented
!                             there yet (TdPA)
!
!     11/10/2022:    V3.0.3 - Added the JKQa argument in RTCoeff and
!                             RTCoeffe, which are passed to
!                             emiss2ord and emiss2ordNB (TdPA)
!
!     10/25/2022:    V3.0.2 - Nullify pointers as when starting each
!                             routine (TdPA)
!                           - Clean pointers at exit (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case the following
!                             changes were needed:
!                              o Atmo%v has changed to Atmo%vx,%vy,
!                                and %vz.
!                              o The dimensionality of the LOS
!                                geometrical tensors has changed.
!                             (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added branching in the magnetic branch
!                             of the b-b RT coefficients when there
!                             is OpenMP (TdPA)
!
!     09/11/2020:    V1.6.3 - Changes to accomodate the new Frec and
!                             Red structures (TdPA)
!
!     06/02/2020:    V1.6.2 - Limited size of directional dimension
!                             of the normalization array in cases
!                             where it is possible (TdPA)
!
!     11/19/2019:    V1.6.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Added stored offsets in frequency for
!                             the use of Voigt files (TdPA)
!
!     11/13/2019:    V1.6.0 - Now admits Voigt profiles stored in a
!                             file (TdPA)
!                           - Photoionization quantities can be in
!                             RAM without storing Voigt profiles for
!                             intensity (TdPA)
!
!     10/03/2019:    V1.5.1 - trano is now indexed in Frec and Red
!                             structures (TdPA)
!
!     05/31/2019:    V1.5.0 - Changed the dimensionality of the
!                             profile variable. Now it runs
!                             sequentially on atoms, transitions and
!                             frequencies to save memory and reduce
!                             the size of data shared through MPI
!                             messages (TdPA)
!                           - Added the routine RTcoeffe for
!                             similarity with the intensity case
!                             and avoid non needed calculations (TdPA)
!
!     05/08/2019:    V1.4.0 - Got rid of the (atomic,transition) pair
!                             of indexes in every radiation tensor and
!                             now they have been compressed in just
!                             one dimension (TdPA)
!                           - Bugfix: There must be a index for
!                             directions to be used by Frec%stype that
!                             is limited to Red%njdir, before I got
!                             out of bounds for some cases (TdPA)
!
!     03/18/2019:    V1.3.1 - Removed searchs by looping (TdPA)
!
!     02/20/2019:    V1.3.0 - Using specific TINY variables (TdPA)
!                           - This routine does not need the
!                             aborted dependence (TdPA)
!
!     01/30/2019:    V1.2.1 - Changed axial by RTaxial (TdPA)
!
!     08/06/2018:    V1.2.0 - Changed how the argument Atom%Norm is
!                             passed in every b-b routine (TdPA)
!                           - Split the call of photoeps depending on
!                             what is stored in RAM (TdPA)
!
!     05/16/2018:    V1.1.2 - Changed the azimuthal dimension of
!                             Stokes from nPh2 to nPh (TdPA)
!                           - Added stype to the calls to emiss2ord
!                             and emiss2ordNB (TdPA)
!
!     10/04/2017:    V1.1.1 - The absorption profile that goes into
!                             the integral is not corrected by stim.
!                             emission (TdPA)
!
!     10/03/2017:    V1.1.0 - Implemented Non magnetic case in b-b
!                             transitions (TdPA)
!                           - Had forgotten to implement it in
!                             RTAbs too
!
!     09/08/2017:    V1.0.6 - Avoid divisions by 0 (TdPA)
!                           - emiss2ord does not need jdir (TdPA)
!
!     08/21/2017:    V1.0.5 - Bugfix: if no PRD, jbdir and jcdir
!                             should not be defined (TdPA)
!
!     06/28/2017:    V1.0.4 - Receives and passes Red (TdPA)
!                           - Passing ia, jdir, jbdir and jcdir to
!                             emiss2 (TdPA)
!                           - Passing the full trano to emiss2 (TdPA)
!                           - The limits of emiss2 now are in other
!                             location in Frec (TdPA)
!
!     06/20/2017:    V1.0.3 - Bugfix: Integral of profile should not
!                             touch boundaries in the loop (TdPA)
!
!     06/19/2017:    V1.0.2 - Changed intent of Frec (TdPA)
!                           - Split one index of the RT variables to
!                             avoid array copies (TdPA)
!                           - Data stored in data1 directly (TdPA)
!
!     06/12/2017:    V1.0.1 - The frequency limits are inputs (TdPA)
!                           - The is no Cont%l anymore (TdPA)
!                           - Passing transition limits into auxiliar
!                             routines (TdPA)
!                           - Stimulated emission is 1st order (TdPA)
!
!     04/27/2017:    V1.0.0 - First version (TdPA)
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
!  RTCoeff
!      This routine calculates the RT coefficients in quadrature
!    directions
!
!  RTCoeffe
!      This routine calculates the RT coefficients in emerging
!    directions
!
!  RTAbs
!    Calculates just intensity absorptivity
!
!  RTCoeff_CLE
!      This routine calculates the RT coefficients for a CLE point
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use math_mod
      use parameters_mod , only : convF , c , c2 , vacuum , TINYB
      use rtcoeffaux_mod
      use rtcoeffiaux_mod , only : absorbILTE , rt1ordILTE
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the radiation transfer coefficients for solver or
      !! solver_serial\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!     JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients\n
      !!        data2(dfloat(:,:)): Line profiles
      subroutine RTCoeff(Frec,Red,Atom,LTElines,Atmo,MPID,Flgsg, &
                         Geom,iz,ith,iph,if0,if1,JKQa,JKQ,JKQC,cdir, &
                         Cont,Bfield,Stokes,data1,data2)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(inout):: Red
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(0:3,if0:if1,0:4):: data1
      double precision, dimension(:,:):: data2
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(in):: JKQC

      ! Local

      integer:: iS,K,iQ,iterml,itermu,ia,jtran,ktran,ltran,icom
      integer:: jdir,jbdir,jcdir,jddir,icdir,ilevell,ilevelu,nodir
      integer:: ifreq,if0l,if1l,if0l2,if1l2,t0,t1,iil,jjl,nf
      integer:: offset,ios,indx,indxf

      double precision:: daux,DwT,pE,absK,Dw,vfac,ct,st,cc,sc,loffset
      double precision, dimension(if0:if1,0:3):: etaA,epsA,rhsA,rhaA
      double precision, dimension(if0:if1):: etmp0,estmp0,es2tmp0
      double precision, dimension(if0:if1):: etmp1,estmp1, &
                                             rstmp1,rtmp1,es2tmp1
      double precision, dimension(if0:if1):: etmp2,estmp2, &
                                             rstmp2,rtmp2,es2tmp2
      double precision, dimension(if0:if1):: etmp3,estmp3, &
                                             rstmp3,rtmp3,es2tmp3
      double precision, dimension(if0:if1):: intgr

      complex(kind=8), dimension(0:3,-2:2,0:2):: TBo,TSo
      complex(kind=8), dimension(:,:), pointer:: prof,p_prof

      type(Frequencyc2_class), pointer:: p_frec
      type(Redc2_class), pointer:: p_red
      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(prof,p_prof,p_frec,p_red,p_Norm,p_LTEprof)


      ! If formal solution direction
      jdir = Geom%i_geom(iph,ith)
      TBo = Geom%TB(:,:,:,iph,ith,iz)
      TSo = Geom%TS(:,:,:,iph,ith)
      icdir = min(jdir,cdir)

      ! If there are dynamics
      if (dyn) then
        nodir = jdir
        ct = Geom%V_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = Geom%v_mux(iph)
        sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)
      else
        nodir = 1
      end if

      if (PRD) then
        jbdir = min(jdir,Frec%ndir)
        jcdir = min(jdir,Red%ndir)
        jddir = min(jdir,Red%njdir)
      end if


      ! If using file, nullify pointer to norm and allocate profile
      if (VPFIL) then

        p_Norm => Atom(1)%Normp(1,1,1)
        allocate(prof(if0:if1,Atom(1)%Mncom))

      ! If not using file
      else

        allocate(p_prof(1,1))

      end if


      !
      ! Initialize to 0 the parts without continuum
      !
      data1(1,:,0) = 0d0

      ! If not axial
      if (.not.RTaxial) then
        data1(2,:,0) = 0d0
        data1(3,:,0) = 0d0
        data1(2,:,1) = 0d0
        data1(3,:,1) = 0d0
        data1(3,:,2) = 0d0
      end if


      !
      ! Continuum contribution
      !

      ! Absorptivity
      data1(0,:,0) = Cont(:,1,icdir)

      ! For each Stokes parameter
      do iS=0,3

        ! if axial skip U and V
        if (RTaxial.and.iS.gt.1) cycle

        ! Reset integral
        intgr = .0D0

        !
        ! Compute the sum over K and Q of TKQ*JKQ(k)
        !

        ! For each K
        do K=0,Krad

          ! For each Q
          do iQ=-K,K

            ! Add contribution to the integral
            intgr = intgr + dble(TSo(iS,iQ,K)*JKQC(iQ,K,if0:if1))

          end do ! Q
        end do ! K

        ! Emissivity by scattering
        data1(iS,:,4) = intgr*Cont(:,2,icdir)

      end do ! Stokes parameters

      ! Add thermal emissivity
      data1(0,:,4) = data1(0,:,4) + cont(:,3,icdir)


      !
      ! Calculate the Doppler shift factor
      !

      ! If there are dynamics
      if (dyn) then

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If not, no need to calculate
      else

        vfac = 1d0

      end if ! dynamics

      !
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If this CPU does not have frequencies in this line, skip
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Photon energy (cgs) and convertion factor
        pE = convF*LTElines(ia)%Dfreq
        absK = 1d21*(2d0*c)*LTElines(ia)%Dfreq**2d0

        ! Add the microt. to Doppler width
        Dw = LTElines(ia)%Dfreq*sqrt(DwT*DwT + &
                                     Atmo%vmi(iz)**2d0)

        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Get Prof
          if (LVPRAM) then
            p_LTEprof => LTElines(ia)%prof(iz,jdir)
          else
            p_LTEprof => LTElines(1)%prof(1,1)
          end if

          call rt1ordLTE(LTElines(ia),TBo,Frec%omega,Flgsg,iz, &
                         if0l,if1l,p_LTEprof,Dw,vfac, &
                         Bfield%Bstrength(iz),pE, &
                         etmp0(if0l:if1l),etmp1(if0l:if1l), &
                         etmp2(if0l:if1l),etmp3(if0l:if1l), &
                         rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                         rtmp3(if0l:if1l), &
                         estmp0(if0l:if1l),estmp1(if0l:if1l), &
                         estmp2(if0l:if1l),estmp3(if0l:if1l), &
                         rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                         rstmp3(if0l:if1l))

          !
          ! Stimulated emission contribution
          !

          ! If there is stimulated emission
          if (stm) then

            ! Correct for stimulated emission
            etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                               estmp0(if0l:if1l)/absK
            etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                               estmp1(if0l:if1l)/absK
            etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                               estmp2(if0l:if1l)/absK
            etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                               estmp3(if0l:if1l)/absK
            rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                               rstmp1(if0l:if1l)/absK
            rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                               rstmp2(if0l:if1l)/absK
            rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                               rstmp3(if0l:if1l)/absK

          endif ! Stimulated emission

          ! Multipli by the population of the atom the contribution of
          ! this atom to the RT coefficients
          data1(0,if0l:if1l,0) = data1(0,if0l:if1l,0) + &
                                 etmp0(if0l:if1l)*LTElines(ia)%n(iz)
          data1(1,if0l:if1l,0) = data1(1,if0l:if1l,0) + &
                                 etmp1(if0l:if1l)*LTElines(ia)%n(iz)
          data1(0,if0l:if1l,4) = data1(0,if0l:if1l,4) + &
                               estmp0(if0l:if1l)*pE*LTElines(ia)%n(iz)
          data1(1,if0l:if1l,4) = data1(1,if0l:if1l,4) + &
                               estmp1(if0l:if1l)*pE*LTElines(ia)%n(iz)

          if (.not.RTaxial) then
            data1(2,if0l:if1l,0) = data1(2,if0l:if1l,0) + &
                                   etmp2(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,0) = data1(3,if0l:if1l,0) + &
                                   etmp3(if0l:if1l)*LTElines(ia)%n(iz)
            data1(2,if0l:if1l,1) = data1(2,if0l:if1l,1) + &
                                   rtmp3(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,1) = data1(3,if0l:if1l,1) - &
                                   rtmp2(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,2) = data1(3,if0l:if1l,2) + &
                                   rtmp1(if0l:if1l)*LTElines(ia)%n(iz)
            data1(2,if0l:if1l,4) = data1(2,if0l:if1l,4) + &
                               estmp2(if0l:if1l)*pE*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,4) = data1(3,if0l:if1l,4) + &
                               estmp3(if0l:if1l)*pE*LTElines(ia)%n(iz)
          end if

        ! No magnetic field
        else

          ! Get Prof
          if (LVPRAM) then
            p_LTEprof => LTElines(ia)%prof(iz,nodir)
          else
            p_LTEprof => LTElines(1)%prof(1,1)
          end if

          ! RT coeffs
          call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_LTEprof,Dw,vfac,pE, &
                          etmp0(if0l:if1l),estmp0(if0l:if1l))

          !
          ! Stimulated emission contribution
          !

          ! If there is stimulated emission
          if (stm) then

            ! Correct for stimulated emission
            etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                               estmp0(if0l:if1l)/absK

          endif ! Stimulated emission

          ! Multipli by the population of the atom the contribution of
          ! this atom to the RT coefficients
          data1(0,if0l:if1l,0) = data1(0,if0l:if1l,0) + &
                                 etmp0(if0l:if1l)*LTElines(ia)%n(iz)
          data1(0,if0l:if1l,4) = data1(0,if0l:if1l,4) + &
                               estmp0(if0l:if1l)*pE*LTElines(ia)%n(iz)

        end if ! Magnetic field

      end do ! LTE lines

      ! Initialize profile coefficient
      iil = 1

      ! For each atom
      do ia=1,nA

        ! Initialize atomic RT coefficients
        epsA = .0D0
        rhsA = .0D0
        etaA = .0D0
        rhaA = .0D0

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))


        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)
            nf = if1l - if0l

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(jdir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tBsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncom
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l,1:Atom(ia)%i_Vind(jtran)%ncom)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,jdir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)
#ifdef RDIPEV
            !
            ! First order RT coefficients
            !
            call rt1ord(Atom(ia),TBo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l), &
                        estmp0(if0l:if1l),estmp1(if0l:if1l), &
                        estmp2(if0l:if1l),estmp3(if0l:if1l), &
                        rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                        rstmp3(if0l:if1l))
#else
#ifdef _OPENMP
            ! If dividing in components
            if (Atom(ia)%omp_comp_1ord(jtran)) then

              !
              ! First order emissivity
              !
              call emiss_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                           jtran,itermu,iterml,iz,if0l,if1l, &
                           p_Norm,Dw,vfac,p_prof, &
                           Atom(ia)%omp_1c(jtran), &
                           estmp0(if0l:if1l),estmp1(if0l:if1l), &
                           estmp2(if0l:if1l),estmp3(if0l:if1l), &
                           rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                           rstmp3(if0l:if1l))

              !
              ! Absorptivity
              !
              call absorb_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            Atom(ia)%omp_1c(jtran), &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l))

            ! Dividing in frequencies
            else
#endif

            !
            ! First order emissivity
            !
            call emiss(Atom(ia),TBo,Frec%omega,Flgsg, &
                       jtran,itermu,iterml,iz,if0l,if1l, &
                       p_Norm,Dw,vfac,p_prof, &
                       estmp0(if0l:if1l),estmp1(if0l:if1l), &
                       estmp2(if0l:if1l),estmp3(if0l:if1l), &
                       rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                       rstmp3(if0l:if1l))

            !
            ! Absorptivity
            !
            call absorb(Atom(ia),TBo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l))
#ifdef _OPENMP
            end if ! Dividing in components or frequencies
#endif
#endif

            ! Store absorption profile
            data2(iil:iil+nf,1) = etmp0(if0l:if1l)

            ! If there is no MPI, normalize it here
            if (.not.MPID%mpi) then
              jjl = iil
              daux = data2(jjl,1)*Atom(ia)%W0(jtran)
              do ifreq=if0l+1,if1l-1
                jjl = jjl + 1
                daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)
              end do
              if (if1l.gt.if0l) then
                jjl = jjl + 1
                daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)
              end if

              if (daux.gt.0d0) then
                daux = 1d0/daux
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux
              end if
            end if


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

              ! Store emission profile
              data2(iil:iil+nf,2) = estmp0(if0l:if1l)

              ! If there is no MPI, normalize it here
              if (.not.MPID%mpi) then
                jjl = iil
                daux = data2(jjl,2)*Atom(ia)%W0(jtran)
                do ifreq=if0l+1,if1l-1
                  jjl = jjl + 1
                  daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)
                end do
                if (if1l.gt.if0l) then
                  jjl = jjl + 1
                  daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)
                end if

                if (daux.gt.0d0) then
                  daux = 1d0/daux
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux
                end if
              end if
            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Trano index
              ltran = Atom(ia)%itrano(jtran)

              ! Frec index
              indxf = Frec%indx(jtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (PRAM) then
                indx = Red%indx(jtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ! If from file, point
                if (vpfil) &
                  p_prof => prof(if0l2:if1l2, &
                                 1:Atom(ia)%i_Vind(jtran)%ncom)

                t0 = Atom(ia)%tshift + 1
                t1 = t0 + Atom(ia)%ntran - 1

                ! Angle-average
                if (AV) then
                  call emiss2ord_AA(Atom(ia),Geom,Atmo%vx(iz), &
                                    Atmo%vy(iz),Atmo%vz(iz), &
                                    Frec%omega,p_red,p_frec, &
                                    Flgsg,p_Norm,jtran, &
                                    itermu,iterml,iz, &
                                    if0l2,if1l2,DwT,Dw,vfac, &
                                    Bfield,Atmo%vmi(iz),TBo, &
                                    Stokes,JKQa,JKQ(:,:,t0:t1), &
                                    JKQC,p_prof, &
                                    es2tmp0(if0l2:if1l2), &
                                    es2tmp1(if0l2:if1l2), &
                                    es2tmp2(if0l2:if1l2), &
                                    es2tmp3(if0l2:if1l2))
                else
                  call emiss2ord_AD(Atom(ia),Geom,.False., &
                                    Atmo%vx(iz),Atmo%vy(iz), &
                                    Atmo%vz(iz),Frec%omega,p_red, &
                                    p_frec,Flgsg,p_Norm,jdir,jtran, &
                                    itermu,iterml,iz,if0l2,if1l2, &
                                    DwT,Dw,vfac,Atmo%vmi(iz),TBo, &
                                    Stokes,JKQ(:,:,t0:t1),p_prof, &
                                    es2tmp0(if0l2:if1l2), &
                                    es2tmp1(if0l2:if1l2), &
                                    es2tmp2(if0l2:if1l2), &
                                    es2tmp3(if0l2:if1l2))
                end if

                !
                ! Total emissivity
                !
                estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
                                       estmp0(if0l2:if1l2)
                estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
                                       estmp1(if0l2:if1l2)
                estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
                                       estmp2(if0l2:if1l2)
                estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
                                       estmp3(if0l2:if1l2)
              end if
            end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        ! No magnetic field
        else

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)
            nf = if1l - if0l

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncomNB
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l, &
                             1:Atom(ia)%i_Vind(jtran)%ncomNB)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

            !
            ! Get first order RT coefficients
            !
            call rt1ordNB(Atom(ia),TBo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK,p_prof, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l), &
                          estmp0(if0l:if1l),estmp1(if0l:if1l), &
                          estmp2(if0l:if1l),estmp3(if0l:if1l), &
                          rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                          rstmp3(if0l:if1l))

            ! Store absorption profile
            data2(iil:iil+nf,1) = etmp0(if0l:if1l)

            ! If there is no MPI, normalize it here
            if (.not.MPID%mpi) then
              jjl = iil
              daux = data2(jjl,1)*Atom(ia)%W0(jtran)
              do ifreq=if0l+1,if1l-1
                jjl = jjl + 1
                daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)
              end do
              if (if1l.gt.if0l) then
                jjl = jjl + 1
                daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)
              end if

              if (daux.gt.0d0) then
                daux = 1d0/daux
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux
              end if
            end if


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

              ! Store emission profile
              data2(iil:iil+nf,2) = estmp0(if0l:if1l)

              ! If there is no MPI, normalize it here
              if (.not.MPID%mpi) then
                jjl = iil
                daux = data2(jjl,2)*Atom(ia)%W0(jtran)
                do ifreq=if0l+1,if1l-1
                  jjl = jjl + 1
                  daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)
                end do
                if (if1l.gt.if0l) then
                  jjl = jjl + 1
                  daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)
                end if

                if (daux.gt.0d0) then
                  daux = 1d0/daux
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux
                end if
              end if
            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Trano index
              ltran = Atom(ia)%itrano(jtran)

              ! Frec index
              indxf = Frec%indx(jtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (PRAM) then
                indx = Red%indx(jtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ! If from file, point
                if (vpfil) &
                  p_prof => prof(if0l2:if1l2, &
                                 1:Atom(ia)%i_Vind(jtran)%ncomNB)

                t0 = Atom(ia)%tshift + 1
                t1 = t0 + Atom(ia)%ntran - 1

                ! Angle-averaged
                if (AV) then
                  call emiss2ordNB_AA(Atom(ia),Geom,Atmo%vx(iz), &
                                      Atmo%vy(iz),Atmo%vz(iz), &
                                      Frec%omega,p_red,p_frec, &
                                      Flgsg,p_Norm,jtran, &
                                      itermu,iterml,iz, &
                                      if0l2,if1l2,DwT,Dw,vfac, &
                                      Atmo%vmi(iz),TBo,Stokes,JKQa, &
                                      JKQ(:,:,t0:t1),JKQC,p_prof, &
                                      es2tmp0(if0l2:if1l2), &
                                      es2tmp1(if0l2:if1l2), &
                                      es2tmp2(if0l2:if1l2), &
                                      es2tmp3(if0l2:if1l2))
                else
                  call emiss2ordNB_AD(Atom(ia),Geom,.False., &
                                      Atmo%vx(iz),Atmo%vy(iz), &
                                      Atmo%vz(iz),Frec%omega,p_red, &
                                      p_frec,Flgsg,p_Norm,jdir, &
                                      jtran,itermu,iterml,iz, &
                                      if0l2,if1l2,DwT,Dw,vfac, &
                                      Atmo%vmi(iz),TBo,Stokes, &
                                      JKQ(:,:,t0:t1),p_prof, &
                                      es2tmp0(if0l2:if1l2), &
                                      es2tmp1(if0l2:if1l2), &
                                      es2tmp2(if0l2:if1l2), &
                                      es2tmp3(if0l2:if1l2))
                end if

                !
                ! Total emissivity
                !
                estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
                                       estmp0(if0l2:if1l2)
                estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
                                       estmp1(if0l2:if1l2)
                estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
                                       estmp2(if0l2:if1l2)
                estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
                                       estmp3(if0l2:if1l2)
              end if
            end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        end if ! Magnetic field


        !
        ! Photoionization
        !

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell


          !
          ! Emissivity
          !
          if (PIRAM) then

            call photoepsS(Atom(ia),Frec%omega3(if0l:if1l), &
                           Frec%exu(if0l:if1l,iz), &
                           Atmo%T(iz),Atmo%ne(iz),jtran, &
                           ilevelu,iz,if0l,if1l,estmp0(if0l:if1l), &
                           rstmp1(if0l:if1l))

          else

            call photoeps(Atom(ia),Frec%omega,Atmo%T(iz), &
                          Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                          estmp0(if0l:if1l),rstmp1(if0l:if1l))

          end if

          ! Add contribution to emissivity
          epsA(if0l:if1l,0) = estmp0(if0l:if1l) + epsA(if0l:if1l,0)


          !
          ! Absorptivity
          !
          call photoabs(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                        etmp0(if0l:if1l))

          ! Remove the stimulated part
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - rstmp1(if0l:if1l)

          ! Add contribution to absorptivity
          etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)

        end do ! b-f transitions

        ! Multipli by the population of the atom the contribution of
        ! this atom to the RT coefficients
        data1(0,:,0) = data1(0,:,0) + etaA(:,0)*Atom(ia)%n(iz)
        data1(1,:,0) = data1(1,:,0) + etaA(:,1)*Atom(ia)%n(iz)
        data1(0,:,4) = data1(0,:,4) + epsA(:,0)*Atom(ia)%n(iz)
        data1(1,:,4) = data1(1,:,4) + epsA(:,1)*Atom(ia)%n(iz)

        if (.not.RTaxial) then
          data1(2,:,0) = data1(2,:,0) + etaA(:,2)*Atom(ia)%n(iz)
          data1(3,:,0) = data1(3,:,0) + etaA(:,3)*Atom(ia)%n(iz)
          data1(2,:,1) = data1(2,:,1) + rhaA(:,3)*Atom(ia)%n(iz)
          data1(3,:,1) = data1(3,:,1) - rhaA(:,2)*Atom(ia)%n(iz)
          data1(3,:,2) = data1(3,:,2) + rhaA(:,1)*Atom(ia)%n(iz)
          data1(2,:,4) = data1(2,:,4) + epsA(:,2)*Atom(ia)%n(iz)
          data1(3,:,4) = data1(3,:,4) + epsA(:,3)*Atom(ia)%n(iz)
        end if

      end do ! Atoms


      !
      ! Transform into the data arrays
      !

      ! Absorption matrix
           !col!row
      data1(1,:,0) = data1(1,:,0)/(data1(0,:,0) + vacuum)

      ! Source function
      data1(0,:,4) = data1(0,:,4)/(data1(0,:,0) + vacuum)
      data1(1,:,4) = data1(1,:,4)/(data1(0,:,0) + vacuum)

      ! If not axial symmetry, opacity matrix and source functions
      if (.not.RTaxial) then

        ! Absorption matrix
             !col!row
        data1(2,:,0) = data1(2,:,0)/(data1(0,:,0) + vacuum)
        data1(3,:,0) = data1(3,:,0)/(data1(0,:,0) + vacuum)
        data1(2,:,1) = data1(2,:,1)/(data1(0,:,0) + vacuum)
        data1(3,:,1) = data1(3,:,1)/(data1(0,:,0) + vacuum)
        data1(3,:,2) = data1(3,:,2)/(data1(0,:,0) + vacuum)

        ! Source function
        data1(2,:,4) = data1(2,:,4)/(data1(0,:,0) + vacuum)
        data1(3,:,4) = data1(3,:,4)/(data1(0,:,0) + vacuum)

      end if ! axial symmetry

      ! Nullify pointers
      if (VPFIL) then
        deallocate(prof)
        nullify(prof)
        if (associated(p_prof)) nullify(p_prof)
      else
        if (associated(prof)) nullify(prof)
        deallocate(p_prof)
        nullify(p_prof)
      end if
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)

      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTCoeff'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      close(200)
      urou = 'RTCoeff'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTCoeff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the radiation transfer coefficients for emergence
      !! or emergence_serial\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!     JKQa(dcomplex(:,:,:)): Ad-hoc asymmetry\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!   Stokes(dfloat(:,:,:,:)): Stokes parameters\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients\n
      subroutine RTCoeffe(Frec,Red,Atom,LTElines,Atmo,Flgsg,Geom,iz, &
                          ith,iph,if0,if1,JKQa,JKQ,JKQC,cdir,Cont, &
                          Bfield,Stokes,data1)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(0:3,if0:if1,0:4):: data1
      complex(kind=8), dimension(:,:,:), intent(in):: JKQa
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(in):: JKQC

      ! Local

      integer:: iS,K,iQ,iterml,itermu,ia,jtran,ktran,ltran,icom
      integer:: jdir,jbdir,jcdir,jddir,icdir,ilevell,ilevelu,nodir
      integer:: if0l,if1l,if0l2,if1l2,t0,t1,offset,ios,indx,indxf

      double precision:: DwT,pE,absK,Dw,vfac,ct,st,cc,sc,loffset
      double precision, dimension(if0:if1,0:3):: etaA,epsA,rhsA,rhaA
      double precision, dimension(if0:if1):: etmp0,estmp0,es2tmp0
      double precision, dimension(if0:if1):: etmp1,estmp1, &
                                             rstmp1,rtmp1,es2tmp1
      double precision, dimension(if0:if1):: etmp2,estmp2, &
                                             rstmp2,rtmp2,es2tmp2
      double precision, dimension(if0:if1):: etmp3,estmp3, &
                                             rstmp3,rtmp3,es2tmp3
      double precision, dimension(if0:if1):: intgr

      complex(kind=8), dimension(0:3,-2:2,0:2):: TBo,TSo
      complex(kind=8), dimension(:,:), pointer:: prof,p_prof

      type(Frequencyc2_class), pointer:: p_frec
      type(Redc2_class), pointer:: p_red
      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(prof,p_prof,p_frec,p_red,p_Norm,p_LTEprof)


      ! Emerging ray, take LOS quantities
      jdir = Geom%i_geom(iph,ith)
      TBo = Geom%TBL(:,:,:,iz)
      TSo = Geom%TSL
      icdir = min(jdir,cdir)

      ! If there are dynamics
      if (dyn) then
        nodir = jdir
        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))
      else
        nodir = 1
      end if

      if (PRD) then
        jbdir = min(jdir,Frec%ndir)
        jcdir = min(jdir,Red%ndir)
        jddir = min(jdir,Red%njdir)
      end if

      ! If using file, nullify pointer to norm and allocate profile
      if (VPFIL) then

        p_Norm => Atom(1)%Normp(1,1,1)
        allocate(prof(if0:if1,Atom(1)%Mncom))

      ! If not using file
      else

        allocate(p_prof(1,1))

      end if


      !
      ! Initialize to 0 the parts without continuum
      !
      data1(1,:,0) = 0d0

      ! If not axial
      if (.not.RTaxial) then
        data1(2,:,0) = 0d0
        data1(3,:,0) = 0d0
        data1(2,:,1) = 0d0
        data1(3,:,1) = 0d0
        data1(3,:,2) = 0d0
      end if


      !
      ! Continuum contribution
      !

      ! Absorptivity
      data1(0,:,0) = Cont(:,1,icdir)

      ! For each Stokes parameter
      do iS=0,3

        ! if axial skip U and V
        if (RTaxial.and.iS.gt.1) cycle

        ! Reset integral
        intgr = .0D0

        !
        ! Compute the sum over K and Q of TKQ*JKQ(k)
        !

        ! For each K
        do K=0,Krad

          ! For each Q
          do iQ=-K,K

            ! Add contribution to the integral
            intgr = intgr + dble(TSo(iS,iQ,K)*JKQC(iQ,K,if0:if1))

          end do ! Q
        end do ! K

        ! Emissivity by scattering
        data1(iS,:,4) = intgr*Cont(:,2,icdir)

      end do ! Stokes parameters

      ! Add thermal emissivity
      data1(0,:,4) = data1(0,:,4) + cont(:,3,icdir)


      !
      ! Calculate the Doppler shift factor
      !

      ! If there are dynamics
      if (dyn) then

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If not, no need to calculate
      else

        vfac = 1d0

      end if ! dynamics

      !
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If this CPU does not have frequencies in this line, skip
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Photon energy (cgs) and convertion factor
        pE = convF*LTElines(ia)%Dfreq
        absK = 1d21*(2d0*c)*LTElines(ia)%Dfreq**2d0

        ! Add the microt. to Doppler width
        Dw = LTElines(ia)%Dfreq*sqrt(DwT*DwT + &
                                     Atmo%vmi(iz)**2d0)

        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Get prof
          if (LVPRAM) then
            p_LTEprof => LTElines(ia)%prof(iz,jdir)
          else
            p_LTEprof => LTElines(1)%prof(1,1)
          end if

          call rt1ordLTE(LTElines(ia),TBo,Frec%omega,Flgsg,iz, &
                         if0l,if1l,p_LTEprof,Dw,vfac, &
                         Bfield%Bstrength(iz),pE, &
                         etmp0(if0l:if1l),etmp1(if0l:if1l), &
                         etmp2(if0l:if1l),etmp3(if0l:if1l), &
                         rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                         rtmp3(if0l:if1l), &
                         estmp0(if0l:if1l),estmp1(if0l:if1l), &
                         estmp2(if0l:if1l),estmp3(if0l:if1l), &
                         rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                         rstmp3(if0l:if1l))

          !
          ! Stimulated emission contribution
          !

          ! If there is stimulated emission
          if (stm) then

            ! Correct for stimulated emission
            etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                               estmp0(if0l:if1l)/absK
            etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                               estmp1(if0l:if1l)/absK
            etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                               estmp2(if0l:if1l)/absK
            etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                               estmp3(if0l:if1l)/absK
            rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                               rstmp1(if0l:if1l)/absK
            rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                               rstmp2(if0l:if1l)/absK
            rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                               rstmp3(if0l:if1l)/absK

          endif ! Stimulated emission

          ! Multipli by the population of the atom the contribution of
          ! this atom to the RT coefficients
          data1(0,if0l:if1l,0) = data1(0,if0l:if1l,0) + &
                                 etmp0(if0l:if1l)*LTElines(ia)%n(iz)
          data1(1,if0l:if1l,0) = data1(1,if0l:if1l,0) + &
                                 etmp1(if0l:if1l)*LTElines(ia)%n(iz)
          data1(0,if0l:if1l,4) = data1(0,if0l:if1l,4) + &
                               estmp0(if0l:if1l)*pE*LTElines(ia)%n(iz)
          data1(1,if0l:if1l,4) = data1(1,if0l:if1l,4) + &
                               estmp1(if0l:if1l)*pE*LTElines(ia)%n(iz)

          if (.not.RTaxial) then
            data1(2,if0l:if1l,0) = data1(2,if0l:if1l,0) + &
                                   etmp2(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,0) = data1(3,if0l:if1l,0) + &
                                   etmp3(if0l:if1l)*LTElines(ia)%n(iz)
            data1(2,if0l:if1l,1) = data1(2,if0l:if1l,1) + &
                                   rtmp3(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,1) = data1(3,if0l:if1l,1) - &
                                   rtmp2(if0l:if1l)*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,2) = data1(3,if0l:if1l,2) + &
                                   rtmp1(if0l:if1l)*LTElines(ia)%n(iz)
            data1(2,if0l:if1l,4) = data1(2,if0l:if1l,4) + &
                               estmp2(if0l:if1l)*pE*LTElines(ia)%n(iz)
            data1(3,if0l:if1l,4) = data1(3,if0l:if1l,4) + &
                               estmp3(if0l:if1l)*pE*LTElines(ia)%n(iz)
          end if

        ! No magnetic field
        else

          ! Get prof
          if (LVPRAM) then
            p_LTEprof => LTElines(ia)%prof(iz,nodir)
          else
            p_LTEprof => LTElines(1)%prof(1,1)
          end if

          ! RT coeffs
          call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_LTEprof,Dw,vfac,pE, &
                          etmp0(if0l:if1l),estmp0(if0l:if1l))

          !
          ! Stimulated emission contribution
          !

          ! If there is stimulated emission
          if (stm) then

            ! Correct for stimulated emission
            etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                               estmp0(if0l:if1l)/absK

          endif ! Stimulated emission

          ! Multipli by the population of the atom the contribution of
          ! this atom to the RT coefficients
          data1(0,if0l:if1l,0) = data1(0,if0l:if1l,0) + &
                                 etmp0(if0l:if1l)*LTElines(ia)%n(iz)
          data1(0,if0l:if1l,4) = data1(0,if0l:if1l,4) + &
                               estmp0(if0l:if1l)*pE*LTElines(ia)%n(iz)

        end if ! Magnetic field

      end do ! LTE lines


      ! For each atom
      do ia=1,nA

        ! Initialize atomic RT coefficients
        epsA = .0D0
        rhsA = .0D0
        etaA = .0D0
        rhaA = .0D0

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))


        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(jdir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tBsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncom
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l,1:Atom(ia)%i_Vind(jtran)%ncom)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,jdir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)
#ifdef RDIPEV
            !
            ! First order RT coefficients
            !
            call rt1ord(Atom(ia),TBo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l), &
                        estmp0(if0l:if1l),estmp1(if0l:if1l), &
                        estmp2(if0l:if1l),estmp3(if0l:if1l), &
                        rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                        rstmp3(if0l:if1l))
#else
#ifdef _OPENMP
            ! If dividing in components
            if (Atom(ia)%omp_comp_1ord(jtran)) then

              !
              ! First order emissivity
              !
              call emiss_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                           jtran,itermu,iterml,iz,if0l,if1l, &
                           p_Norm,Dw,vfac,p_prof, &
                           Atom(ia)%omp_1c(jtran), &
                           estmp0(if0l:if1l),estmp1(if0l:if1l), &
                           estmp2(if0l:if1l),estmp3(if0l:if1l), &
                           rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                           rstmp3(if0l:if1l))

              !
              ! Absorptivity
              !
              call absorb_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            Atom(ia)%omp_1c(jtran), &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l))

            ! Dividing in frequencies
            else
#endif
            !
            ! First order emissivity
            !
            call emiss(Atom(ia),TBo,Frec%omega,Flgsg, &
                       jtran,itermu,iterml,iz,if0l,if1l, &
                       p_Norm,Dw,vfac,p_prof, &
                       estmp0(if0l:if1l),estmp1(if0l:if1l), &
                       estmp2(if0l:if1l),estmp3(if0l:if1l), &
                       rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                       rstmp3(if0l:if1l))

            !
            ! Absorptivity
            !
            call absorb(Atom(ia),TBo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l))
#ifdef _OPENMP
            end if ! Dividing in components or frequencies
#endif
#endif


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Trano index
              ltran = Atom(ia)%itrano(jtran)

              ! Frec index
              indxf = Frec%indx(jtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (PRAM) then
                indx = Red%indx(jtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ! If from file, point
                if (vpfil) &
                  p_prof => prof(if0l2:if1l2, &
                                 1:Atom(ia)%i_Vind(jtran)%ncom)

                t0 = Atom(ia)%tshift + 1
                t1 = t0 + Atom(ia)%ntran - 1

                ! Angle-average
                if (AV) then
                  call emiss2ord_AA(Atom(ia),Geom,Atmo%vx(iz), &
                                    Atmo%vy(iz),Atmo%vz(iz), &
                                    Frec%omega,p_red,p_frec, &
                                    Flgsg,p_Norm,jtran, &
                                    itermu,iterml,iz, &
                                    if0l2,if1l2,DwT,Dw,vfac, &
                                    Bfield,Atmo%vmi(iz),TBo, &
                                    Stokes,JKQa,JKQ(:,:,t0:t1), &
                                    JKQC,p_prof, &
                                    es2tmp0(if0l2:if1l2), &
                                    es2tmp1(if0l2:if1l2), &
                                    es2tmp2(if0l2:if1l2), &
                                    es2tmp3(if0l2:if1l2))
                else
                  call emiss2ord_AD(Atom(ia),Geom,.True., &
                                    Atmo%vx(iz),Atmo%vy(iz), &
                                    Atmo%vz(iz),Frec%omega,p_red, &
                                    p_frec,Flgsg,p_Norm,1,jtran, &
                                    itermu,iterml,iz,if0l2,if1l2, &
                                    DwT,Dw,vfac,Atmo%vmi(iz),TBo, &
                                    Stokes,JKQ(:,:,t0:t1),p_prof, &
                                    es2tmp0(if0l2:if1l2), &
                                    es2tmp1(if0l2:if1l2), &
                                    es2tmp2(if0l2:if1l2), &
                                    es2tmp3(if0l2:if1l2))
                end if

                !
                ! Total emissivity
                !
                estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
                                       estmp0(if0l2:if1l2)
                estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
                                       estmp1(if0l2:if1l2)
                estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
                                       estmp2(if0l2:if1l2)
                estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
                                       estmp3(if0l2:if1l2)
              end if
            end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

          end do ! b-b transitions

        ! No magnetic field
        else

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncomNB
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l, &
                             1:Atom(ia)%i_Vind(jtran)%ncomNB)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

            !
            ! Get first order RT coefficients
            !
            call rt1ordNB(Atom(ia),TBo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK,p_prof, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l), &
                          estmp0(if0l:if1l),estmp1(if0l:if1l), &
                          estmp2(if0l:if1l),estmp3(if0l:if1l), &
                          rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                          rstmp3(if0l:if1l))

            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Trano index
              ltran = Atom(ia)%itrano(jtran)

              ! Frec index
              indxf = Frec%indx(jtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (PRAM) then
                indx = Red%indx(jtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ! If from file, point
                if (vpfil) &
                  p_prof => prof(if0l2:if1l2, &
                                 1:Atom(ia)%i_Vind(jtran)%ncomNB)

                t0 = Atom(ia)%tshift + 1
                t1 = t0 + Atom(ia)%ntran - 1

                if (AV) then
                  call emiss2ordNB_AA(Atom(ia),Geom,Atmo%vx(iz), &
                                      Atmo%vy(iz),Atmo%vz(iz), &
                                      Frec%omega,p_red,p_frec, &
                                      Flgsg,p_Norm,jtran, &
                                      itermu,iterml,iz, &
                                      if0l2,if1l2,DwT,Dw,vfac, &
                                      Atmo%vmi(iz),TBo,Stokes,JKQa, &
                                      JKQ(:,:,t0:t1),JKQC,p_prof, &
                                      es2tmp0(if0l2:if1l2), &
                                      es2tmp1(if0l2:if1l2), &
                                      es2tmp2(if0l2:if1l2), &
                                      es2tmp3(if0l2:if1l2))
                else
                  call emiss2ordNB_AD(Atom(ia),Geom,.True., &
                                      Atmo%vx(iz),Atmo%vy(iz), &
                                      Atmo%vz(iz),Frec%omega,p_red, &
                                      p_frec,Flgsg,p_Norm,1, &
                                      jtran,itermu,iterml, &
                                      iz,if0l2,if1l2,DwT,Dw,vfac, &
                                      Atmo%vmi(iz),TBo,Stokes, &
                                      JKQ(:,:,t0:t1),p_prof, &
                                      es2tmp0(if0l2:if1l2), &
                                      es2tmp1(if0l2:if1l2), &
                                      es2tmp2(if0l2:if1l2), &
                                      es2tmp3(if0l2:if1l2))
                end if

                !
                ! Total emissivity
                !
                estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
                                       estmp0(if0l2:if1l2)
                estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
                                       estmp1(if0l2:if1l2)
                estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
                                       estmp2(if0l2:if1l2)
                estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
                                       estmp3(if0l2:if1l2)
              end if
            end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

          end do ! b-b transitions

        end if ! Magnetic field


        !
        ! Photoionization
        !

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell


          !
          ! Emissivity
          !
          if (PIRAM) then

            call photoepsS(Atom(ia),Frec%omega3(if0l:if1l), &
                           Frec%exu(if0l:if1l,iz), &
                           Atmo%T(iz),Atmo%ne(iz),jtran, &
                           ilevelu,iz,if0l,if1l,estmp0(if0l:if1l), &
                           rstmp1(if0l:if1l))

          else

            call photoeps(Atom(ia),Frec%omega,Atmo%T(iz), &
                          Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                          estmp0(if0l:if1l),rstmp1(if0l:if1l))

          end if

          ! Add contribution to emissivity
          epsA(if0l:if1l,0) = estmp0(if0l:if1l) + epsA(if0l:if1l,0)


          !
          ! Absorptivity
          !
          call photoabs(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                        etmp0(if0l:if1l))

          ! Remove the stimulated part
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - rstmp1(if0l:if1l)

          ! Add contribution to absorptivity
          etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)

        end do ! b-f transitions

        ! Multipli by the population of the atom the contribution of
        ! this atom to the RT coefficients
        data1(0,:,0) = data1(0,:,0) + etaA(:,0)*Atom(ia)%n(iz)
        data1(1,:,0) = data1(1,:,0) + etaA(:,1)*Atom(ia)%n(iz)
        data1(0,:,4) = data1(0,:,4) + epsA(:,0)*Atom(ia)%n(iz)
        data1(1,:,4) = data1(1,:,4) + epsA(:,1)*Atom(ia)%n(iz)

        if (.not.RTaxial) then
          data1(2,:,0) = data1(2,:,0) + etaA(:,2)*Atom(ia)%n(iz)
          data1(3,:,0) = data1(3,:,0) + etaA(:,3)*Atom(ia)%n(iz)
          data1(2,:,1) = data1(2,:,1) + rhaA(:,3)*Atom(ia)%n(iz)
          data1(3,:,1) = data1(3,:,1) - rhaA(:,2)*Atom(ia)%n(iz)
          data1(3,:,2) = data1(3,:,2) + rhaA(:,1)*Atom(ia)%n(iz)
          data1(2,:,4) = data1(2,:,4) + epsA(:,2)*Atom(ia)%n(iz)
          data1(3,:,4) = data1(3,:,4) + epsA(:,3)*Atom(ia)%n(iz)
        end if

      end do ! Atoms


      !
      ! Transform into the data arrays
      !

      ! Absorption matrix
           !col!row
      data1(1,:,0) = data1(1,:,0)/(data1(0,:,0) + vacuum)

      ! Source function
      data1(0,:,4) = data1(0,:,4)/(data1(0,:,0) + vacuum)
      data1(1,:,4) = data1(1,:,4)/(data1(0,:,0) + vacuum)

      ! If not axial symmetry, opacity matrix and source functions
      if (.not.RTaxial) then

        ! Absorption matrix
             !col!row
        data1(2,:,0) = data1(2,:,0)/(data1(0,:,0) + vacuum)
        data1(3,:,0) = data1(3,:,0)/(data1(0,:,0) + vacuum)
        data1(2,:,1) = data1(2,:,1)/(data1(0,:,0) + vacuum)
        data1(3,:,1) = data1(3,:,1)/(data1(0,:,0) + vacuum)
        data1(3,:,2) = data1(3,:,2)/(data1(0,:,0) + vacuum)

        ! Source function
        data1(2,:,4) = data1(2,:,4)/(data1(0,:,0) + vacuum)
        data1(3,:,4) = data1(3,:,4)/(data1(0,:,0) + vacuum)

      end if ! axial symmetry

      ! Nullify pointers
      if (VPFIL) then
        deallocate(prof)
        nullify(prof)
        if (associated(p_prof)) nullify(p_prof)
      else
        if (associated(prof)) nullify(prof)
        deallocate(p_prof)
        nullify(p_prof)
      end if
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)

      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTCoeffe'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTCoeffe'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTCoeffe

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the intensity absorption coefficients\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!           etaI(dfloat(:)): Absorption coefficient
      subroutine RTAbs(Frec,Atom,LTElines,Atmo,Flgsg,Geom,iz, &
                       ith,iph,if0,if1,cdir,Cont,Bfield,etaI)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      type(Frequency_class), intent(in):: Frec
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(if0:if1), intent(out):: etaI

      ! Local

      integer:: iterml,itermu,ia,jtran,jdir,icdir,ilevell,ilevelu
      integer:: if0l,if1l,offset,icom,ios,nodir

      double precision:: DwT,pE,absK,Dw,vfac,ct,st,cc,sc,loffset
      double precision, dimension(if0:if1):: etaA
      double precision, dimension(if0:if1):: etmp0,estmp0
      double precision, dimension(if0:if1):: etmp1,estmp1, &
                                             rstmp1,rtmp1
      double precision, dimension(if0:if1):: etmp2,estmp2, &
                                             rstmp2,rtmp2
      double precision, dimension(if0:if1):: etmp3,estmp3, &
                                             rstmp3,rtmp3

      complex(kind=8), dimension(0:3,-2:2,0:2):: TBo,TSo
      complex(kind=8), dimension(:,:), pointer:: prof,p_prof

      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(prof,p_prof,p_Norm,p_LTEprof)


      ! Emerging ray, take LOS quantities
      jdir = Geom%i_geom(iph,ith)
      TBo = Geom%TBL(:,:,:,iz)
      TSo = Geom%TSL
      icdir = min(jdir,cdir)

      ! If there are dynamics
      if (dyn) then
        nodir = jdir
        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))
      else
        nodir = 1
      end if

      ! If using file, nullify pointer to norm and allocate profile
      if (VPFIL) then

        p_Norm => Atom(1)%Normp(1,1,1)
        allocate(prof(if0:if1,Atom(1)%Mncom))

      ! If not using file
      else

        allocate(p_prof(1,1))

      end if


      !
      ! Continuum contribution
      !

      ! If there is continuum
      etaI = Cont(:,1,icdir)


      !
      ! Calculate the Doppler shift factor
      !

      ! If there are dynamics
      if (dyn) then

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If not, no need to calculate
      else

        vfac = 1d0

      end if ! dynamics


      !
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If this CPU does not have frequencies in this line, skip
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Photon energy (cgs) and convertion factor
        pE = convF*LTElines(ia)%Dfreq
        absK = 1d21*(2d0*c)*LTElines(ia)%Dfreq**2d0

        ! Add the microt. to Doppler width
        Dw = LTElines(ia)%Dfreq*sqrt(DwT*DwT + &
                                     Atmo%vmi(iz)**2d0)

        !
        ! Check if magnetic field
        !

        ! If stimulated
        if (stm) then

          ! If magnetic field
          if (Bfield%Bstrength(iz).gt.TINYB) then

            ! Get prof
            if (LVPRAM) then
              p_LTEprof => LTElines(ia)%prof(iz,jdir)
            else
              p_LTEprof => LTElines(1)%prof(1,1)
            end if

            call rt1ordLTE(LTElines(ia),TBo,Frec%omega,Flgsg,iz, &
                           if0l,if1l,p_LTEprof,Dw,vfac, &
                           Bfield%Bstrength(iz),pE, &
                           etmp0(if0l:if1l),etmp1(if0l:if1l), &
                           etmp2(if0l:if1l),etmp3(if0l:if1l), &
                           rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                           rtmp3(if0l:if1l), &
                           estmp0(if0l:if1l),estmp1(if0l:if1l), &
                           estmp2(if0l:if1l),estmp3(if0l:if1l), &
                           rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                           rstmp3(if0l:if1l))

          ! No magnetic field
          else

            ! Get prof
            if (LVPRAM) then 
              p_LTEprof => LTElines(ia)%prof(iz,nodir)
            else
              p_LTEprof => LTElines(1)%prof(1,1)
            end if

            ! RT coeffs
            call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                            p_LTEprof,Dw,vfac,pE, &
                            etmp0(if0l:if1l),estmp0(if0l:if1l))

          end if ! Magnetic field

          ! Correct for stimulated emission
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                             estmp0(if0l:if1l)/absK
          etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                             estmp1(if0l:if1l)/absK
          etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                             estmp2(if0l:if1l)/absK
          etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                             estmp3(if0l:if1l)/absK
          rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                             rstmp1(if0l:if1l)/absK
          rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                             rstmp2(if0l:if1l)/absK
          rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                             rstmp3(if0l:if1l)/absK

        ! If no stimulated
        else

          ! If magnetic field
          if (Bfield%Bstrength(iz).gt.TINYB) then

            ! Get prof
            if (LVPRAM) then
              p_LTEprof => LTElines(ia)%prof(iz,jdir)
            else
              p_LTEprof => LTElines(1)%prof(1,1)
            end if

            call absorbLTE(LTElines(ia),TBo,Frec%omega,Flgsg,iz, &
                           if0l,if1l,p_LTEprof,Dw,vfac, &
                           Bfield%Bstrength(iz),pE, &
                           etmp0(if0l:if1l),etmp1(if0l:if1l), &
                           etmp2(if0l:if1l),etmp3(if0l:if1l), &
                           rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                           rtmp3(if0l:if1l))

          ! No magnetic field
          else

            ! Get prof
            if (LVPRAM) then
              p_LTEprof => LTElines(ia)%prof(iz,nodir)
            else
              p_LTEprof => LTElines(1)%prof(1,1)
            end if

            ! RT coeffs
            call absorbILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                            p_LTEprof,Dw,vfac,pE,etmp0(if0l:if1l))

          end if ! Magnetic field

        endif ! Stimulated emission

        ! Absorptivity
        etaI(if0l:if1l) = etaI(if0l:if1l) + etmp0(if0l:if1l)* &
                                            LTElines(ia)%n(iz)

      end do ! LTE lines

      ! For each atom
      do ia=1,nA

        ! Initialize atomic RT coefficients
        etaA = .0D0

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))


        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(jdir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tBsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncom
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l,1:Atom(ia)%i_Vind(jtran)%ncom)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,jdir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

#ifdef RDIPEV
            ! If stimulated
            if (stm) then

              !
              ! First order RT coefficients
              !
              call rt1ord(Atom(ia),TBo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK,p_prof, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l), &
                          estmp0(if0l:if1l),estmp1(if0l:if1l), &
                          estmp2(if0l:if1l),estmp3(if0l:if1l), &
                          rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                          rstmp3(if0l:if1l))

              ! Correction
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK

            ! Not stimulated
            else

#ifdef _OPENMP
              ! If dividing in components
              if (Atom(ia)%omp_comp_1ord(jtran)) then

                !
                ! Absorptivity
                !
                call absorb_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                              jtran,itermu,iterml,iz,if0l,if1l, &
                              p_Norm,Dw,vfac,absK,p_prof, &
                              Atom(ia)%omp_1c(jtran), &
                              etmp0(if0l:if1l),etmp1(if0l:if1l), &
                              etmp2(if0l:if1l),etmp3(if0l:if1l), &
                              rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                              rtmp3(if0l:if1l))

              ! Dividing in frequencies
              else
#endif
              !
              ! Absorptivity
              !
              call absorb(Atom(ia),TBo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK,p_prof, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l))
#ifdef _OPENMP
              end if ! Dividing in components or frequencies
#endif
            end if ! Stimulated

! No energy eigenbasis
#else
#ifdef _OPENMP
            ! If dividing in components
            if (Atom(ia)%omp_comp_1ord(jtran)) then

              !
              ! Absorptivity
              !
              call absorb_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            Atom(ia)%omp_1c(jtran), &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l))
            ! Dividing in frequencies
            else
#endif
            !
            ! Absorptivity
            !
            call absorb(Atom(ia),TBo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l))
#ifdef _OPENMP
            end if ! Dividing in components or frequencies
#endif

            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then
#ifdef _OPENMP
              ! If dividing in components
              if (Atom(ia)%omp_comp_1ord(jtran)) then

                !
                ! First order emissivity
                !
                call emiss_c(Atom(ia),TBo,Frec%omega,Flgsg, &
                             jtran,itermu,iterml,iz,if0l,if1l, &
                             p_Norm,Dw,vfac,p_prof, &
                             Atom(ia)%omp_1c(jtran), &
                             estmp0(if0l:if1l),estmp1(if0l:if1l), &
                             estmp2(if0l:if1l),estmp3(if0l:if1l), &
                             rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                             rstmp3(if0l:if1l))
              ! Dividing in frequencies
              else
#endif
              !
              ! First order emissivity
              !
              call emiss(Atom(ia),TBo,Frec%omega,Flgsg, &
                         jtran,itermu,iterml,iz,if0l,if1l, &
                         p_Norm,Dw,vfac,p_prof, &
                         estmp0(if0l:if1l),estmp1(if0l:if1l), &
                         estmp2(if0l:if1l),estmp3(if0l:if1l), &
                         rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                         rstmp3(if0l:if1l))
#ifdef _OPENMP
              end if ! Dividing in components or frequencies
#endif

                etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                   estmp0(if0l:if1l)/absK
            endif ! Stimulated emission
! energy eigenbasis
#endif

            ! Add the contribution to the absorptivity of this atom
            etaA(if0l:if1l) = etmp0(if0l:if1l) + etaA(if0l:if1l)

          end do ! b-b transitions

        ! No magnetic field
        else

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line,
            ! skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! If reading from file
            if (vpfil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(jtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              ! For each component
              do icom=1,Atom(ia)%i_Vind(jtran)%ncomNB
                offset = int(Atom(ia)%f0size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
                read(200, err=1100) prof(if0l:if1l,icom)
                offset = int(Atom(ia)%f1size(jtran))
                if (offset.gt.0) call fseek(200, offset, 1)
              end do

              close(200)

              ! Point
              p_prof => prof(if0l:if1l, &
                             1:Atom(ia)%i_Vind(jtran)%ncomNB)

            ! Not reading from file
            else

              ! Get norm
              p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

            end if

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

            ! If stimulated
            if (stm) then

              !
              ! Get first order RT coefficients
              !
              call rt1ordNB(Atom(ia),TBo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l), &
                            estmp0(if0l:if1l),estmp1(if0l:if1l), &
                            estmp2(if0l:if1l),estmp3(if0l:if1l), &
                            rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                            rstmp3(if0l:if1l))

              ! Correction
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK

            ! No stimulated
            else

              !
              ! Absorptivity
              !
              call absorbNB(Atom(ia),TBo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l))

            end if ! Stimulated

            ! Add the contribution to the absorptivity of this atom
            etaA(if0l:if1l) = etmp0(if0l:if1l) + etaA(if0l:if1l)

          end do ! b-b transitions

        end if ! Magnetic field



        !
        ! Photoionization
        !

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell


          !
          ! Emissivity
          !
          if (PIRAM) then

            call photoepsS(Atom(ia),Frec%omega3(if0l:if1l), &
                           Frec%exu(if0l:if1l,iz),Atmo%T(iz), &
                           Atmo%ne(iz),jtran,ilevelu,iz, &
                           if0l,if1l,estmp0(if0l:if1l), &
                           rstmp1(if0l:if1l))

          else

            call photoeps(Atom(ia),Frec%omega,Atmo%T(iz), &
                          Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                          estmp0(if0l:if1l),rstmp1(if0l:if1l))


          end if

          !
          ! Absorptivity
          !
          call photoabs(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                        etmp0(if0l:if1l))

          ! Remove the stimulated part
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - rstmp1(if0l:if1l)

          ! Add contribution to absorptivity
          etaA(if0l:if1l) = etmp0(if0l:if1l) + etaA(if0l:if1l)

        end do ! b-f transitions

        ! Multipli by the population of the atom the contribution of
        ! this atom to the RT coefficients
        etaI = etaI + etaA*Atom(ia)%n(iz)

      end do

      ! Nullify pointers
      if (VPFIL) then
        deallocate(prof)
        nullify(prof)
        if (associated(p_prof)) nullify(p_prof)
      else
        if (associated(prof)) nullify(prof)
        deallocate(p_prof)
        nullify(p_prof)
      end if
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)

      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTAbs'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTAbs'
      close(200)
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTAbs

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the radiation transfer coefficients for getRTCLE\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Geom(Geometry_class): Structure with quadrature data\n
      !!     GeomP(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!              TS(dcomplex): Geometrical tensor for the LOS in
      !!                            the vertical reference frame\n
      !!              TB(dcomplex): Geometrical tensor for the LOS in
      !!                            the magnetic reference frame\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!        spect(spect_class): Structure with the input spectra
      !!                            data\n
      !!         Cont(dfloat(:,:)): Background opacity data\n
      !!     add_cont_cle(logical): Include continuum contributions\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients
      subroutine RTCoeff_CLE(Frec,Atom,Atmo,MPID,Flgsg,Geom,GeomP, &
                             Bfield,TS,TB,if0,if1,JKQ,JKQC,spect, &
                             Cont,add_cont_cle,data1)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(inout):: Frec
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Coronapoint_class), intent(in):: GeomP
      type(Bfield_class), intent(in):: Bfield
      type(Spect_class), intent(in):: spect
      logical, intent(in):: add_cont_cle
      integer, intent(in):: if0,if1
      double precision, dimension(if0:if1,3):: Cont
      double precision, dimension(0:3,if0:if1,0:4):: data1
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TS,TB
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,if0:if1), intent(in):: JKQC

      ! Local

      integer:: iS,K,iQ,ia,ilevell,ilevelu,iterml,itermu,jtran,ktran
      integer:: ifreq,if0l,if1l!,if0l2,if1l2

      double precision:: vfac,DwT,Dw,pE,absK,ct,st,cc,sc
      double precision, dimension(if0:if1,0:3):: etaA,epsA,rhsA,rhaA
      double precision, dimension(if0:if1):: etmp0,estmp0!,es2tmp0
      double precision, dimension(if0:if1):: etmp1,estmp1, &
                                             rstmp1,rtmp1!,es2tmp1
      double precision, dimension(if0:if1):: etmp2,estmp2, &
                                             rstmp2,rtmp2!,es2tmp2
      double precision, dimension(if0:if1):: etmp3,estmp3, &
                                             rstmp3,rtmp3!,es2tmp3
      double precision, dimension(if0:if1):: intgr

      complex(kind=8), dimension(:,:), pointer:: p_prof

      type(Nindex_class), pointer:: p_Norm


      !
      ! Initialize
      !
      nullify(p_Norm,p_prof)

      ! Dummy
      allocate(p_prof(1,1))

      ! If there are dynamics
      if (dyn) then
        ct = cos(GeomP%geom(1))
        st = sqrt(1d0 - ct*ct)
        cc = cos(GeomP%geom(2))
        sc = sin(GeomP%geom(2))
      end if

      !
      ! Continuum contribution
      !

      ! If including continuum
      if (add_cont_cle) then

        ! Absorptivity
        data1(0,:,0) = Cont(:,1)

        ! For each Stokes parameter
        do iS=0,3

          ! Reset integral
          intgr = .0D0

          !
          ! Compute the sum over K and Q of TKQ*JKQ(k)
          !

          ! For each K
          do K=0,Krad

            ! For each Q
            do iQ=-K,K

              ! Add contribution to the integral
              intgr = intgr + dble(TS(iS,iQ,K)*JKQC(iQ,K,if0:if1))

            end do ! Q
          end do ! K

          ! Emissivity by scattering
          data1(iS,:,4) = intgr*Cont(:,2)

        end do ! Stokes parameters

        ! Add thermal emissivity
        data1(0,:,4) = data1(0,:,4) + cont(:,3)

      ! No continuum
      else

        ! Zero
        data1(0,:,0) = 0d0
        data1(:,:,4) = 0d0

      end if


      !
      ! Calculate the Doppler shift factor
      !

      ! If there are dynamics
      if (dyn) then

        vfac = 1d0 - Atmo%vx(1)*st*cc - atmo%vy(1)*st*sc - &
                     Atmo%vz(1)*ct

      ! If not, no need to calculate
      else

        vfac = 1d0

      end if ! dynamics

      !
      ! Fill exu
      !
      if (associated(Frec%exu)) then

        ! Initialize argument
        Frec%exu(:,1) = Frec%omega_ou(Frec%pif0:Frec%pif1)* &
                        (c2*1d4/Atmo%T(1))

        ! Get inverse exponential
        do ifreq=Frec%pif0,Frec%pif1

          ! Call math
          Frec%exu(ifreq,1) =  diexp(Frec%exu(ifreq,1))

        end do

      end if ! Filling exu

      ! For each atom
      do ia=1,nA

        ! Initialize atomic RT coefficients
        epsA = .0D0
        rhsA = .0D0
        etaA = .0D0
        rhaA = .0D0

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(1))

        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(1).gt.TINYB) then

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! Get norm
            p_Norm => Atom(ia)%Normp(jtran,1,1)

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(1)**2d0)
#ifdef RDIPEV
            !
            ! First order RT coefficients
            !
            call rt1ord(Atom(ia),TB,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,1,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l), &
                        estmp0(if0l:if1l),estmp1(if0l:if1l), &
                        estmp2(if0l:if1l),estmp3(if0l:if1l), &
                        rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                        rstmp3(if0l:if1l))
#else
#ifdef _OPENMP
            ! If dividing in components
            if (Atom(ia)%omp_comp_1ord(jtran)) then

              !
              ! First order emissivity
              !
              call emiss_c(Atom(ia),TB,Frec%omega_ou,Flgsg, &
                           jtran,itermu,iterml,1,if0l,if1l, &
                           p_Norm,Dw,vfac,p_prof, &
                           Atom(ia)%omp_1c(jtran), &
                           estmp0(if0l:if1l),estmp1(if0l:if1l), &
                           estmp2(if0l:if1l),estmp3(if0l:if1l), &
                           rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                           rstmp3(if0l:if1l))

              !
              ! Absorptivity
              !
              call absorb_c(Atom(ia),TB,Frec%omega_ou,Flgsg, &
                            jtran,itermu,iterml,1,if0l,if1l, &
                            p_Norm,Dw,vfac,absK,p_prof, &
                            Atom(ia)%omp_1c(jtran), &
                            etmp0(if0l:if1l),etmp1(if0l:if1l), &
                            etmp2(if0l:if1l),etmp3(if0l:if1l), &
                            rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                            rtmp3(if0l:if1l))

            ! Dividing in frequencies
            else
#endif

            !
            ! First order emissivity
            !
            call emiss(Atom(ia),TB,Frec%omega_ou,Flgsg, &
                       jtran,itermu,iterml,1,if0l,if1l, &
                       p_Norm,Dw,vfac,p_prof, &
                       estmp0(if0l:if1l),estmp1(if0l:if1l), &
                       estmp2(if0l:if1l),estmp3(if0l:if1l), &
                       rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                       rstmp3(if0l:if1l))

            !
            ! Absorptivity
            !
            call absorb(Atom(ia),TB,Frec%omega_ou,Flgsg, &
                        jtran,itermu,iterml,1,if0l,if1l, &
                        p_Norm,Dw,vfac,absK,p_prof, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l))
#ifdef _OPENMP
            end if ! Dividing in components or frequencies
#endif
#endif

            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
!           if(PRD.and.Atom(ia)%lemiss2(jtran))then

!             ! TODO
!             call emiss2ord(Atom(ia),Geom,.False., &
!                            Atmo%vx(iz),Atmo%vy(iz),Atmo%vz(iz), &
!                            Frec%omega,p_red, &
!                            Frec%stype(:,:,jddir),p_frec, &
!                            Frec%nth,Frec%nph,Frec%nfs(jddir), &
!                            Flgsg,p_Norm,jtran,itermu,iterml, &
!                            iz,iph,ith,if0l2,if1l2,DwT,Dw,vfac, &
!                            Bfield,Atmo%vmi(iz),TBo,Stokes,JKQa, &
!                            JKQ(:,:,t0:t1),JKQC,p_prof, &
!                            es2tmp0(if0l2:if1l2), &
!                            es2tmp1(if0l2:if1l2), &
!                            es2tmp2(if0l2:if1l2), &
!                            es2tmp3(if0l2:if1l2))

!             !
!             ! Total emissivity
!             !
!             estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
!                                    estmp0(if0l2:if1l2)
!             estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
!                                    estmp1(if0l2:if1l2)
!             estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
!                                    estmp2(if0l2:if1l2)
!             estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
!                                      estmp3(if0l2:if1l2)
!           end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

          end do ! b-b transitions

        ! No magnetic field
        else

          !
          ! Transition lines
          !

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! If this CPU does not have frequencies in this line, skip
            if (Atom(ia)%fflag(jtran)%absent) cycle

            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

            ! Get norm
            p_Norm => Atom(ia)%Normp(jtran,1,1)

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Photon energy (cgs) and convertion factor
            pE = convF*Atom(ia)%Dfreq(jtran)
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(1)**2d0)

            !
            ! Get first order RT coefficients
            !
            call rt1ordNB(Atom(ia),TB,Frec%omega_ou,Flgsg, &
                          jtran,itermu,iterml,1,if0l,if1l, &
                          p_Norm,Dw,vfac,absK,p_prof, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l), &
                          estmp0(if0l:if1l),estmp1(if0l:if1l), &
                          estmp2(if0l:if1l),estmp3(if0l:if1l), &
                          rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                          rstmp3(if0l:if1l))

            !
            ! First order emissivity
            !
           !call emissNB(Atom(ia),TB,Frec%omega_ou,Flgsg, &
           !             jtran,itermu,iterml,1,if0l,if1l, &
           !             p_Norm,Dw,vfac,p_prof, &
           !             estmp0(if0l:if1l),estmp1(if0l:if1l), &
           !             estmp2(if0l:if1l),estmp3(if0l:if1l), &
           !             rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
           !             rstmp3(if0l:if1l))

            !
            ! Absorptivity
            !
           !call absorbNB(Atom(ia),TB,Frec%omega_ou,Flgsg, &
           !              jtran,itermu,iterml,1,if0l,if1l, &
           !              p_Norm,Dw,vfac,absK,p_prof, &
           !              etmp0(if0l:if1l),etmp1(if0l:if1l), &
           !              etmp2(if0l:if1l),etmp3(if0l:if1l), &
           !              rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
           !              rtmp3(if0l:if1l))

            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK
              etmp1(if0l:if1l) = etmp1(if0l:if1l) - &
                                 estmp1(if0l:if1l)/absK
              etmp2(if0l:if1l) = etmp2(if0l:if1l) - &
                                 estmp2(if0l:if1l)/absK
              etmp3(if0l:if1l) = etmp3(if0l:if1l) - &
                                 estmp3(if0l:if1l)/absK
              rtmp1(if0l:if1l) = rtmp1(if0l:if1l) - &
                                 rstmp1(if0l:if1l)/absK
              rtmp2(if0l:if1l) = rtmp2(if0l:if1l) - &
                                 rstmp2(if0l:if1l)/absK
              rtmp3(if0l:if1l) = rtmp3(if0l:if1l) - &
                                 rstmp3(if0l:if1l)/absK

            endif ! Stimulated emission

            ! Add the contribution to the absorptivity and dispersion
            ! of this atom
            etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)
            etaA(if0l:if1l,1) = etmp1(if0l:if1l) + etaA(if0l:if1l,1)
            etaA(if0l:if1l,2) = etmp2(if0l:if1l) + etaA(if0l:if1l,2)
            etaA(if0l:if1l,3) = etmp3(if0l:if1l) + etaA(if0l:if1l,3)
            rhaA(if0l:if1l,1) = rtmp1(if0l:if1l) + rhaA(if0l:if1l,1)
            rhaA(if0l:if1l,2) = rtmp2(if0l:if1l) + rhaA(if0l:if1l,2)
            rhaA(if0l:if1l,3) = rtmp3(if0l:if1l) + rhaA(if0l:if1l,3)

            !
            ! Second order emissivity
            !
!           if(PRD.and.Atom(ia)%lemiss2(jtran))then

!             ! TODO
!             call emiss2ordNB(Atom(ia),Geom,.False.,Atmo%vx(iz), &
!                              Atmo%vy(iz),Atmo%vz(iz), &
!                              Frec%omega,p_red, &
!                              Frec%stype(:,:,jddir),p_frec, &
!                              Frec%nth,Frec%nph,Frec%nfs(jddir), &
!                              Flgsg,p_Norm, &
!                              jtran,itermu,iterml, &
!                              iz,iph,ith,if0l2,if1l2,DwT,Dw,vfac, &
!                              Atmo%vmi(iz),TBo,Stokes,JKQa, &
!                              JKQ(:,:,t0:t1),JKQC,p_prof, &
!                              es2tmp0(if0l2:if1l2), &
!                              es2tmp1(if0l2:if1l2), &
!                              es2tmp2(if0l2:if1l2), &
!                              es2tmp3(if0l2:if1l2))

!             !
!             ! Total emissivity
!             !
!             estmp0(if0l2:if1l2) = es2tmp0(if0l2:if1l2) + &
!                                    estmp0(if0l2:if1l2)
!             estmp1(if0l2:if1l2) = es2tmp1(if0l2:if1l2) + &
!                                    estmp1(if0l2:if1l2)
!             estmp2(if0l2:if1l2) = es2tmp2(if0l2:if1l2) + &
!                                    estmp2(if0l2:if1l2)
!             estmp3(if0l2:if1l2) = es2tmp3(if0l2:if1l2) + &
!                                    estmp3(if0l2:if1l2)
!           end if

            ! Add the contribution to the emissivity and dispersion
            ! of this atom
            epsA(if0l:if1l,0) = estmp0(if0l:if1l)*pE + &
                                epsA(if0l:if1l,0)
            epsA(if0l:if1l,1) = estmp1(if0l:if1l)*pE + &
                                epsA(if0l:if1l,1)
            epsA(if0l:if1l,2) = estmp2(if0l:if1l)*pE + &
                                epsA(if0l:if1l,2)
            epsA(if0l:if1l,3) = estmp3(if0l:if1l)*pE + &
                                epsA(if0l:if1l,3)

          end do ! b-b transitions

        end if ! Magnetic field


        !
        ! Photoionization
        !

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell


          !
          ! Emissivity
          !
          call photoepsS(Atom(ia),Frec%omega3_ou(if0l:if1l), &
                         Frec%exu(if0l:if1l,1), &
                         Atmo%T(1),Atmo%ne(1),jtran, &
                         ilevelu,1,if0l,if1l,estmp0(if0l:if1l), &
                         rstmp1(if0l:if1l))

          ! Add contribution to emissivity
          epsA(if0l:if1l,0) = estmp0(if0l:if1l) + epsA(if0l:if1l,0)


          !
          ! Absorptivity
          !
          call photoabs(Atom(ia),jtran,ilevell,1,if0l,if1l, &
                        etmp0(if0l:if1l))

          ! Remove the stimulated part
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - rstmp1(if0l:if1l)

          ! Add contribution to absorptivity
          etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)

        end do ! b-f transitions

        ! Multipli by the population of the atom the contribution of
        ! this atom to the RT coefficients
        data1(0,:,0) = data1(0,:,0) + etaA(:,0)*Atom(ia)%n(1)
        data1(1,:,0) = data1(1,:,0) + etaA(:,1)*Atom(ia)%n(1)
        data1(0,:,4) = data1(0,:,4) + epsA(:,0)*Atom(ia)%n(1)
        data1(1,:,4) = data1(1,:,4) + epsA(:,1)*Atom(ia)%n(1)
        data1(2,:,0) = data1(2,:,0) + etaA(:,2)*Atom(ia)%n(1)
        data1(3,:,0) = data1(3,:,0) + etaA(:,3)*Atom(ia)%n(1)
        data1(2,:,1) = data1(2,:,1) + rhaA(:,3)*Atom(ia)%n(1)
        data1(3,:,1) = data1(3,:,1) - rhaA(:,2)*Atom(ia)%n(1)
        data1(3,:,2) = data1(3,:,2) + rhaA(:,1)*Atom(ia)%n(1)
        data1(2,:,4) = data1(2,:,4) + epsA(:,2)*Atom(ia)%n(1)
        data1(3,:,4) = data1(3,:,4) + epsA(:,3)*Atom(ia)%n(1)

      end do ! Atoms


      !
      ! Transform into the data arrays
      !

      ! Absorption matrix
           !col!row
      data1(1,:,0) = data1(1,:,0)/(data1(0,:,0) + vacuum)
      data1(2,:,0) = data1(2,:,0)/(data1(0,:,0) + vacuum)
      data1(3,:,0) = data1(3,:,0)/(data1(0,:,0) + vacuum)
      data1(2,:,1) = data1(2,:,1)/(data1(0,:,0) + vacuum)
      data1(3,:,1) = data1(3,:,1)/(data1(0,:,0) + vacuum)
      data1(3,:,2) = data1(3,:,2)/(data1(0,:,0) + vacuum)

      ! Source function
      data1(0,:,4) = data1(0,:,4)/(data1(0,:,0) + vacuum)
      data1(1,:,4) = data1(1,:,4)/(data1(0,:,0) + vacuum)
      data1(2,:,4) = data1(2,:,4)/(data1(0,:,0) + vacuum)
      data1(3,:,4) = data1(3,:,4)/(data1(0,:,0) + vacuum)

      ! Nullify pointers
      deallocate(p_prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)

      return

      end subroutine RTCoeff_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeff_mod
