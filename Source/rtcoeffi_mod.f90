      !> Compute radiation transfer coefficients for NLTE problem of
      !! the first kind
      module rtcoeffi_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     04/20/2017
!  Last version:
!     10/31/2023 V3.0.8
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     10/31/2023:    V3.0.8 - Split the calls to emissI2ord into AA
!                             and AD (TdPA)
!                           - The index for Red%indx should be for
!                             direction index jddir, and not the
!                             jcdir one (TdPA)
!
!     10/16/2023:    V3.0.7 - Made LTElines allocatable to satisfy
!                             memory warnings (TdPA)
!
!     10/04/2023:    V3.0.6 - Bugfix: When not keeping the LTE lines
!                             profiles, point to a dummy (TdPA)
!
!     09/25/2023:    V3.0.5 - Always trim Voigt p. file name (TdPA)
!
!     08/17/2023:    V3.0.4 - Added pointer cleaning when there are
!                             file errors (TdPA)
!
!     08/07/2023:    V3.0.3 - Added the contribution of LTE lines to
!                             the RT coefficients (TdPA)
!                           - Combined together the absorbI and
!                             emissI calls, as there is no need to
!                             compute profiles twice (TdPA)
!
!     10/25/2022:    V3.0.2 - Nullify pointers as when starting each
!                             routine (TdPA)
!                           - Clean pointers at exit (TdPA)
!
!     07/27/2022:    V3.0.1 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.0 - To implement the 1.5D case Atmo%v has
!                             changed to Atmo%vx,%vy, and %vz (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!                           - Added OpenMP (TdPA)
!
!     09/11/2020:    V1.7.6 - Changes to accomodate the new Frec and
!                             Red structures (TdPA)
!
!     07/22/2020:    V1.7.5 - Bugfix: The modification in V1.7.4 did
!                             not include Termprof (TdPA)
!
!     06/02/2020:    V1.7.4 - Limited size of directional dimension
!                             of the normalization array in cases
!                             where it is possible (TdPA)
!
!     03/18/2020:    V1.7.3 - Bugfix: In rtcoeffie, the extra
!                             directional indexes were being computed
!                             even in CRD, accessing an undefined
!                             memory location (TdPA)
!
!     12/10/2019:    V1.7.2 - Now only passes to emissI2ord the JKQC
!                             that will be needed (TdPA)
!
!     11/19/2019:    V1.7.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!                           - Added stored offsets in frequency for
!                             the use of Voigt files (TdPA)
!
!     11/13/2019:    V1.7.0 - Now admits Voigt profiles stored in a
!                             file (TdPA)
!                           - Photoionization quantities can be in
!                             RAM without storing Voigt profiles for
!                             intensity (TdPA)
!
!     10/03/2019:    V1.6.2 - trano is now indexed in Frec and Red
!                             structures (TdPA)
!
!     08/08/2019:    V1.6.1 - Now it can use directly the stored
!                             profiles for the data2 absorption array,
!                             and the data2 emission is copied (TdPA)
!
!     06/04/2019:    V1.6.0 - Simplified the iterative scheme (TdPA)
!                           - Now there is no search of FS transitions
!                             because they are indexed (TdPA)
!
!     05/31/2019:    V1.5.0 - Changed the dimensionality of the
!                             profile and ratio variableis. Now it
!                             runs sequentially on atoms, transitions
!                             and frequencies to save memory and
!                             reduce the size of data shared through
!                             MPI messages (TdPA)
!                           - No check for emergence in TermProf,
!                             because it never is (TdPA)
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
!     03/12/2019:    V1.3.1 - Removed iexu from photoepsIS (TdPA)
!                           - Removed searches by looping (TdPA)
!                           - Added a very tiny denominator to avoid
!                             dividing by zero (TdPA)
!
!     02/20/2019:    V1.3.0 - Using specific TINY variables (TdPA)
!                           - This routine does not need the
!                             aborted dependence (TdPA)
!
!     09/04/2018:    V1.2.1 - Frec%iexu is passed to photoepsIS (TdPA)
!
!     08/06/2018:    V1.2.0 - Changed how the argument Atom%Norm is
!                             passed in every b-b routine (TdPA)
!                           - Split the call of photoeps depending on
!                             what is stored in RAM (TdPA)
!
!     05/16/2018:    V1.1.2 - Changed the azimuthal dimension of
!                             Stokes from nPh2 to nPh (TdPA)
!                           - Added stype variable to the call to
!                             emissI2ord (TdPA)
!
!     10/04/2017:    V1.1.1 - The absorption profile that goes into
!                             the integral is not corrected by stim.
!                             emission (TdPA)
!
!     10/03/2017:    V1.1.0 - Non magnetic case in Termprof (TdPA)
!
!     09/08/2017:   V1.0.12 - Avoid divisions by 0 (TdPA)
!                           - emiss2ord does not need jdir (TdPA)
!
!     08/21/2017:   V1.0.11 - Bugfix: if no PRD, jbdir and jcdir
!                             should not be defined (TdPA)
!
!     06/28/2017:   V1.0.10 - Receives and passes Red (TdPA)
!                           - Passing ia, jdir ,jbdir and jcdir into
!                             emiss2 (TdPA)
!                           - Passing the full trano to emiss2 (TdPA)
!                           - The limits of emiss2 now are in other
!                             location in Frec (TdPA)
!
!     06/22/2017:    V1.0.9 - Added Termprof (TdPA)
!
!     06/20/2017:    V1.0.8 - Bugfix: Integral of profile should not
!                             touch boundaries in the loop (TdPA)
!                           - The PRD factor in rLine was being
!                             applied for all lines given that at
!                             least one was in PRD (nosense), and
!                             the frequency variables in emiss2ord
!                             call had the full limits (TdPA)
!
!     06/19/2017:    V1.0.7 - Data stored in output directly (TdPA)
!
!     06/12/2017:    V1.0.6 - The frequency limits are inputs (TdPA)
!                           - The is no Cont%l anymore (TdPA)
!                           - Passing transition limits into auxiliar
!                             routines (TdPA)
!                           - Stimulated emission is 1st order (TdPA)
!
!     06/09/2017:    V1.0.5 - Introduced RTCoeffJ (TdPA)
!
!     06/08/2017:    V1.0.4 - Removed ALI flag (TdPA)
!                           - Added ALI for second order (TdPA)
!
!     06/06/2017:    V1.0.3 - The stimulated emission correction in
!                             the lambda operator was always taken
!                             into account, even when it was off for
!                             the absorptivity and the statistical
!                             equilibrium equations (TdPA)
!
!     05/12/2017:    V1.0.2 - Now for PRD the input frequency goes
!                             with FS transitions, not with term
!                             transitions (TdPA)
!                           - Added ALI flag (TdPA)
!
!     05/03/2017:    V1.0.1 - data1 was defined as 0:3, instead of
!                             0:2 (TdPA)
!
!     04/20/2017:    V1.0.0 - First version (TdPA)
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
!  RTCoeffI
!    This routine calculates the RT coefficients (only intensity)
!
!  RTCoeffI
!    This routine calculates the RT coefficients for emergence LOS
!    (only intensity)
!
!  RTAbsI
!    Calculates just intensity absorptivity
!
!  RTCoeffJ
!    This routine calculates the RT coefficients (only continuum)
!
!  Termprof
!    This routine calculates the RT coefficients (only intensity but
!    with full Stokes routines)
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use omp_mod
      use parameters_mod , only : convF , c , vacuum , TINYB , &
                                  VTINY , IPI
      use rtcoeffaux_mod
      use rtcoeffiaux_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the radiation transfer coefficients for solveri or
      !! solveri_serial\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!            JKQ(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         JKQC(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!     Stokes(dfloat(:,:,:)): Intensity\n
      !!             rl(dfloat(:)): Bound-bound transition strengths\n
      !!             rp(dfloat(:)): Bound-free transition strengths\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients\n
      !!        data2(dfloat(:,:)): Line profiles
      subroutine RTCoeffI(Frec,Red,Atom,LTElines,Atmo,MPID,Geom,iz, &
                          ith,iph,if0,if1,JKQ,JKQC,cdir,Cont, &
                          Stokes,rl,rp,data1,data2)

      ! I/O

      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(inout):: Red
      type(MPI_class), intent(in):: MPID
      type(Geometry_class), intent(in):: Geom
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(nxt), intent(in):: JKQ
      double precision, dimension(nfreq), intent(in):: JKQC
      double precision, dimension(:):: rl
      double precision, dimension(:):: rp
      double precision, dimension(if0:if1,0:1):: data1
      double precision, dimension(:,:):: data2

      ! Local

      integer:: iterml,itermu,iJl,iJu,ilevell,ilevelu,ifreq,nf2
      integer:: ia,jtran,ktran,fjtran,ffktran,ffjtran,ffltran
      integer:: jcdir,icdir,jddir,if0l,if1l,if0l2,if1l2,t0,t1,nodir
      integer:: iil,jjl,iip,jdir,jbdir,nf,offset,ios,ggf0,ggf1
      integer:: indx,indxf

      double precision:: daux,DwT,Dfreq,pE,absK,Dw,vfac
      double precision:: ct,st,cc,sc,pop
      double precision:: rhou,loffset
      double precision, dimension(if0:if1):: rpf,etmp
      double precision, dimension(if0:if1):: estmp,rstmp,es2tmp
      double precision, dimension(if0:if1):: prof

      ! Pointer
      type(Frequencyc2_class), pointer:: p_frec
      type(Redc2_class), pointer:: p_red
      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(p_frec,p_red,p_Norm,p_LTEprof)


      !
      ! Get the direction index for the continuum and calculate
      ! the Doppler shift if necessary
      !
      jdir = Geom%i_geom(iph,ith)
      if (PRD) then
        jbdir = min(jdir,Frec%ndir)
        jcdir = min(jdir,Red%ndir)
        jddir = min(jdir,Red%njdir)
      end if

      ! If there are dynamics
      if (dyn) then

        icdir = min(jdir,cdir)

        nodir = jdir

        ct = Geom%V_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = Geom%v_mux(iph)
        sc = Geom%v_muy(iph)*sqrt(1d0 - cc*cc)

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If there are not dynamics
      else

        icdir = 1
        nodir = 1

        vfac = 1d0

      end if

      ! If using file, nullify pointer
      if (VIFIL) p_Norm => Atom(1)%Normp(1,1,1)

      !
      ! Continuum contribution
      !

      ! Absorptivity
      data1(:,0) = Cont(:,1,icdir)

      ! Emissivity
      data1(:,1) = Cont(:,3,icdir) + Cont(:,2,icdir)*JKQC(if0:if1)

      !
      ! Calculate RT coefficients
      !

      ! LTE lines
      do ia=1,nLTEl

        ! Skip if absent
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = LTElines(ia)%n(iz)

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Get Norm
        if (LVIRAM) then
          p_LTEprof => LTElines(ia)%prof(iz,nodir)
        else
          p_LTEprof => LTElines(1)%prof(1,1)
        end if

        ! Get frequency of FS transition
        Dfreq = LTElines(ia)%Eu - LTElines(ia)%El

        ! Photon energy (cgs) and convertion factor
        pE = convF*Dfreq
        absK = 1d21*(2d0*c)*Dfreq**2d0

        ! Add the microturbulence to Doppler width
        Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        ! RT coeffs
        call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                        p_LTEprof,Dw,vfac,pE, &
                        etmp(if0l:if1l),estmp(if0l:if1l))

        !
        ! Stimulated emission contribution
        !
!$omp parallel default(none) &
!$omp shared(stm,etmp,if0l,if1l,estmp,absK,data1,pop)
        ! If there is stimulated emission
        if (stm) then
!$omp workshare
          etmp(if0l:if1l) = etmp(if0l:if1l) - &
                            estmp(if0l:if1l)/absK
!$omp end workshare
        end if

        ! Add the contribution to the absorptivity of this atom
!$omp workshare
        data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                             etmp(if0l:if1l)*pop
!$omp end workshare
!$omp end parallel

        ! Add the contribution to the emissivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(estmp,data1,if0l,if1l,pE,pop)
        data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                             estmp(if0l:if1l)*(pE*pop)
!$omp end parallel workshare

      end do ! LTE lines


      ! Initialize profile coefficient
      iil = 1
      iip = 1

      !
      ! For each atom
      do ia=1,nA

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = Atom(ia)%n(iz)


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

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          iterml = Atom(ia)%fst(jtran)%iterml

          ! Get Norm
          if (.not.VIFIL) p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

            ffktran = ffjtran + Atom(ia)%tfshift

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJl,iterml)

            ! Photon energy (cgs) and convertion factor
            pE = convF*Dfreq
            absK = 1d21*(2d0*c)*Dfreq**2d0

            ! Add the microturbulence to Doppler width
            Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


            !
            ! If reading from file
            !

            if (vifil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(ffjtran) + &
                        Atom(ia)%f0size(ffjtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              read(200, err=1100) prof(if0l:if1l)

              close(200)

            end if


            !
            ! RT coeffs
            !
            call rt1ordI(Atom(ia),Frec%omega,jtran,fjtran, &
                         itermu,iterml,iJu,iJl,iz, &
                         if0l,if1l,p_Norm,Dw,vfac,pE, &
                         prof(if0l:if1l),etmp(if0l:if1l), &
                         estmp(if0l:if1l),rhou)


            !
            ! Store absorption profile
            !

            ! If in file
            if (vifil) then

!$omp parallel workshare default(none) &
!$omp shared(data2,iil,nf,prof,if0l,if1l,Dw)
              data2(iil:iil+nf,1) = prof(if0l:if1l)*1d-5*sqrt(IPI)/Dw
!$omp end parallel workshare

            ! Storing
            else if (p_Norm%VRAM) then

!$omp parallel workshare default(none) &
!$omp shared(data2,iil,nf,p_Norm,fjtran,Dw)
              data2(iil:iil+nf,1) = p_Norm%prof(fjtran,1,1,1)%p* &
                                    1d-5*sqrt(IPI)/Dw
!$omp end parallel workshare

            ! Computing
            else

              ! If MPI
              if (MPID%mpi) then

                ! Copy absorbtivity in data2
!$omp parallel workshare default(none) &
!$omp shared(data2,iil,nf,etmp,if0l,if1l)
                data2(iil:iil+nf,1) = etmp(if0l:if1l)
!$omp end parallel workshare

              ! If there is no MPI, normalize it here
              else

                ! Initialize
                daux = 0d0

!$omp parallel default(none) &
!$omp private(jjl,ifreq) &
!$omp shared(data2,iil,nf,etmp,if0l,if1l,MPID,Atom,Frec,ia,jtran) &
!$omp reduction(+: daux)
                ! Copy absorbtivity in data2
!$omp workshare
                data2(iil:iil+nf,1) = etmp(if0l:if1l)
!$omp end workshare

                ! For each line frequency
!$omp do
                do ifreq=if0l,if1l

                  ! Left limit
                  if (ifreq.eq.if0l) then

                    jjl = iil
                    daux = daux + data2(jjl,1)*Atom(ia)%W0(jtran)

                  ! Not left limit
                  else

                    ! Get index
                    jjl = iil + ifreq - if0l

                    ! Right limit
                    if (ifreq.eq.if1l) then

                      daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)

                    ! No limit
                    else

                      daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)

                    end if ! Right limit
                  end if ! Left limit

                end do ! Frequencies
!$omp end do
!$omp end parallel
                if (daux.gt.0d0) then
                  daux = 1d0/daux
!$omp parallel workshare default(none) shared(data2,iil,nf,daux)
                  data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux
!$omp end parallel workshare
                end if ! Non zero normalization
              end if ! Not MPI
            end if ! Type of profile


            !
            ! Stimulated emission contribution
            !
!$omp parallel default(none) &
!$omp shared(etmp,if0l,if1l,iil,nf,estmp,stm,absK,pop,pE,rhou) &
!$omp shared(data1,data2,rL)
            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
!$omp workshare
              etmp(if0l:if1l) = etmp(if0l:if1l) - &
                                estmp(if0l:if1l)/absK

              !
              ! Store emission profile
              !

              ! For intensity it is the same than absorption
              data2(iil:iil+nf,2) = data2(iil:iil+nf,1)
!$omp end workshare
             !data2(iil:iil+nf,2) = estmp(if0l:if1l)

             !! If there is no MPI, normalize it here
             !if (.not.MPID%mpi) then
             !  jjl = iil
             !  daux = data2(jjl,2)*Atom(ia)%W0(jtran)
             !  do ifreq=if0l+1,if1l-1
             !    jjl = jjl + 1
             !    daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)
             !  end do
             !  if (if1l.gt.if0l) then
             !    jjl = jjl + 1
             !    daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)
             !  end if

             !  if (daux.gt.0d0) then
             !    daux = 1d0/daux
             !    data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux
             !  end if
             !end if
            endif ! Stimulated emission

!$omp workshare
            ! Add the contribution to the absorptivity of this atom
            data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                                 etmp(if0l:if1l)*pop

            ! Store the numerator of the opacity fraction of this
            ! line
            rL(iil:iil+nf) = estmp(if0l:if1l)*pE*pop/rhou
!$omp end workshare
!$omp end parallel

            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Real index in trano
              ffltran = Atom(ia)%itrano(ffjtran)

              ! Frec index
              indxf = Frec%indx(ffjtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (IRAM) then
                indx = Red%indx(ffjtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ggf0 = Frec%dzao(indxf)%ggf0
                ggf1 = Frec%dzao(indxf)%ggf1

                t0 = Atom(ia)%tfshift + 1
                t1 = t0 + Atom(ia)%nftran - 1

                ! Angle-average
                if (AVI) then
                  call emissI2ord_AA(Atom(ia),Geom, &
                                     Atmo%vx(iz),Atmo%vy(iz), &
                                     Atmo%vz(iz),Frec%omega,p_red, &
                                     p_frec, &
                                     p_Norm,jtran,fjtran,itermu, &
                                     iterml,iJu,iJl,iz, &
                                     if0l2,if1l2,DwT,Dw,vfac, &
                                     Atmo%vmi(iz),Stokes,JKQ(t0:t1), &
                                     JKQC(ggf0:ggf1), &
                                     prof(if0l2:if1l2), &
                                     es2tmp(if0l2:if1l2), &
                                     rpf(if0l2:if1l2))
                ! Angle-dependent
                else
                  call emissI2ord_AD(Atom(ia),Geom,.False., &
                                     Atmo%vx(iz),Atmo%vy(iz), &
                                     Atmo%vz(iz),Frec%omega,p_red, &
                                     p_frec, &
                                     p_Norm,jdir,jtran,fjtran, &
                                     itermu,iterml,iJu,iJl,iz, &
                                     if0l2,if1l2,DwT,Dw,vfac, &
                                     Atmo%vmi(iz),Stokes,JKQ(t0:t1), &
                                     prof(if0l2:if1l2), &
                                     es2tmp(if0l2:if1l2), &
                                     rpf(if0l2:if1l2))
                end if


                !
                ! Total emissivity
                !
!$omp parallel default(none) &
!$omp private(jjl,nf2) &
!$omp shared(estmp,if0l2,if1l2,es2tmp,iil,if0l,rL,rpf)
!$omp workshare
                estmp(if0l2:if1l2) = es2tmp(if0l2:if1l2) + &
                                     estmp(if0l2:if1l2)
!$omp end workshare

                ! Introduce the PRD factor
                jjl = iil + if0l2 - if0l
                nf2 = if1l2 - if0l2
!$omp workshare
                rL(jjl:jjl+nf2) = rL(jjl:jjl+nf2)*rpf(if0l2:if1l2)
!$omp end workshare
!$omp end parallel

              end if
            end if


            ! Add the contribution to the emissivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(data1,if0l,if1l,estmp,pE,pop)
            data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                                 estmp(if0l:if1l)*(pE*pop)
!$omp end parallel workshare

            ! Advance the index
            iil = iil + nf + 1

          end do ! FS transition
        end do ! transition


        !
        ! Photoionization
        !

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ktran = jtran + Atom(ia)%pshift

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1
          nf = if1l - if0l

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell


          !
          ! Emissivity
          !
          if (PIRAM) then

            call photoepsIS(Atom(ia),Frec%omega3(if0l:if1l), &
                            Frec%exu(if0l:if1l,iz), &
                            Atmo%T(iz),Atmo%ne(iz), &
                            jtran,ilevelu,iz,if0l,if1l, &
                            estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          else

            call photoepsI(Atom(ia),Frec%omega, &
                           Atmo%T(iz),Atmo%ne(iz), &
                           jtran,ilevelu,iz,if0l,if1l, &
                           estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          end if


          !
          ! Absorptivity
          !
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                         etmp(if0l:if1l))
!$omp parallel default(none) &
!$omp shared(rP,iip,nf,estmp,if0l,if1l,pop,rhou,data1,etmp,rstmp)
!$omp workshare
          ! ALI ratio
          rP(iip:iip+nf) = estmp(if0l:if1l)*pop/rhou

          ! Add contribution to emissivity
          data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                               estmp(if0l:if1l)*pop

          ! Remove the stimulated part
          etmp(if0l:if1l) = etmp(if0l:if1l) - rstmp(if0l:if1l)
!$omp end workshare
!$omp workshare
          ! Add contribution to absorptivity
          data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                               etmp(if0l:if1l)*pop
!$omp end workshare
!$omp end parallel
          ! Update the index
          iip = iip + nf + 1

        end do ! b-f transitions
      end do ! Atom


      !
      ! Finish the construction of the opacity fractions
      !
!$omp parallel default(none) &
!$omp private(ia,jtran,if0l,if1l,nf,fjtran,iil,iip,ktran) &
!$omp shared(nA,Atom,rL,data1,rP)
      ! Initialize
      iil = 1
      iip = 1

      ! For each atom
      do ia=1,nA

        ! For each FS b-b transition
        do jtran=1,Atom(ia)%ntran

          ! If this CPU does not have frequencies in this line, skip
          if (Atom(ia)%fflag(jtran)%absent) cycle

          ! Store frequency limits
          if0l = Atom(ia)%if0(jtran)
          if1l = Atom(ia)%if1(jtran)
          nf = if1l - if0l

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Get the opacity fraction
!$omp workshare
            rL(iil:iil+nf) = rL(iil:iil+nf)/ &
                             (data1(if0l:if1l,0) + vacuum)
!$omp end workshare

            ! Update index
            iil = iil + nf + 1

          end do ! b-b FS transitions
        end do ! b-b transitions

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! If this CPU does not have frequencies in this transition,
          ! skip
          if (Atom(ia)%phot(jtran)%absent) cycle

          ktran = jtran + Atom(ia)%pshift

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1
          nf = if1l - if0l

          ! Get the opacity fraction
!$omp workshare
          rP(iip:iip+nf) = rP(iip:iip+nf)/ &
                           (data1(if0l:if1l,0) + vacuum)
!$omp end workshare

          ! Update index
          iip = iip + nf + 1

        end do ! b-f transitions
      end do ! Atoms


      !
      ! Transform into the data arrays
      !

      ! Source function
!$omp workshare
      data1(:,1) = data1(:,1)/(data1(:,0) + vacuum)
!$omp end workshare
!$omp end parallel

      ! Nullify pointers
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)


      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTCoeffI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      close(200)
      urou = 'RTCoeffI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTCoeffI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the radiation transfer coefficients for emergencei
      !! or emergencei_serial\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!            JKQ(dfloat(:)): Mean intensity integrated over
      !!                            absorption profile\n
      !!         JKQC(dfloat(:,:)): Mean intensity with frequency
      !!                            dependence\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!     Stokes(dfloat(:,:,:)): Intensity\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients
      subroutine RTCoeffIe(Frec,Red,Atom,LTElines,Atmo,Geom,iz, &
                           ith,iph,if0,if1,JKQ,JKQC,cdir,Cont, &
                           Stokes,data1)

      ! I/O

      type(Frequency_class), intent(inout):: Frec
      type(Red_class), intent(inout):: Red
      type(Geometry_class), intent(in):: Geom
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh), &
                        intent(in):: Stokes
      double precision, dimension(nxt), intent(in):: JKQ
      double precision, dimension(nfreq), intent(in):: JKQC
      double precision, dimension(if0:if1,0:1):: data1

      ! Local

      integer:: iterml,itermu,iJl,iJu,ilevell,ilevelu,indx,indxf
      integer:: ia,jtran,fjtran,ffjtran,ffktran,ggf0,ggf1
      integer:: jdir,jbdir,jcdir,jddir,icdir,ios,nodir
      integer:: if0l,if1l,if0l2,if1l2,t0,t1,offset

      double precision:: DwT,Dfreq,pE,absK,Dw,vfac,pop
      double precision:: ct,st,cc,sc,rhou,loffset
      double precision, dimension(if0:if1):: rpf,etmp,prof, &
                                             estmp,rstmp,es2tmp

      type(Frequencyc2_class), pointer:: p_frec
      type(Redc2_class), pointer:: p_red
      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(p_frec,p_red,p_Norm,p_LTEprof)


      !
      ! Get the direction index for the continuum and calculate
      ! the Doppler shift if necessary
      !
      jdir = Geom%i_geom(iph,ith)
      if (PRD) then
        jbdir = min(jdir,Frec%ndir)
        jcdir = min(jdir,Red%ndir)
        jddir = min(jdir,Red%njdir)
      end if

      ! If there are dynamics
      if (dyn) then

        icdir = min(jdir,cdir)

        nodir = jdir

        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If there are not dynamics
      else

        icdir = 1

        nodir = 1

        vfac = 1d0

      end if

      ! If using file, nullify pointer
      if (VIFIL) p_Norm => Atom(1)%Normp(1,1,1)


      !
      ! Continuum contribution
      !

      ! Absorptivity
      data1(:,0) = Cont(:,1,icdir)

      ! Emissivity
      data1(:,1) = Cont(:,3,icdir) + Cont(:,2,icdir)*JKQC(if0:if1)

      !
      ! Calculate RT coefficients
      !

      ! LTE lines
      do ia=1,nLTEl

        ! Skip if absent
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = LTElines(ia)%n(iz)

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Get Norm
        if (LVIRAM) then
          p_LTEprof => LTElines(ia)%prof(iz,nodir)
        else
          p_LTEprof => LTElines(1)%prof(1,1)
        end if

        ! Get frequency of FS transition
        Dfreq = LTElines(ia)%Eu - LTElines(ia)%El

        ! Photon energy (cgs) and convertion factor
        pE = convF*Dfreq
        absK = 1d21*(2d0*c)*Dfreq**2d0

        ! Add the microturbulence to Doppler width
        Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        ! RT coeffs
        call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                        p_LTEprof,Dw,vfac,pE, &
                        etmp(if0l:if1l),estmp(if0l:if1l))


        !
        ! Stimulated emission contribution
        !
!$omp parallel default(none) &
!$omp shared(stm,etmp,if0l,if1l,estmp,absK,data1,pop)
        ! If there is stimulated emission
        if (stm) then
!$omp workshare
          etmp(if0l:if1l) = etmp(if0l:if1l) - &
                            estmp(if0l:if1l)/absK
!$omp end workshare
        end if

        ! Add the contribution to the absorptivity of this atom
!$omp workshare
        data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                             etmp(if0l:if1l)*pop
!$omp end workshare
!$omp end parallel

        ! Add the contribution to the emissivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(estmp,data1,if0l,if1l,pE,pop)
        data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                             estmp(if0l:if1l)*(pE*pop)
!$omp end parallel workshare

      end do ! LTE lines


      !
      ! For each atom
      do ia=1,nA

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = Atom(ia)%n(iz)


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

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          iterml = Atom(ia)%fst(jtran)%iterml

          ! Get Norm
          if (.not.VIFIL) p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJl,iterml)

            ! Photon energy (cgs) and convertion factor
            pE = convF*Dfreq
            absK = 1d21*(2d0*c)*Dfreq**2d0

            ! Add the microturbulence to Doppler width
            Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


            !
            ! If reading from file
            !

            if (vifil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(ffjtran) + &
                        Atom(ia)%f0size(ffjtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              read(200, err=1100) prof(if0l:if1l)

              close(200)

            end if


            ! RT coeffs
            call rt1ordI(Atom(ia),Frec%omega, &
                         jtran,fjtran,itermu,iterml,iJu,iJl,iz, &
                         if0l,if1l,p_Norm,Dw,vfac,pE, &
                         prof(if0l:if1l),etmp(if0l:if1l), &
                         estmp(if0l:if1l),rhou)


            !
            ! Stimulated emission contribution
            !

!$omp parallel default(none) &
!$omp shared(stm,etmp,if0l,if1l,estmp,absK,data1,pop)
            ! If there is stimulated emission
            if (stm) then
!$omp workshare
              etmp(if0l:if1l) = etmp(if0l:if1l) - &
                                estmp(if0l:if1l)/absK
!$omp end workshare
            end if

            ! Add the contribution to the absorptivity of this atom
!$omp workshare
            data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                                 etmp(if0l:if1l)*pop
!$omp end workshare
!$omp end parallel


            !
            ! Second order emissivity
            !
            if(PRD.and.Atom(ia)%lemiss2(jtran))then

              ! Real index
              ffktran = Atom(ia)%itrano(ffjtran)

              ! Frec index
              indxf = Frec%indx(ffjtran,ia,iz,jbdir)
              p_frec => Frec%dzao(indxf)

              ! If storing
              if (IRAM) then
                indx = Red%indx(ffjtran,ia,iz,jddir)
                p_red => Red%dzao(indx)
              else
                p_red => Red%dzao(1)
              end if

              if (Frec%dzao(indxf)%nran.gt.0) then

                if0l2 = Frec%dzao(indxf)%gf0
                if1l2 = Frec%dzao(indxf)%gf1

                ggf0 = Frec%dzao(indxf)%ggf0
                ggf1 = Frec%dzao(indxf)%ggf1

                t0 = Atom(ia)%tfshift + 1
                t1 = t0 + Atom(ia)%nftran - 1

                ! Angle-average
                if (AVI) then
                  call emissI2ord_AA(Atom(ia),Geom, &
                                     Atmo%vx(iz),Atmo%vy(iz), &
                                     Atmo%vz(iz),Frec%omega,p_red, &
                                     p_frec, &
                                     p_Norm,jtran, &
                                     fjtran,itermu,iterml,iJu,iJl, &
                                     iz,if0l2,if1l2,DwT, &
                                     Dw,vfac,Atmo%vmi(iz),Stokes, &
                                     JKQ(t0:t1),JKQC(ggf0:ggf1), &
                                     prof(if0l2:if1l2), &
                                     es2tmp(if0l2:if1l2), &
                                     rpf(if0l2:if1l2))
                ! Angle-dependent
                else
                  call emissI2ord_AD(Atom(ia),Geom,.True., &
                                     Atmo%vx(iz),Atmo%vy(iz), &
                                     Atmo%vz(iz),Frec%omega,p_red, &
                                     p_frec, &
                                     p_Norm,1,jtran, &
                                     fjtran,itermu,iterml,iJu,iJl, &
                                     iz,if0l2,if1l2,DwT, &
                                     Dw,vfac,Atmo%vmi(iz),Stokes, &
                                     JKQ(t0:t1), &
                                     prof(if0l2:if1l2), &
                                     es2tmp(if0l2:if1l2), &
                                     rpf(if0l2:if1l2))
                end if


                !
                ! Total emissivity
                !
!$omp parallel workshare default(none) &
!$omp shared(estmp,es2tmp,if0l2,if1l2)
                estmp(if0l2:if1l2) = es2tmp(if0l2:if1l2) + &
                                      estmp(if0l2:if1l2)
!$omp end parallel workshare

              end if
            end if

            ! Add the contribution to the emissivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(estmp,data1,if0l,if1l,pE,pop)
            data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                                 estmp(if0l:if1l)*(pE*pop)
!$omp end parallel workshare

          end do ! FS transition
        end do ! transition


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

            call photoepsIS(Atom(ia),Frec%omega3(if0l:if1l), &
                            Frec%exu(if0l:if1l,iz), &
                            Atmo%T(iz), &
                            Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                            estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          else

            call photoepsI(Atom(ia),Frec%omega, &
                           Atmo%T(iz),Atmo%ne(iz),jtran,ilevelu, &
                           iz,if0l,if1l,estmp(if0l:if1l), &
                           rstmp(if0l:if1l),rhou)

          end if


          !
          ! Absorptivity
          !
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                         etmp(if0l:if1l))

!$omp parallel default(none) &
!$omp shared(data1,if0l,if1l,estmp,etmp,pop,rstmp)
          ! Add contribution to emissivity
!$omp workshare
          data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                               estmp(if0l:if1l)*pop

          ! Remove the stimulated part
          etmp(if0l:if1l) = etmp(if0l:if1l) - rstmp(if0l:if1l)
!$omp end workshare

          ! Add contribution to absorptivity
!$omp workshare
          data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                               etmp(if0l:if1l)*pop
!$omp end workshare
!$omp end parallel

        end do ! b-f transitions

      end do ! Atom


      !
      ! Transform into the data arrays
      !

      ! Source function
!$omp parallel workshare default(none) shared(data1)
      data1(:,1) = data1(:,1)/(data1(:,0) + vacuum)
!$omp end parallel workshare

      ! Nullify pointers
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)

      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTCoeffIe'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      close(200)
      urou = 'RTCoeffIe'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_frec)) nullify(p_frec)
      if (associated(p_red)) nullify(p_red)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTCoeffIe

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the intensity absorption coefficients\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!   LTElines(LTEline_class): Structure with the LTE line data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
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
      !!           etaI(dfloat(:)): Absorption coefficient
      subroutine RTAbsI(Frec,Atom,LTElines,Atmo,Geom,iz, &
                        ith,iph,if0,if1,cdir,Cont,etaI)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: Geom
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(LTEline_class), dimension(:), intent(in), allocatable:: &
                                                              LTElines
      integer, intent(in):: iz,ith,iph,cdir,if0,if1
      double precision, &
                  dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(if0:if1):: etaI

      ! Local

      integer:: iterml,itermu,iJl,iJu,ilevell,ilevelu,offset,ios
      integer:: ia,jtran,fjtran,ffjtran,jdir,icdir,if0l,if1l,nodir

      double precision:: DwT,Dfreq,pE,absK,Dw,vfac,pop
      double precision:: ct,st,cc,sc,rhou,loffset
      double precision, dimension(if0:if1):: prof,etmp,estmp,rstmp

      type(Nindex_class), pointer:: p_Norm
      type(LTEprof_class), pointer:: p_LTEprof

      !
      ! Initialize
      !
      nullify(p_Norm,p_LTEprof)


      !
      ! Get the direction index for the continuum and calculate
      ! the Doppler shift if necessary
      !
      jdir = Geom%i_geom(iph,ith)

      ! If there are dynamics
      if (dyn) then

        icdir = min(jdir,cdir)

        nodir = jdir

        ct = Geom%L_mu(ith)
        st = sqrt(1d0 - ct*ct)
        cc = cos(Geom%L_phi(iph))
        sc = sin(Geom%L_phi(iph))

        vfac = 1d0 - atmo%vx(iz)*st*cc - atmo%vy(iz)*st*sc - &
                     atmo%vz(iz)*ct

      ! If there are not dynamics
      else

        icdir = 1

        nodir = 1

        vfac = 1d0

      end if

      ! If using file, nullify pointer
      if (VIFIL) p_Norm => Atom(1)%Normp(1,1,1)


      !
      ! Continuum contribution
      !

      etaI = Cont(:,1,icdir)


      !
      ! Calculate RT coefficients
      !

      ! LTE lines
      do ia=1,nLTEl

        ! Skip if absent
        if (LTElines(ia)%absent) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = LTElines(ia)%n(iz)

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Get Norm
        if (LVIRAM) then
          p_LTEprof => LTElines(ia)%prof(iz,nodir)
        else
          p_LTEprof => LTElines(1)%prof(1,1)
        end if

        ! Get frequency of FS transition
        Dfreq = LTElines(ia)%Eu - LTElines(ia)%El

        ! Photon energy (cgs) and convertion factor
        pE = convF*Dfreq
        absK = 1d21*(2d0*c)*Dfreq**2d0

        ! Add the microturbulence to Doppler width
        Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        !
        ! If stimulated emission
        if (stm) then

          ! RT coeffs
          call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_LTEprof,Dw,vfac,pE, &
                          etmp(if0l:if1l),estmp(if0l:if1l))

          !
          ! Stimulated emission contribution
          !
!$omp parallel workshare default(none) &
!$omp shared(etmp,if0l,if1l,estmp,absK)
          etmp(if0l:if1l) = etmp(if0l:if1l) - &
                            estmp(if0l:if1l)/absK
!$omp end parallel workshare

        ! No stimulated emission
        else

          !
          ! Absorptivity
          !
          call absorbILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_LTEprof,Dw,vfac,pE,etmp(if0l:if1l))

        end if

        ! Add the contribution to the absorptivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(etmp,if0l,if1l,etaI,pop)
        etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)
!$omp end parallel workshare

      end do ! LTE lines


      !
      ! For each atom
      do ia=1,nA

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store population in double precision
        pop = Atom(ia)%n(iz)


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

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          iterml = Atom(ia)%fst(jtran)%iterml

          ! Get Norm
          if (.not.VIFIL) p_Norm => Atom(ia)%Normp(jtran,iz,nodir)

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJl,iterml)

            ! Photon energy (cgs) and convertion factor
            pE = convF*Dfreq
            absK = 1d21*(2d0*c)*Dfreq**2d0

            ! Add the microturbulence to Doppler width
            Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


            !
            ! If reading from file
            !

            if (vifil) then

              ! Open files
              open(200, file=trim(Atom(ia)%vfile), status='unknown', &
                   iostat=ios, err=1000, access='stream', &
                   action='read', form='unformatted')

              ! Jump
              loffset = dble(Atom(ia)%hvifil) + &
                        Atom(ia)%dsize(nodir) + &
                        Atom(ia)%zsize(iz) + &
                        Atom(ia)%tsize(ffjtran) + &
                        Atom(ia)%f0size(ffjtran)
              do while(loffset.gt.offlimit)
                offset = int(offlimit)
                call fseek(200, offset, 1)
                loffset = loffset - offlimit
              end do
              offset = int(loffset)
              call fseek(200, offset, 1)

              read(200, err=1100) prof(if0l:if1l)

              close(200)

            end if

            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! RT coeffs
              call rt1ordI(Atom(ia),Frec%omega, &
                           jtran,fjtran,itermu,iterml,iJu,iJl,iz, &
                           if0l,if1l,p_Norm,Dw,vfac,pE, &
                           prof(if0l:if1l),etmp(if0l:if1l), &
                           estmp(if0l:if1l),rhou)

!$omp parallel workshare default(none) &
!$omp shared(etmp,if0l,if1l,estmp,absK)
              etmp(if0l:if1l) = etmp(if0l:if1l) - &
                                estmp(if0l:if1l)/absK
!$omp end parallel workshare

            ! No stimulated emission
            else

              ! Absorptivity
              call absorbI(Atom(ia),Frec%omega, &
                           jtran,fjtran,itermu,iterml,iJu,iJl,iz, &
                           if0l,if1l,p_Norm,Dw,vfac,pE, &
                           prof(if0l:if1l),etmp(if0l:if1l))

            end if

            ! Add the contribution to the absorptivity of this atom
!$omp parallel workshare default(none) &
!$omp shared(etmp,if0l,if1l,etaI,pop)
            etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)
!$omp end parallel workshare

          end do ! FS transition
        end do ! transition


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

            call photoepsIS(Atom(ia),Frec%omega3(if0l:if1l), &
                            Frec%exu(if0l:if1l,iz), &
                            Atmo%T(iz),Atmo%ne(iz), &
                            jtran,ilevelu,iz,if0l,if1l, &
                            estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          else

            call photoepsI(Atom(ia),Frec%omega,Atmo%T(iz), &
                           Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                           estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          end if

          !
          ! Absorptivity
          !
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                         etmp(if0l:if1l))

!$omp parallel default(none) &
!$omp shared(etmp,rstmp,etaI,if0l,if1l,pop)
          ! Remove the stimulated part
!$omp workshare
          etmp(if0l:if1l) = etmp(if0l:if1l) - rstmp(if0l:if1l)
!$omp end workshare

          ! Add contribution to absorptivity
!$omp workshare
          etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)
!$omp end workshare
!$omp end parallel

        end do ! b-f transitions

      end do ! Atom

      ! Nullify pointers
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)

      return

1000  umsg = 'Error opening '//trim(Atom(ia)%vfile)//' file'
      urou = 'RTAbsI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      close(200)
      urou = 'RTAbsI'
      call abortedS(umsg,urou,-1,.True.,.True.)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine RTAbsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the continuum radiation transfer coefficients\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!           JKQC(dfloat(:)): Mean intensity with frequency
      !!                            dependence\n
      !!             cdir(integer): Direction index for background
      !!                            opacities\n
      !!       Cont(dfloat(:,:,:)): Background opacity data\n
      !!        data1(dfloat(:,:)): Radiation transfer coefficients
      subroutine RTCoeffJ(Geom,ith,iph,if0,if1,JKQC,cdir,Cont, &
                          data1)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      integer, intent(in):: ith,iph,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir):: Cont
      double precision, dimension(nfreq), intent(in):: JKQC
      double precision, dimension(if0:if1,0:1):: data1

      ! Local

      integer:: jdir,icdir


      !
      ! Get the direction index for the continuum and calculate
      ! the Doppler shift if necessary
      !
      jdir = Geom%i_geom(iph,ith)

      ! Index of continuum direction
      icdir = 1

      ! If there are dynamics
      if (dyn) icdir = min(jdir,cdir)

      !
      ! Calculate RT coefficients
      !

      ! Absorptivity
!$omp parallel default(none) &
!$omp shared(data1,Cont,icdir,JKQC,if0,if1)
!$omp workshare
      data1(:,0) = Cont(:,1,icdir)

      ! Emissivity
      data1(:,1) = Cont(:,3,icdir) + Cont(:,2,icdir)*JKQC(if0:if1)
!$omp end workshare

      ! Source function
!$omp workshare
      data1(:,1) = data1(:,1)/(data1(:,0) + vacuum)
!$omp end workshare
!$omp end parallel

      return

      end subroutine RTCoeffJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes the line profiles for JKQgen or JKQgen_serial\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!        Flgsg(Fctsg_class): Structure with factorials and
      !!                            signs\n
      !!      Geom(Geometry_class): Structure with geometry data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!               iz(integer): Height index\n
      !!              ith(integer): Output direction polar index\n
      !!              iph(integer): Output direction azimuth index\n
      !!              if0(integer): First frequency index for this
      !!                            CPU\n
      !!              if1(integer): Last frequency index for this
      !!                            CPU\n
      !!        data2(dfloat(:,:)): Line profiles
      subroutine Termprof(Frec,Atom,Atmo,MPID,Flgsg,Geom, &
                          Bfield,iz,ith,iph,if0,if1,data2)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), intent(in):: Atom
      type(Frequency_class), intent(inout):: Frec
      type(MPI_class), intent(in):: MPID
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,ith,iph,if0,if1
      double precision, dimension(:,:):: data2

      ! Local

      integer:: iterml,itermu,ia,jtran,ktran,jdir,icom,ios
      integer:: ifreq,if0l,if1l,iil,jjl,nf,offset,nodir

      double precision:: daux,DwT,absK,Dw,vfac,ct,st,cc,sc,loffset
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
      nullify(p_Norm,prof,p_prof,p_LTEprof)


      ! Formal solution direction
      jdir = Geom%i_geom(iph,ith)
      TBo = Geom%TB(:,:,:,iph,ith,iz)
      TSo = Geom%TS(:,:,:,iph,ith)

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


      ! If using file, nullify pointer to norm and allocate profile
      if (VPFIL) then

        p_Norm => Atom(1)%Normp(1,1,1)
        allocate(prof(if0:if1,Atom(1)%Mncom))

      ! If not using file
      else

        allocate(p_prof(1,1))

      end if


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

      ! Initialize index
      iil = 1

      ! For each atom
      do ia=1,nA

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
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

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

            ! If MPI
            if (MPID%mpi) then

              ! Copy absorbtivity in data2
!$omp parallel workshare default(none) &
!$omp shared(data2,etmp0,iil,nf,if0l,if1l)
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)
!$omp end parallel workshare

            ! No MPI
            else

              ! Initialize
              daux = 0d0

!$omp parallel default(none) &
!$omp private(jjl,ifreq) &
!$omp shared(data2,etmp0,Atom,ia,jtran,iil,if0l,if1l,Frec,nf) &
!$omp reduction(+: daux)

              ! Store absorption profile
!$omp workshare
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)
!$omp end workshare

              ! For each line frequency
!$omp do
              do ifreq=if0l,if1l

                ! Left limit
                if (ifreq.eq.if0l) then

                  jjl = iil
                  daux = daux + data2(jjl,1)*Atom(ia)%W0(jtran)

                ! Not left limit
                else

                  ! Get index
                  jjl = iil + ifreq - if0l

                  ! Right limit
                  if (ifreq.eq.if1l) then

                    daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)

                  ! No limit
                  else

                    daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)

                  end if ! Right limit
                end if ! Left limit

              end do ! Frequencies
!$omp end do
!$omp end parallel
              if (daux.gt.0d0) then
                daux = 1d0/daux
!$omp parallel workshare default(none) shared(data2,iil,nf,daux)
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux
!$omp end parallel workshare
              end if
            end if


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! If MPI
              if (MPID%mpi) then
!$omp parallel workshare default(none) &
!$omp shared(etmp0,if0l,if1l,iil,nf,estmp0,absK,data2)
                ! Correct for stimulated emission
                etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                   estmp0(if0l:if1l)/absK

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)
!$omp end parallel workshare
              ! No MPI
              else

                ! Initialize
                daux = 0d0

!$omp parallel default(none) &
!$omp private(ifreq,jjl) &
!$omp shared(etmp0,if0l,if1l,estmp0,data2,nf,iil,Frec,Atom,ia,jtran) &
!$omp shared(absK) &
!$omp reduction(+: daux)
!$omp workshare
                ! Correct for stimulated emission
                etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                   estmp0(if0l:if1l)/absK

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)
!$omp end workshare

                ! For each line frequency
!$omp do
                do ifreq=if0l,if1l

                  ! Left limit
                  if (ifreq.eq.if0l) then

                    jjl = iil
                    daux = daux + data2(jjl,2)*Atom(ia)%W0(jtran)

                  ! Not left limit
                  else

                    ! Get index
                    jjl = iil + ifreq - if0l

                    ! Right limit
                    if (ifreq.eq.if1l) then

                      daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)

                    ! No limit
                    else

                      daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)

                    end if ! Right limit
                  end if ! Left limit

                end do ! Frequencies
!$omp end do
!$omp end parallel
                if (daux.gt.0d0) then
                  daux = 1d0/daux
!$omp parallel workshare default(none) shared(data2,iil,nf,daux)
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux
!$omp end parallel workshare
                end if
              end if
            endif ! Stimulated emission

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        ! No field
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
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)

            !
            ! First order emissivity
            !
            call emissNB(Atom(ia),TBo,Frec%omega,Flgsg, &
                         jtran,itermu,iterml,iz,if0l,if1l, &
                         p_Norm,Dw,vfac,p_prof, &
                         estmp0(if0l:if1l),estmp1(if0l:if1l), &
                         estmp2(if0l:if1l),estmp3(if0l:if1l), &
                         rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                         rstmp3(if0l:if1l))

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

            ! If MPI
            if (MPID%mpi) then

              ! Copy absorbtivity in data2
!$omp parallel workshare default(none) &
!$omp shared(data2,etmp0,iil,nf,if0l,if1l)
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)
!$omp end parallel workshare

            ! No MPI
            else

              ! Initialize
              daux = 0d0

!$omp parallel default(none) &
!$omp private(jjl,ifreq) &
!$omp shared(data2,etmp0,Atom,ia,jtran,iil,if0l,if1l,Frec,nf) &
!$omp reduction(+: daux)

              ! Store absorption profile
!$omp workshare
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)
!$omp end workshare

              ! For each line frequency
!$omp do
              do ifreq=if0l,if1l

                ! Left limit
                if (ifreq.eq.if0l) then

                  jjl = iil
                  daux = daux + data2(jjl,1)*Atom(ia)%W0(jtran)

                ! Not left limit
                else

                  ! Get index
                  jjl = iil + ifreq - if0l

                  ! Right limit
                  if (ifreq.eq.if1l) then

                    daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)

                  ! No limit
                  else

                    daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)

                  end if ! Right limit
                end if ! Left limit

              end do ! Frequencies
!$omp end do
!$omp end parallel
              if (daux.gt.0d0) then
                daux = 1d0/daux
!$omp parallel workshare default(none) shared(data2,iil,nf,daux)
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux
!$omp end parallel workshare
              end if
            end if


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! If MPI
              if (MPID%mpi) then
!$omp parallel workshare default(none) &
!$omp shared(etmp0,if0l,if1l,iil,nf,estmp0,absK,data2)
                ! Correct for stimulated emission
                etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                   estmp0(if0l:if1l)/absK

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)
!$omp end parallel workshare
              ! No MPI
              else

                ! Initialize
                daux = 0d0

!$omp parallel default(none) &
!$omp private(ifreq,jjl) &
!$omp shared(etmp0,if0l,if1l,estmp0,data2,iil,nf,Frec,Atom,ia,jtran) &
!$omp shared(absK) &
!$omp reduction(+: daux)
!$omp workshare
                ! Correct for stimulated emission
                etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                   estmp0(if0l:if1l)/absK

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)
!$omp end workshare

                ! For each line frequency
!$omp do
                do ifreq=if0l,if1l

                  ! Left limit
                  if (ifreq.eq.if0l) then

                    jjl = iil
                    daux = daux + data2(jjl,2)*Atom(ia)%W0(jtran)

                  ! Not left limit
                  else

                    ! Get index
                    jjl = iil + ifreq - if0l

                    ! Right limit
                    if (ifreq.eq.if1l) then

                      daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)

                    ! No limit
                    else

                      daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)

                    end if ! Right limit
                  end if ! Left limit

                end do ! Frequencies
!$omp end do
!$omp end parallel
                if (daux.gt.0d0) then
                  daux = 1d0/daux
!$omp parallel workshare default(none) shared(data2,iil,nf,daux)
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux
!$omp end parallel workshare
                end if
              end if
            endif ! Stimulated emission

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        end if ! Magnetic field

      end do ! Atoms

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
      urou = 'Termprof'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return
1100  umsg = 'Error reading '//trim(Atom(ia)%vfile)//' file'
      close(200)
      urou = 'Termprof'
      call abortedS(umsg,urou,-1,.True.,.True.)
      deallocate(prof)
      nullify(prof)
      if (associated(p_prof)) nullify(p_prof)
      if (associated(p_norm)) nullify(p_norm)
      if (associated(p_LTEprof)) nullify(p_LTEprof)
      return

      end subroutine Termprof

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffi_mod
