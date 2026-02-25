      !> Compute radiation transfer coefficients for the NLTE problem
      !! of the second kind
      module rtcoeff_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     27/04/2017
!  Last version:
!     20/02/2026 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     20/02/2026:    V4.0.3 - Bugfix: the render limit of LTE lines
!                             was not being accounted for (TdPA)
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
!    - Second order emissivity for the CLE case
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  RTCoeff
!    Calculate the radiation transfer coefficients for a given height
!  and direction for the polarized problem
!
!  RTAbs
!    Calculate the opacity for a given height and direction for the
!  polarized problem
!
!  RTCoeff_CLE
!    Calculate the RT coefficients at a given point in a CLE
!  synthesis
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
      use rtcoeffiaux_mod , only : absorbILTE , rt1ordILTE , &
                                   photoabsI
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the radiation transfer coefficients for a given
      !! height and direction for the polarized problem\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!        Geom(Geometry_class): Structure with geometric data\n
      !!                vfac(double): Doppler shift factor\n
      !!                 iz(integer): Height index\n
      !!               jdir(integer): Current direction index\n
      !!                if0(integer): Lower limit index for
      !!                              frequency\n
      !!                if1(integer): Upper limit index for
      !!                              frequency\n
      !!       JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                              frequency dependence\n
      !!               cdir(integer): Direction index for background
      !!                              opacities\n
      !!         Cont(double(:,:,:)): Background opacity data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!         TSo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                              vertical reference frame\n
      !!        TKQo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                              suitable reference frame\n
      !!          data1(double(:,:)): Radiation transfer
      !!                              coefficients\n
      !!          data2(double(:,:)): Line profiles\n
      !!          iterating(logical): If solving self-consistent
      !!                              problem
      subroutine RTCoeff(Frec,Red,Atom,LTElines,Atmo,Flgsg, &
                         Geom,vfac,iz,jdir,if0,if1,JKQC, &
                         cdir,Cont,Bfield,TSo,TKQo, &
                         data1,data2,iterating)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      logical, intent(in):: iterating
      integer, intent(in):: iz,jdir,cdir,if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(if0:if1,3,cdir), intent(in):: Cont
      double precision, dimension(0:3,if0:if1,0:4), &
                        intent(out):: data1
      double precision, dimension(:,:), intent(inout):: data2
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(in):: JKQC
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TSo
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQo

      ! Local

      integer:: iS,K,iQ,iterml,itermu,ia,jtran,ktran
      integer:: icdir,ilevell,ilevelu,nodir,indx,jndx
      integer:: ifreq,if0l,if1l,if0l2,if1l2,iil,jjl,nf

      double precision:: daux,DwT,pE,absK,Dw
      double precision, dimension(if0:if1,0:3):: etaA,epsA,rhsA,rhaA
      double precision, dimension(if0:if1):: etmp0,estmp0
      double precision, dimension(if0:if1):: etmp1,estmp1
      double precision, dimension(if0:if1):: rstmp1,rtmp1
      double precision, dimension(if0:if1):: etmp2,estmp2
      double precision, dimension(if0:if1):: rstmp2,rtmp2
      double precision, dimension(if0:if1):: etmp3,estmp3
      double precision, dimension(if0:if1):: rstmp3,rtmp3
      double precision, dimension(if0:if1):: intgr

      ! Pointer

      type(Prof_class), pointer:: p_Norm


      ! Initialize
      nullify(p_Norm)

      ! Initialize continuum index
      icdir = 1

      ! Correct direction index
      nodir = min(jdir,ubound(Red%idzao,3),Geom%njdir)

      ! If there are dynamics, correct continuum index
      if (dyn) icdir = min(jdir,cdir)

      ! Initialize to 0 the parts without continuum
      data1(1,:,0) = 0d0

      ! If not axial
      if (.not.RTaxial) then

        ! Initialize to 0 the parts without continuum
        data1(2,:,0) = 0d0
        data1(3,:,0) = 0d0
        data1(2,:,1) = 0d0
        data1(3,:,1) = 0d0
        data1(3,:,2) = 0d0

      end if ! Non-axial RT


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
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If this CPU does not have frequencies in this line, skip
        if (LTElines(ia)%absent) cycle

        ! Skip if limited
        if (iz.lt.LTElines(ia)%Rz0) cycle

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

        ! Point to this line norm
        indx = Red%idzao(nxtran+ia,iz,nodir)
        p_Norm => Red%dzao(indx)

        !
        ! Check if magnetic field
        !

        ! If magnetic field
        if (Bfield%Bstrength(iz).gt.TINYB) then

          ! Get RT coeffs.
          call rt1ordLTE(LTElines(ia),TKQo,Frec%omega,Flgsg,iz, &
                         if0l,if1l,p_Norm,Dw,vfac, &
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

          ! Non-axial radiation transfer
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

          ! RT coeffs
          call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,pE, &
                          etmp0(if0l:if1l),estmp0(if0l:if1l))

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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Point to this line norm
            indx = Red%idzao(ktran,iz,nodir)
            p_Norm => Red%dzao(indx)

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)
            nf = if1l - if0l

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
            ! First order RT coefficients
            !
            call rt1ord(Atom(ia),TKQo,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,absK, &
                        etmp0(if0l:if1l),etmp1(if0l:if1l), &
                        etmp2(if0l:if1l),etmp3(if0l:if1l), &
                        rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                        rtmp3(if0l:if1l), &
                        estmp0(if0l:if1l),estmp1(if0l:if1l), &
                        estmp2(if0l:if1l),estmp3(if0l:if1l), &
                        rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                        rstmp3(if0l:if1l))

            ! If self-consistent solution
            if (iterating) then

              ! Store absorption profile
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)

              ! If there is no MPI, normalize it here
              if (nproc.le.1) then

                ! Get index
                jjl = iil

                ! Add to norm
                daux = data2(jjl,1)*Atom(ia)%W0(jtran)

                ! For non-boundary frequencies
                do ifreq=if0l+1,if1l-1

                  ! Get index
                  jjl = jjl + 1

                  ! Add to norm
                  daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)

                end do ! Internal frequencies

                ! If not a single point
                if (if1l.gt.if0l) then

                  ! Get index
                  jjl = jjl + 1

                  ! Add to norm
                  daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)

                end if ! Not a single point

                ! Non-zero norm
                if (daux.gt.0d0) then

                  ! Normalize profile
                  daux = 1d0/daux
                  data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux

                end if ! Non-zero norm
              end if ! MPI
            end if ! Iterating


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

              ! If self-consistent solution
              if (iterating) then

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)

                ! If there is no MPI, normalize it here
                if (nproc.le.1) then

                  ! Get index
                  jjl = iil

                  ! Add to norm
                  daux = data2(jjl,2)*Atom(ia)%W0(jtran)

                  ! For non-boundary frequencies
                  do ifreq=if0l+1,if1l-1

                    ! Get index
                    jjl = jjl + 1

                    ! Add to norm
                    daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)

                  end do ! non-boundary frequencies

                  ! If not a single frequency
                  if (if1l.gt.if0l) then

                    ! Get index
                    jjl = jjl + 1

                    ! Add to norm
                    daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)

                  end if ! Single frequency

                  ! Non-zero norm
                  if (daux.gt.0d0) then

                    ! Normalize
                    daux = 1d0/daux
                    data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux

                  end if ! Non-zero norm
                end if ! MPI
              end if ! Self-consistent solution
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
            if(PRD.and.Atom(ia)%lemiss2(jtran).and. &
               iz.lt.Rz1_PRD)then

              ! Red index
              jndx = Red%izao(jtran,ia,Rz0)
              indx = Red%izao(jtran,ia,iz)

              ! There are valid ranges
              if (Red%ao(jndx)%nran.gt.0) then

                ! Get limits for emissivity
                if0l2 = Red%ao(jndx)%gf0
                if1l2 = Red%ao(jndx)%gf1

                ! Valid range
                if (if1l2.ge.if0l2) then

                  ! If this height has only one dimension (LOS)
                  if (size(Red%zao(indx)%eps20,1).eq.1) then

                    ! Get values from corresponding direction
                    estmp0(if0l2:if1l2) = Red%zao(indx)%eps20(1,:)
                    estmp1(if0l2:if1l2) = Red%zao(indx)%eps21(1,:)
                    estmp2(if0l2:if1l2) = Red%zao(indx)%eps22(1,:)
                    estmp3(if0l2:if1l2) = Red%zao(indx)%eps23(1,:)

                  ! If this height has more than one direction
                  else

                    ! Get values from corresponding direction
                    estmp0(if0l2:if1l2) = Red%zao(indx)%eps20(jdir,:)
                    estmp1(if0l2:if1l2) = Red%zao(indx)%eps21(jdir,:)
                    estmp2(if0l2:if1l2) = Red%zao(indx)%eps22(jdir,:)
                    estmp3(if0l2:if1l2) = Red%zao(indx)%eps23(jdir,:)

                  end if ! Number of directions
                end if ! Valid range
              end if ! Valid ranges
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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Point to this line norm
            indx = Red%idzao(ktran,iz,nodir)
            p_Norm => Red%dzao(indx)

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)
            nf = if1l - if0l

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
            call rt1ordNB(Atom(ia),TKQo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l), &
                          estmp0(if0l:if1l),estmp1(if0l:if1l), &
                          estmp2(if0l:if1l),estmp3(if0l:if1l), &
                          rstmp1(if0l:if1l),rstmp2(if0l:if1l), &
                          rstmp3(if0l:if1l))

            ! If self-consistent solution
            if (iterating) then

              ! Store absorption profile
              data2(iil:iil+nf,1) = etmp0(if0l:if1l)

              ! If there is no MPI, normalize it here
              if (nproc.le.1) then

                ! Get index
                jjl = iil

                ! Add to norm
                daux = data2(jjl,1)*Atom(ia)%W0(jtran)

                ! For non-boundary frequencies
                do ifreq=if0l+1,if1l-1

                  ! Get index
                  jjl = jjl + 1

                  ! Add to norm
                  daux = daux + data2(jjl,1)*Frec%W_freq(ifreq)

                end do ! Non-boundary frequencies

                ! Not single frequency
                if (if1l.gt.if0l) then

                  ! Get index
                  jjl = jjl + 1

                  ! Add to norm
                  daux = daux + data2(jjl,1)*Atom(ia)%W1(jtran)

                end if ! single frequency

                ! If non-zero norm
                if (daux.gt.0d0) then

                  ! Normalize
                  daux = 1d0/daux
                  data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux

                end if ! Non-zero norm
              end if ! MPI
            end if ! Iterating


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

              ! If iterating
              if (iterating) then

                ! Store emission profile
                data2(iil:iil+nf,2) = estmp0(if0l:if1l)

                ! If there is no MPI, normalize it here
                if (nproc.le.1) then

                  ! Get index
                  jjl = iil

                  ! Add to norm
                  daux = data2(jjl,2)*Atom(ia)%W0(jtran)

                  ! Non-boundary frequencies
                  do ifreq=if0l+1,if1l-1

                    ! Get index
                    jjl = jjl + 1

                    ! Add to norm
                    daux = daux + data2(jjl,2)*Frec%W_freq(ifreq)

                  end do ! Non-boundary index

                  ! Single frequency
                  if (if1l.gt.if0l) then

                    ! Get index
                    jjl = jjl + 1

                    ! Add to norm
                    daux = daux + data2(jjl,2)*Atom(ia)%W1(jtran)

                  end if ! Single frequency

                  ! Non-zero norm
                  if (daux.gt.0d0) then

                    ! Normalize
                    daux = 1d0/daux
                    data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux

                  end if ! Non-zero norm
                end if ! MPI
              end if ! Iterating
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
            if (PRD.and.Atom(ia)%lemiss2(jtran).and. &
                iz.lt.Rz1_PRD) then

              ! Frec index
              jndx = Red%izao(jtran,ia,Rz0)
              indx = Red%izao(jtran,ia,iz)

              ! There are valid ranges
              if (Red%ao(jndx)%nran.gt.0) then

                ! Get limits for emissivity
                if0l2 = Red%ao(jndx)%gf0
                if1l2 = Red%ao(jndx)%gf1

                ! Valid range
                if (if1l2.ge.if0l2) then

                  ! If this height has only one dimension (LOS)
                  if (size(Red%zao(indx)%eps20,1).eq.1) then

                    ! Get values from corresponding direction
                    estmp0(if0l2:if1l2) = Red%zao(indx)%eps20(1,:)
                    estmp1(if0l2:if1l2) = Red%zao(indx)%eps21(1,:)
                    estmp2(if0l2:if1l2) = Red%zao(indx)%eps22(1,:)
                    estmp3(if0l2:if1l2) = Red%zao(indx)%eps23(1,:)

                  ! If this height has more than one direction
                  else

                    ! Get values from corresponding direction
                    estmp0(if0l2:if1l2) = Red%zao(indx)%eps20(jdir,:)
                    estmp1(if0l2:if1l2) = Red%zao(indx)%eps21(jdir,:)
                    estmp2(if0l2:if1l2) = Red%zao(indx)%eps22(jdir,:)
                    estmp3(if0l2:if1l2) = Red%zao(indx)%eps23(jdir,:)

                  end if ! Number of directions
                end if ! Valid ranges
              end if ! Valid ranges
            end if ! PRD

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

          ! If storing in RAM
          if (PIRAM) then

            ! Get b-f emissivity
            call photoepsS(Atom(ia),Frec%omega3(if0l:if1l), &
                           Frec%exu(if0l:if1l,iz), &
                           Atmo%T(iz),Atmo%ne(iz),jtran, &
                           ilevelu,iz,if0l,if1l,estmp0(if0l:if1l), &
                           rstmp1(if0l:if1l))

          ! Not storing in RAM
          else

            ! Get b-f emissivity
            call photoeps(Atom(ia),Frec%omega,Atmo%T(iz), &
                          Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                          estmp0(if0l:if1l),rstmp1(if0l:if1l))

          end if ! Storing in RAM

          ! Add contribution to emissivity
          epsA(if0l:if1l,0) = estmp0(if0l:if1l) + epsA(if0l:if1l,0)


          ! Get b-f Absorptivity
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
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

        ! Non-axial RT
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
      if (associated(p_Norm)) nullify(p_Norm)

      return

      end subroutine RTCoeff

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the opacity for a given height and direction for
      !! the polarized problem\n
      !!       Frec(Frequency_class): Structure with frequency data\n
      !!              Red(Red_class): Structure with redistribution
      !!                              input frequency data,
      !!                              redistribution function data,
      !!                              and profile or normalization
      !!                              data\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!                vfac(double): Doppler shift factor\n
      !!          Flgsg(Fctsg_class): Structure with factorials,
      !!                              signs, and J-symbols\n
      !!                 iz(integer): Height index\n
      !!               jdir(integer): Current direction index\n
      !!                if0(integer): Lower limit index for
      !!                              frequency\n
      !!                if1(integer): Upper limit index for
      !!                              frequency\n
      !!               cdir(integer): Direction index for background
      !!                              opacities\n
      !!         Cont(double(:,:,:)): Background opacity data\n
      !!        Bfield(Bfield_class): Structure with magnetic field
      !!                              data\n
      !!        TKQo(dcomplx(:,:,:)): Geometrical tensors in the
      !!                              suitable reference frame\n
      !!             etaI(double(:)): Opacity
      subroutine RTAbs(Frec,Red,Atom,LTElines,Atmo,vfac,Flgsg, &
                       iz,jdir,if0,if1,cdir,Cont,Bfield,TKQo, &
                       etaI)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(inout):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,cdir,if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(if0:if1,3,cdir), intent(in):: Cont
      double precision, dimension(if0:if1), intent(out):: etaI
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQo

      ! Local

      integer:: iterml,itermu,ia,jtran,jdir,icdir,ilevell,ilevelu
      integer:: if0l,if1l,indx,ktran

      double precision:: DwT,pE,absK,Dw
      double precision, dimension(if0:if1):: etaA
      double precision, dimension(if0:if1):: etmp0,estmp0
      double precision, dimension(if0:if1):: etmp1,estmp1
      double precision, dimension(if0:if1):: rstmp1,rtmp1
      double precision, dimension(if0:if1):: etmp2,estmp2
      double precision, dimension(if0:if1):: rstmp2,rtmp2
      double precision, dimension(if0:if1):: etmp3,estmp3
      double precision, dimension(if0:if1):: rstmp3,rtmp3

      ! Pointer

      type(Prof_class), pointer:: p_Norm


      ! Initialize
      nullify(p_Norm)

      ! Initialize continuum index
      icdir = 1

      ! If there are dynamics, correct continuum index
      if (dyn) icdir = min(jdir,cdir)


      !
      ! Continuum contribution
      !

      ! Continuum absorptivity
      etaI = Cont(:,1,icdir)


      !
      ! LTE lines
      !

      ! For each LTE line
      do ia=1,nLTEl

        ! If this CPU does not have frequencies in this line, skip
        if (LTElines(ia)%absent) cycle

        ! Skip if limited
        if (iz.lt.LTElines(ia)%Rz0) cycle

        ! Thermal part of the Doppler width
        DwT = LTElines(ia)%cDopp*sqrt(Atmo%T(iz))

        ! Store frequency limits
        if0l = LTElines(ia)%if0
        if1l = LTElines(ia)%if1

        ! Point to this line norm
        indx = Red%idzao(nxtran+ia,iz,1)
        p_Norm => Red%dzao(indx)

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

            ! Get RT coefficients
            call rt1ordLTE(LTElines(ia),TKQo,Frec%omega,Flgsg,iz, &
                           if0l,if1l,p_Norm,Dw,vfac, &
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

            ! RT coeffs
            call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,pE, &
                            etmp0(if0l:if1l),estmp0(if0l:if1l))

          end if ! Magnetic field

          ! Correct for stimulated emission
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                             estmp0(if0l:if1l)/absK

        ! If no stimulated
        else

          ! If magnetic field
          if (Bfield%Bstrength(iz).gt.TINYB) then

            ! Get absorption
            call absorbLTE(LTElines(ia),TKQo,Frec%omega,Flgsg,iz, &
                           if0l,if1l,p_Norm,Dw,vfac, &
                           Bfield%Bstrength(iz),pE, &
                           etmp0(if0l:if1l),etmp1(if0l:if1l), &
                           etmp2(if0l:if1l),etmp3(if0l:if1l), &
                           rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                           rtmp3(if0l:if1l))

          ! No magnetic field
          else

            ! RT coeffs
            call absorbILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,pE,etmp0(if0l:if1l))

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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Point to this line norm
            indx = Red%idzao(ktran,iz,1)
            p_Norm => Red%dzao(indx)

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

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
              ! First order RT coefficients
              !
              call rt1ord(Atom(ia),TKQo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK, &
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

              !
              ! Absorptivity
              !
              call absorb(Atom(ia),TKQo,Frec%omega,Flgsg, &
                          jtran,itermu,iterml,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,absK, &
                          etmp0(if0l:if1l),etmp1(if0l:if1l), &
                          etmp2(if0l:if1l),etmp3(if0l:if1l), &
                          rtmp1(if0l:if1l),rtmp2(if0l:if1l), &
                          rtmp3(if0l:if1l))

            end if ! Stimulated

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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Point to this line norm
            indx = Red%idzao(ktran,iz,1)
            p_Norm => Red%dzao(indx)

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

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
              call rt1ordNB(Atom(ia),TKQo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK, &
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
              call absorbNB(Atom(ia),TKQo,Frec%omega,Flgsg, &
                            jtran,itermu,iterml,iz,if0l,if1l, &
                            p_Norm,Dw,vfac,absK, &
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

          ! If storing in RAM
          if (PIRAM) then

            ! Get b-f emissivity
            call photoepsS(Atom(ia),Frec%omega3(if0l:if1l), &
                           Frec%exu(if0l:if1l,iz),Atmo%T(iz), &
                           Atmo%ne(iz),jtran,ilevelu,iz, &
                           if0l,if1l,estmp0(if0l:if1l), &
                           rstmp1(if0l:if1l))

          ! Not storing in RAM
          else

            ! Get b-f emissivity
            call photoeps(Atom(ia),Frec%omega,Atmo%T(iz), &
                          Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                          estmp0(if0l:if1l),rstmp1(if0l:if1l))

          end if ! Storing in RAM

          ! Get b-f Absorptivity
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
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
      if (associated(p_Norm)) nullify(p_Norm)

      return

      end subroutine RTAbs

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the RT coefficients at a given point in a CLE
      !! synthesis\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!            Red(Red_class): Structure with redistribution
      !!                            input frequency data,
      !!                            redistribution function data,
      !!                            and profile or normalization
      !!                            data\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!  GeomP(Coronapoint_class): Structure with geometric data for
      !!                            a CLE node\n
      !!              vfac(double): Doppler shift factor\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        TS(dcomplx(:,:,:)): Geometrical tensors in the
      !!                            vertical reference frame\n
      !!        TB(dcomplx(:,:,:)): Geometrical tensors in the
      !!                            suitable reference frame\n
      !!              if0(integer): Lower limit index for frequency\n
      !!              if1(integer): Upper limit index for frequency\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors integrated
      !!                            over the absorption profile\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!        spect(spect_class): Structure with the input spectra
      !!                            data\n
      !!       Cont(double(:,:,:)): Background opacity data\n
      !!     add_cont_cle(logical): If including continuum
      !!                            contributions\n
      !!        data1(double(:,:)): Radiation transfer coefficients
      subroutine RTCoeff_CLE(Frec,Red,Atom,Atmo,Flgsg,Geom,GeomP, &
                             vfac,Bfield,TS,TB,if0,if1,JKQ,JKQC, &
                             spect,Cont,add_cont_cle,data1)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Coronapoint_class), intent(in):: GeomP
      type(Bfield_class), intent(in):: Bfield
      type(Spect_class), intent(in):: spect
      logical, intent(in):: add_cont_cle
      integer, intent(in):: if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(if0:if1,3), intent(in):: Cont
      double precision, dimension(0:3,if0:if1,0:4), &
                        intent(out):: data1
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TS,TB
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(in):: JKQ
      complex(kind=8), dimension(-2:2,0:2,if0:if1), intent(in):: JKQC

      ! Local

      integer:: iS,K,iQ,ia,ilevell,ilevelu,iterml,itermu,jtran,ktran
      integer:: ifreq,if0l,if1l!,if0l2,if1l2

      double precision:: DwT,Dw,pE,absK
      double precision, dimension(if0:if1,0:3):: etaA,epsA,rhsA,rhaA
      double precision, dimension(if0:if1):: etmp0,estmp0!,es2tmp0
      double precision, dimension(if0:if1):: etmp1,estmp1
      double precision, dimension(if0:if1):: rstmp1,rtmp1!,es2tmp1
      double precision, dimension(if0:if1):: etmp2,estmp2
      double precision, dimension(if0:if1):: rstmp2,rtmp2!,es2tmp2
      double precision, dimension(if0:if1):: etmp3,estmp3
      double precision, dimension(if0:if1):: rstmp3,rtmp3!,es2tmp3
      double precision, dimension(if0:if1):: intgr

      ! Pointer

      type(Prof_class), pointer:: p_Norm


      !
      ! Initialize
      !
      nullify(p_Norm)

      ! Get dummy norm
      p_Norm => Red%dzao(1)


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

        ! Zero out
        data1(0,:,0) = 0d0
        data1(:,:,4) = 0d0

      end if ! Continuum contribution


      ! If exu is associated
      if (associated(Frec%exu)) then

        ! Initialize argument
        Frec%exu(:,1) = Frec%omega_ou(Frec%pif0:Frec%pif1)* &
                        (c2*1d4/Atmo%T(1))

        ! For every relevant frequency
        do ifreq=Frec%pif0,Frec%pif1

          ! Get inverse exponential
          Frec%exu(ifreq,1) =  diexp(Frec%exu(ifreq,1))

        end do ! For every relevant frequency

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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

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
            ! First order RT coefficients
            !
            call rt1ord(Atom(ia),TB,Frec%omega,Flgsg, &
                        jtran,itermu,iterml,1,if0l,if1l, &
                        p_Norm,Dw,vfac,absK, &
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
            ! TODO TODO

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

            ! Rolling index
            ktran = jtran + Atom(ia)%tshift

            ! Store frequency limits
            if0l = Atom(ia)%if0(jtran)
            if1l = Atom(ia)%if1(jtran)

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
                          p_Norm,Dw,vfac,absK, &
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
            ! TODO TODO

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
          call photoabsI(Atom(ia),jtran,ilevell,1,if0l,if1l, &
                         etmp0(if0l:if1l))

          ! Remove the stimulated part
          etmp0(if0l:if1l) = etmp0(if0l:if1l) - rstmp1(if0l:if1l)

          ! Add contribution to absorptivity
          etaA(if0l:if1l,0) = etmp0(if0l:if1l) + etaA(if0l:if1l,0)

        end do ! b-f transitions

        ! Multiply by the population of the atom the contribution of
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
      if (associated(p_Norm)) nullify(p_Norm)

      return

      ! Fake assign to deceive the compiler
      iS = Geom%nTh
      DwT = GeomP%theta
      DwT = dble(JKQ(-2,0,1))
      iS = Spect%nfreq

      end subroutine RTCoeff_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeff_mod
