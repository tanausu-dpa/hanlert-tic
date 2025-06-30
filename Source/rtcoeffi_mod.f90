      !> Compute radiation transfer coefficients for NLTE problem of
      !! the first kind
      module rtcoeffi_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!     Roberto Casini (HAO)
!  Start:
!     20/04/2017
!  Last version:
!     26/06/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     26/06/2025:    V4.0.3 - Updated due to changes in the Red_class
!                             and Redb_class structures (TdPA)
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
!  RTCoeffI
!    Calculate the radiation transfer coefficients for a given height
!  and direction for the unpolarized problem
!
!  RTAbsI
!    Calculate the opacity for a given height and direction for the
!  unpolarized problem
!
!  RTCoeffJ
!    Calculate the continuum radiation transfer coefficients for a
!  given height and direction for the unpolarized problem
!
!  Termprof
!    Calculate the transition profiles for a given height and
!  direction for the polarized problem without PRD
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use parameters_mod , only : convF , c , vacuum , TINYB , &
                                  VTINY , IPI
      use rtcoeffaux_mod
      use rtcoeffiaux_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the radiation transfer coefficients for a given
      !! height and direction for the unpolarized problem\n
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
      !!                 iz(integer): Height index\n
      !!               jdir(integer): Current direction index\n
      !!              nodir(integer): Corrected direction index for
      !!                              Red pointers\n
      !!                if0(integer): Lower limit index for
      !!                              frequency\n
      !!                if1(integer): Upper limit index for
      !!                              frequency\n
      !!             J00C(double(:)): Mean intensity with frequency
      !!                              dependence\n
      !!               cdir(integer): Direction index for background
      !!                              opacities\n
      !!         Cont(double(:,:,:)): Background opacity data\n
      !!               rl(double(:)): Bound-bound transition
      !!                              strengths\n
      !!               rp(double(:)): Bound-free transition
      !!                              strengths\n
      !!          data1(double(:,:)): Radiation transfer
      !!                              coefficients\n
      !!          data2(double(:,:)): Line profiles\n
      !!          iterating(logical): If solving self-consistent
      !!                              problem
      subroutine RTCoeffI(Frec,Red,Atom,LTElines,Atmo,vfac,iz, &
                          jdir,nodir,if0,if1,J00C,cdir,Cont,rl,rp, &
                          data1,data2,iterating)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      logical, intent(in):: iterating
      integer, intent(in):: iz,jdir,nodir,cdir,if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(if0:if1,3,cdir), intent(in):: Cont
      double precision, dimension(nfreq), intent(in):: J00C
      double precision, dimension(:), intent(inout):: rl
      double precision, dimension(:), intent(inout):: rp
      double precision, dimension(if0:if1,0:1), intent(out):: data1
      double precision, dimension(:,:), intent(inout):: data2

      ! Local

      integer:: iterml,itermu,iJl,iJu,ilevell,ilevelu,ifreq,nf2
      integer:: ia,jtran,ktran,fjtran,ffktran,ffjtran,indx,jndx
      integer:: icdir,if0l,if1l,if0l2,if1l2,iil,jjl,iip,nf

      double precision:: daux,DwT,Dfreq,pE,absK,Dw,pop,rhou
      double precision, dimension(if0:if1):: rpf,etmp,estmp,rstmp

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

      ! Absorptivity
      data1(:,0) = Cont(:,1,icdir)

      ! Emissivity
      data1(:,1) = Cont(:,3,icdir) + Cont(:,2,icdir)*J00C(if0:if1)

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

        ! Point to this line norm
        indx = Red%idzao(nxt+ia,iz,nodir)
        p_Norm => Red%dzao(indx)

        ! Get frequency of FS transition
        Dfreq = LTElines(ia)%Eu - LTElines(ia)%El

        ! Photon energy (cgs) and convertion factor
        pE = convF*Dfreq
        absK = 1d21*(2d0*c)*Dfreq**2d0

        ! Add the microturbulence to Doppler width
        Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)

        ! RT coeffs
        call rt1ordILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                        p_Norm,Dw,vfac,pE, &
                        etmp(if0l:if1l),estmp(if0l:if1l))

        !
        ! Stimulated emission contribution
        !

        ! If there is stimulated emission, subtract from
        ! absorption
        if (stm) &
          etmp(if0l:if1l) = etmp(if0l:if1l) - &
                            estmp(if0l:if1l)/absK

        ! Add the contribution to the absorptivity of this atom
        data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                             etmp(if0l:if1l)*pop

        ! Add the contribution to the emissivity of this atom
        data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                             estmp(if0l:if1l)*(pE*pop)

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

          ! Rolling index
          ktran = jtran + Atom(ia)%tshift

          ! Store frequency limits
          if0l = Atom(ia)%if0(jtran)
          if1l = Atom(ia)%if1(jtran)
          nf = if1l - if0l

          ! Identify involved terms
          itermu = Atom(ia)%fst(jtran)%itermu
          iterml = Atom(ia)%fst(jtran)%iterml

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
            ffktran = ffjtran + Atom(ia)%tfshift

            ! Point to this line norm
            indx = Red%idzao(ffktran,iz,nodir)
            p_Norm => Red%dzao(indx)

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJl,iterml)

            ! Photon energy (cgs) and convertion factor
            pE = convF*Dfreq
            absK = 1d21*(2d0*c)*Dfreq**2d0

            ! Add the microturbulence to Doppler width
            Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


            !
            ! RT coeffs
            !
            call rt1ordI(Atom(ia),Frec%omega,jtran, &
                         itermu,iterml,iJu,iJl,iz, &
                         if0l,if1l,p_Norm,Dw,vfac,pE, &
                         etmp(if0l:if1l),estmp(if0l:if1l),rhou)


            !
            ! Store absorption profile
            !
            if (iterating) then

              ! Storing
              if (p_Norm%VRAM) then

                data2(iil:iil+nf,1) = p_Norm%p*1d-5*sqrt(IPI)/Dw

              ! Computing
              else

                ! If MPI
                if (nproc.gt.1) then

                  ! Copy absorbtivity in data2
                  data2(iil:iil+nf,1) = etmp(if0l:if1l)

                ! If there is no MPI, normalize it here
                else

                  ! Initialize
                  daux = 0d0

                  ! Copy absorbtivity in data2
                  data2(iil:iil+nf,1) = etmp(if0l:if1l)

                  ! For each line frequency
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

                  ! If non-zero norm
                  if (daux.gt.0d0) then

                    ! Normalize
                    daux = 1d0/daux
                    data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux

                  end if ! Non zero normalization
                end if ! Not MPI
              end if ! Type of profile
            end if ! Iterating


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp(if0l:if1l) = etmp(if0l:if1l) - &
                                estmp(if0l:if1l)/absK

              !
              ! Store emission profile
              !

              ! Self-consistent solution
              if (iterating) then

                ! For intensity it is the same than absorption
                data2(iil:iil+nf,2) = data2(iil:iil+nf,1)

              end if ! Self-consistent solution
            endif ! Stimulated emission

            ! Add the contribution to the absorptivity of this atom
            data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                                 etmp(if0l:if1l)*pop

            ! Store the numerator of the opacity fraction of this
            ! line if self-consistent solution
            if (iterating) &
              rL(iil:iil+nf) = estmp(if0l:if1l)*pE*pop/rhou

            !
            ! Second order emissivity
            !
            if (PRD.and.Atom(ia)%lemiss2(jtran).and. &
                iz.lt.Rz1_PRD) then

              ! Red index
              jndx = Red%izao(ffjtran,ia,Rz0)
              indx = Red%izao(ffjtran,ia,iz)

              ! There are valid ranges
              if (Red%ao(jndx)%nran.gt.0) then

                ! Get limits for emissivity
                if0l2 = Red%ao(jndx)%gf0
                if1l2 = Red%ao(jndx)%gf1

                ! If valid range
                if (if1l2.ge.if0l2) then

                  ! If this height has only one dimension
                  if (size(Red%zao(indx)%eps20,1).eq.1) then

                    ! Get values from direction 1
                    estmp(if0l2:if1l2) = Red%zao(indx)%eps20(1,:)

                    ! If self-consistent solution, get correction
                    if (iterating) &
                      rpf(if0l2:if1l2) = Red%zao(indx)%rpf(1,:)

                  ! If this height has more than one direction
                  else

                    ! Get values from corresponding direction
                    estmp(if0l2:if1l2) = Red%zao(indx)%eps20(jdir,:)

                    ! If self-consistent solution, get correction
                    if (iterating) &
                      rpf(if0l2:if1l2) = Red%zao(indx)%rpf(jdir,:)

                  end if ! Number of directions

                  ! If doing self-consistent solution
                  if (iterating) then

                    ! Introduce the PRD factor
                    jjl = iil + if0l2 - if0l
                    nf2 = if1l2 - if0l2
                    rL(jjl:jjl+nf2) = rL(jjl:jjl+nf2)*rpf(if0l2:if1l2)

                  end if ! Self-consistent solution
                end if ! Valid range
              end if ! Valid ranges
            end if ! PRD

            ! Add the contribution to the emissivity of this atom
            data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                                 estmp(if0l:if1l)*pE*pop

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

          ! Store frequency limits
          if0l = Atom(ia)%phot(jtran)%if0
          if1l = Atom(ia)%phot(jtran)%if1
          nf = if1l - if0l

          ! Identify the levels involved
          ilevelu = Atom(ia)%phot(jtran)%ilevelu
          ilevell = Atom(ia)%phot(jtran)%ilevell

          ! If storing in RAM
          if (PIRAM) then

            ! Get b-f emissivity
            call photoepsIS(Atom(ia),Frec%omega3(if0l:if1l), &
                            Frec%exu(if0l:if1l,iz), &
                            Atmo%T(iz),Atmo%ne(iz), &
                            jtran,ilevelu,iz,if0l,if1l, &
                            estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          ! Not storing in RAM
          else

            ! Get b-f emissivity
            call photoepsI(Atom(ia),Frec%omega, &
                           Atmo%T(iz),Atmo%ne(iz), &
                           jtran,ilevelu,iz,if0l,if1l, &
                           estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          end if ! Storing in RAM

          ! Get b-f Absorptivity
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                         etmp(if0l:if1l))

          ! ALI ratio
          if (iterating) &
            rP(iip:iip+nf) = estmp(if0l:if1l)*pop/rhou

          ! Add contribution to emissivity
          data1(if0l:if1l,1) = data1(if0l:if1l,1) + &
                               estmp(if0l:if1l)*pop

          ! Remove the stimulated part
          etmp(if0l:if1l) = etmp(if0l:if1l) - rstmp(if0l:if1l)

          ! Add contribution to absorptivity
          data1(if0l:if1l,0) = data1(if0l:if1l,0) + &
                               etmp(if0l:if1l)*pop

          ! Update the index
          iip = iip + nf + 1

        end do ! b-f transitions
      end do ! Atom


      !
      ! Finish the construction of the opacity fractions
      !
      if (iterating) then

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
              rL(iil:iil+nf) = rL(iil:iil+nf)/ &
                               (data1(if0l:if1l,0) + vacuum)

              ! Update index
              iil = iil + nf + 1

            end do ! b-b FS transitions
          end do ! b-b transitions

          ! For each b-f transition
          do jtran=1,Atom(ia)%nphot

            ! If this CPU does not have frequencies in this
            ! transition, skip
            if (Atom(ia)%phot(jtran)%absent) cycle

            ! Store frequency limits
            if0l = Atom(ia)%phot(jtran)%if0
            if1l = Atom(ia)%phot(jtran)%if1
            nf = if1l - if0l

            ! Get the opacity fraction
            rP(iip:iip+nf) = rP(iip:iip+nf)/ &
                             (data1(if0l:if1l,0) + vacuum)

            ! Update index
            iip = iip + nf + 1

          end do ! b-f transitions
        end do ! Atoms

      end if ! Formal solution


      !
      ! Transform into the data arrays
      !

      ! Source function
      data1(:,1) = data1(:,1)/(data1(:,0) + vacuum)

      ! Nullify pointers
      if (associated(p_Norm)) nullify(p_Norm)

      return

      end subroutine RTCoeffI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the opacity for a given height and direction for
      !! the unpolarized problem\n
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
      !!                 iz(integer): Height index\n
      !!               jdir(integer): Current direction index\n
      !!                if0(integer): Lower limit index for
      !!                              frequency\n
      !!                if1(integer): Upper limit index for
      !!                              frequency\n
      !!             J00C(double(:)): Mean intensity with frequency
      !!                              dependence\n
      !!               cdir(integer): Direction index for background
      !!                              opacities\n
      !!         Cont(double(:,:,:)): Background opacity data\n
      !!             etaI(double(:)): Opacity
      subroutine RTAbsI(Frec,Red,Atom,LTElines,Atmo,vfac,iz, &
                        jdir,if0,if1,cdir,Cont,etaI)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      integer, intent(in):: iz,jdir,cdir,if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(if0:if1,3,cdir), intent(in):: Cont
      double precision, dimension(if0:if1), intent(out):: etaI

      ! Local

      integer:: iterml,itermu,iJl,iJu,ilevell,ilevelu,indx
      integer:: ia,jtran,fjtran,ffjtran,ffktran,icdir,if0l,if1l

      double precision:: DwT,Dfreq,pE,absK,Dw,pop,rhou
      double precision, dimension(if0:if1):: etmp,estmp,rstmp

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

        ! Point to this line norm
        indx = Red%idzao(nxt+ia,iz,1)
        p_Norm => Red%dzao(indx)

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
                          p_Norm,Dw,vfac,pE, &
                          etmp(if0l:if1l),estmp(if0l:if1l))

          !
          ! Stimulated emission contribution
          !
          etmp(if0l:if1l) = etmp(if0l:if1l) - &
                            estmp(if0l:if1l)/absK

        ! No stimulated emission
        else

          !
          ! Absorptivity
          !
          call absorbILTE(LTElines(ia),Frec%omega,iz,if0l,if1l, &
                          p_Norm,Dw,vfac,pE,etmp(if0l:if1l))

        end if

        ! Add the contribution to the absorptivity of this atom
        etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)

      end do ! LTE lines


      !
      ! NLTE lines
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

          ! For each FS transition
          do fjtran=1,Atom(ia)%fst(jtran)%nt

            ! Idenfity involved levels
            iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
            iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

            ! Get the sequential index of this FS transition
            ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)
            ffktran = ffjtran + Atom(ia)%tfshift

            ! Point to this line norm
            indx = Red%idzao(ffktran,iz,1)
            p_Norm => Red%dzao(indx)

            ! Get frequency of FS transition
            Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                    Atom(ia)%FSfreq(iJl,iterml)

            ! Photon energy (cgs) and convertion factor
            pE = convF*Dfreq
            absK = 1d21*(2d0*c)*Dfreq**2d0

            ! Add the microturbulence to Doppler width
            Dw = Dfreq*sqrt(DwT*DwT + Atmo%vmi(iz)**2d0)


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! RT coeffs
              call rt1ordI(Atom(ia),Frec%omega, &
                           jtran,itermu,iterml,iJu,iJl,iz, &
                           if0l,if1l,p_Norm,Dw,vfac,pE, &
                           etmp(if0l:if1l),estmp(if0l:if1l),rhou)

              ! Subtract from absorption
              etmp(if0l:if1l) = etmp(if0l:if1l) - &
                                estmp(if0l:if1l)/absK

            ! No stimulated emission
            else

              ! Absorptivity
              call absorbI(Atom(ia),Frec%omega, &
                           jtran,itermu,iterml,iJu,iJl,iz, &
                           if0l,if1l,p_Norm,Dw,vfac,pE, &
                           etmp(if0l:if1l))

            end if

            ! Add the contribution to the absorptivity of this atom
            etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)

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

          ! If storing photoionization quantities in RAM
          if (PIRAM) then

            ! Get b-f emissivity
            call photoepsIS(Atom(ia),Frec%omega3(if0l:if1l), &
                            Frec%exu(if0l:if1l,iz), &
                            Atmo%T(iz),Atmo%ne(iz), &
                            jtran,ilevelu,iz,if0l,if1l, &
                            estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          ! Not storing in RAM
          else

            ! Get b-f emissivity
            call photoepsI(Atom(ia),Frec%omega,Atmo%T(iz), &
                           Atmo%ne(iz),jtran,ilevelu,iz,if0l,if1l, &
                           estmp(if0l:if1l),rstmp(if0l:if1l),rhou)

          end if ! Storing in RAM

          ! Get b-f absorptivity
          call photoabsI(Atom(ia),jtran,ilevell,iz,if0l,if1l, &
                         etmp(if0l:if1l))

          ! Remove the stimulated part
          etmp(if0l:if1l) = etmp(if0l:if1l) - rstmp(if0l:if1l)

          ! Add contribution to absorptivity
          etaI(if0l:if1l) = etmp(if0l:if1l)*pop + etaI(if0l:if1l)

        end do ! b-f transitions
      end do ! Atom

      ! Nullify pointers
      if (associated(p_Norm)) nullify(p_Norm)

      return

      end subroutine RTAbsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the continuum radiation transfer coefficients for
      !! a given height and direction for the unpolarized problem\n
      !!        jdir(integer): Current direction index\n
      !!         if0(integer): Lower limit index for frequency\n
      !!         if1(integer): Upper limit index for frequency\n
      !!      J00C(double(:)): Mean intensity with frequency
      !!                       dependence\n
      !!        cdir(integer): Direction index for background
      !!                       opacities\n
      !!  Cont(double(:,:,:)): Background opacity data\n
      !!   data1(double(:,:)): Radiation transfer coefficients\n
      subroutine RTCoeffJ(jdir,if0,if1,J00C,cdir,Cont,data1)

      ! I/O

      integer, intent(in):: jdir,cdir,if0,if1
      double precision, dimension(if0:if1,3,cdir), intent(in):: Cont
      double precision, dimension(nfreq), intent(in):: J00C
      double precision, dimension(if0:if1,0:1), intent(out):: data1

      ! Local

      integer:: icdir


      ! Index of continuum direction
      icdir = 1

      ! If there are dynamics
      if (dyn) icdir = min(jdir,cdir)

      !
      ! Calculate RT coefficients
      !

      ! Absorptivity
      data1(:,0) = Cont(:,1,icdir)

      ! Emissivity
      data1(:,1) = Cont(:,3,icdir) + Cont(:,2,icdir)*J00C(if0:if1)

      ! Source function
      data1(:,1) = data1(:,1)/(data1(:,0) + vacuum)

      return

      end subroutine RTCoeffJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the transition profiles for a given height and
      !! direction for the polarized problem without PRD\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!         Red(Red_class): Structure with redistribution input
      !!                         frequency data, redistribution
      !!                         function data, and profile or
      !!                         normalization data\n
      !!    Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!           vfac(double): Doppler shift factor\n
      !!     Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                         J-symbols\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!   Bfield(Bfield_class): Structure with magnetic field data\n
      !!            iz(integer): Height index\n
      !!          jdir(integer): Current direction index\n
      !!           if0(integer): Lower limit index for frequency\n
      !!           if1(integer): Upper limit index for frequency\n
      !!    TSo(dcomplx(:,:,:)): Geometrical tensors in the vertical
      !!                         reference frame\n
      !!   TKQo(dcomplx(:,:,:)): Geometrical tensors in the suitable
      !!                         reference frame\n
      !!     data2(double(:,:)): Line profiles
      subroutine Termprof(Frec,Red,Atom,Atmo,vfac,Flgsg,Geom, &
                          Bfield,iz,jdir,if0,if1,TSo,TKQo, &
                          data2)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Frequency_class), intent(in):: Frec
      type(Red_class), intent(in):: Red
      type(Fctsg_class), intent(in):: Flgsg
      type(Geometry_class), intent(in):: Geom
      type(Bfield_class), intent(in):: Bfield
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: vfac
      double precision, dimension(:,:), intent(inout):: data2
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TSo
      complex(kind=8), dimension(0:3,-2:2,0:2), intent(in):: TKQo

      ! Local

      integer:: iterml,itermu,ia,jtran,ktran,jdir
      integer:: ifreq,if0l,if1l,iil,jjl,nf,nodir,indx

      double precision:: daux,DwT,absK,Dw
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

      ! Correct direction index
      nodir = min(jdir,ubound(Red%idzao,3),Geom%njdir)

      ! Initialize index
      iil = 1

      ! For each atom
      do ia=1,nA

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(Atmo%T(iz))

        !
        ! If magnetic field
        !
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

            ! Copy absorbtivity in data2
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

              ! If non-zero norm
              if (daux.gt.0d0) then

                ! Normalize
                daux = 1d0/daux
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux

              end if ! Non-zero norm
            end if ! MPI


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK

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

                ! If non-zero norm
                if (daux.gt.0d0) then

                  ! Normalize
                  daux = 1d0/daux
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux

                end if ! Non-zero norm
              end if ! MPI
            endif ! Stimulated emission

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        !
        ! No field
        !
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
            absK = 1d21*(2d0*c)*Atom(ia)%Dfreq(jtran)**2d0

            ! Add the microt. to Doppler width
            Dw = Atom(ia)%Dfreq(jtran)*sqrt(DwT*DwT + &
                                            Atmo%vmi(iz)**2d0)
            !
            ! First order RT coefficients
            !

            ! Get absorption and emission
            call rt1ordNB(Atom(ia),TSo,Frec%omega,Flgsg, &
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

            ! Copy absorbtivity in data2
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

              ! Non-zero norm
              if (daux.gt.0d0) then

                ! Normalize profile
                daux = 1d0/daux
                data2(iil:iil+nf,1) = data2(iil:iil+nf,1)*daux

              end if ! Non-zero norm
            end if ! MPI


            !
            ! Stimulated emission contribution
            !

            ! If there is stimulated emission
            if (stm) then

              ! Correct for stimulated emission
              etmp0(if0l:if1l) = etmp0(if0l:if1l) - &
                                 estmp0(if0l:if1l)/absK

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

                  ! Normalize profile
                  daux = 1d0/daux
                  data2(iil:iil+nf,2) = data2(iil:iil+nf,2)*daux

                end if ! Non-zero norm
              end if ! MPI
            endif ! Stimulated emission

            ! Advance the index
            iil = iil + nf + 1

          end do ! b-b transitions

        end if ! Magnetic field

      end do ! Atoms

      ! Nullify pointers
      if (associated(p_Norm)) nullify(p_Norm)

      return

      end subroutine Termprof

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffi_mod
