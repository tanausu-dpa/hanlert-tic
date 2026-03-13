      !> Radiation transfer coefficients
      module rtcoeffiaux_mod
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
!     12/03/2026 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     12/03/2026:    V4.0.2 - Bugfix: Added a missing normalization
!                             factor in rt1ordI (TdPA)
!                           - Added the sourI subroutine (TdPA)
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
!  Info:
!
!   - Units for absorption are cm^-1 (true absorption coefficient).
!     Needs to be multiplied by the actual atomic density.
!
!   - Units for emissivity are given in number of photons per unit
!     interval of time (s) and normalized frequency (in units of
!     Doppler width), emitted by a unit volume of gas (cm^-3) of unit
!     atomic density, within one steradian. In order to compute
!     photometric values of the intensity, the output needs be
!     multiplied by the actual atomic density.
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  rt1ordI
!    Calculate the absorptivity and emissivity for the intensity for
!  a given atomic line
!
!  absorbI
!    Calculate the absorptivity for the intensity for a given atomic
!  line
!
!  emissI
!    Calculate the emissivity for the intensity for a given atomic
!  line in the comoving frame
!
!  sourI
!    Calculate the emissivity and absorptivity for the intensity for
!  a given atomic line in the comoving frame
!
!  emissI2ord
!    Calculate the second order emissivity of a given atomic line.
!  This subroutine only computes the positive contribution of the
!  coherent scattering and the negative contribution for the
!  flat-spectrum, i.e., the result needs to be added to the product
!  of emissI to get the actual emissivity
!
!  get_WarrI
!    Calculate the redistribution function
!
!  rt1ordILTE
!    Calculate the absorptivity and emissivity for the intensity for
!  a given LTE line
!
!  absorbILTE
!    Calculate the absorptivity for the intensity for a given LTE line
!
!  emissILTE
!    Calculate the emissivity for the intensity for a given LTE line
!
!  photoabsI
!    Calculate the absorptivity of a given photoionization transition
!
!  photoepsI
!    Calculate the emissivity of a given recombination transition
!
!  photoepsIS
!    Calculate the emissivity of a given recombination transition with
!  frequency quantities stored in RAM
!
!  getJMV
!    Calculate the frequency dependent mean intensity for the
!  angle-averaged second order emissivity in the presence of
!  velocities in the comoving frame
!
!  getStkinInu
!    Interpolate the intensity into the requested frequency
!
!  getStkinI
!    Interpolate the Stokes parameters into the input frequency axis
!
!  getJin
!    Interpolate the frequency dependent mean intensity into the input
!  frequency axis
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use math_mod
      use parameters_mod , only : IPI, c, c2, convF, cSaha, kb, &
                                  fktoJ, cZero, bigexp, IPI41, &
                                  IPI42, IPI2, IPI4, TINYO, TINYWAR, &
                                  vrfrac
      use profile_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity and emissivity for the intensity
      !! for a given atomic line\n
      !!   Atom(Atom_class): Structure with atomic data\n
      !!   omega(double(:)): Frequency array\n
      !!     itran(integer): Index of transition to compute\n
      !!    itermu(integer): Upper term of the transition\n
      !!    iterml(integer): Lower term of the transition\n
      !!       iJu(integer): Upper level of the transition\n
      !!       iJl(integer): Lower level of the transition\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!  Norma(Prof_class): Normalization factors for Voigt profiles
      !!                     or Voigt profiles\n
      !!         Dw(double): Doppler width of the transition\n
      !!       vfac(double): Doppler shift factor\n
      !!         pE(double): Unit transformation factor\n
      !!     eta(double(:)): Intensity absorptivity\n
      !!     eps(double(:)): Intensity emissivity\n
      !!       rhou(double): Factor for Lambda operator
      subroutine rt1ordI(Atom,omega,itran,itermu,iterml, &
                         iJu,iJl,iz,if0,if1,Norma,Dw,vfac,pE, &
                         eta,eps,rhou)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw,pE,vfac
      double precision, intent(out):: rhou
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq,i,iR

      double precision:: eu,el,au,al,aul,rho,at,feta,feps
      double precision:: iDw,Dfreqw,vfacw,prof


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      !
      ! Lower

      ! Level index
      i = Atom%irho(iterml)%irho_ij(iJl)

      ! Population
      rho = dble(Atom%popu(i,iz))

      ! Absorptibity factor
      feta = rho*1d3*IPI41*Atom%fst(itran)%Blu(iJl,iJu)*pE*iDw

      !
      ! Upper

      ! Level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Population
      rho = Atom%popu(i,iz)

      ! Emissivity factor
      feps = rho*1d3*IPI41*Atom%fst(itran)%Aul(iJu,iJl)*iDw

      ! If stored in RAM
      if (Norma%VRAM) then

        ! Use profile
        eps = feps*Norma%p(if0:if1)
        eta = feta*Norma%p(if0:if1)

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)

        ! Level quantities

        ! Damping parameter
        au = Atom%damp(itermu,iz)

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Damping parameter
        al = Atom%damp(iterml,iz)

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Intermediate quantities
        at = (au+al+aul)*iDw
        Dfreqw = (eu - el)*iDw
        vfacw = vfac*iDw
        feps = feps*Norma%Norm(1)
        feta = feta*Norma%Norm(1)

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          ! Save contribution
          eps(ifreq) = feps*prof
          eta(ifreq) = feta*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine rt1ordI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity for the intensity for a given
      !! atomic line\n
      !!   Atom(Atom_class): Structure with atomic data\n
      !!   omega(double(:)): Frequency array\n
      !!     itran(integer): Index of transition to compute\n
      !!    itermu(integer): Upper term of the transition\n
      !!    iterml(integer): Lower term of the transition\n
      !!       iJu(integer): Upper level of the transition\n
      !!       iJl(integer): Lower level of the transition\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!  Norma(Prof_class): Normalization factors for Voigt profiles
      !!                     or Voigt profiles\n
      !!         Dw(double): Doppler width of the transition\n
      !!       vfac(double): Doppler shift factor\n
      !!         pE(double): Unit transformation factor\n
      !!     eta(double(:)): Intensity absorptivity\n
      subroutine absorbI(Atom,omega,itran,itermu,iterml, &
                         iJu,iJl,iz,if0,if1,Norma,Dw,vfac,pE,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw,pE,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq,i

      double precision:: eu,el,au,al,aul,rho,at,feta
      double precision:: iDw,Dfreqw,vfacw,prof


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      ! Level index
      i = Atom%irho(iterml)%irho_ij(iJl)

      ! Population
      rho = dble(Atom%popu(i,iz))

      ! Absorptibity factor
      feta = rho*1d3*IPI41*Atom%fst(itran)%Blu(iJl,iJu)*pE*iDw

      ! If stored in RAM
      if (Norma%VRAM) then

        eta = feta*Norma%p(if0:if1)

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        al = Atom%damp(iterml,iz)

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Damping parameter
        au = Atom%damp(itermu,iz)

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)

        ! Intermediate quantities
        at = (au+al+aul)*iDw
        Dfreqw = (eu - el)*iDw
        vfacw = vfac*iDw
        feta = feta*Norma%Norm(1)

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          ! Save contribution
          eta(ifreq) = feta*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine absorbI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity for the intensity for a given atomic
      !! line in the comoving frame\n
      !!   Atom(Atom_class): Structure with atomic data\n
      !!   omega(double(:)): Frequency array\n
      !!     itran(integer): Index of transition to compute\n
      !!    itermu(integer): Upper term of the transition\n
      !!    iterml(integer): Lower term of the transition\n
      !!       iJu(integer): Upper level of the transition\n
      !!       iJl(integer): Lower level of the transition\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!  Norma(Prof_class): Normalization factors for Voigt profiles
      !!                     or Voigt profiles\n
      !!         Dw(double): Doppler width of the transition\n
      !!     eps(double(:)): Intensity emissivity
      subroutine emissI(Atom,omega,itran,itermu,iterml, &
                        iJu,iJl,iz,if0,if1,Norma,Dw,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps

      ! Local

      integer:: ifreq,i

      double precision:: el,eu,al,au,aul,rho,at,iDw,Dfreqw,prof,feps


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      ! Level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! Population
      rho = Atom%popu(i,iz)

      ! Emissivity factor
      feps = rho*1d3*IPI41*Atom%fst(itran)%Aul(iJu,iJl)*iDw

      ! If stored
      if (Norma%VRAM) then

        ! Get profile
        eps = feps*Norma%p(if0:if1)

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)

        ! Level quantities

        ! Damping parameter
        au = Atom%damp(itermu,iz)

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Damping parameter
        al = Atom%damp(iterml,iz)

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Intermediate quantities
        at = (au+al+aul)*iDw
        Dfreqw = (eu - el)*iDw
        feps = feps*Norma%Norm(1)

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*iDw,at,prof)

          ! Save contribution
          eps(ifreq) = feps*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine emissI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity and absorptivity for the intensity
      !! for a given atomic line in the comoving frame\n
      !!   Atom(Atom_class): Structure with atomic data\n
      !!   omega(double(:)): Frequency array\n
      !!     itran(integer): Index of transition to compute\n
      !!    itermu(integer): Upper term of the transition\n
      !!    iterml(integer): Lower term of the transition\n
      !!       iJu(integer): Upper level of the transition\n
      !!       iJl(integer): Lower level of the transition\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!  Norma(Prof_class): Normalization factors for Voigt profiles
      !!                     or Voigt profiles\n
      !!         Dw(double): Doppler width of the transition\n
      !!         pE(double): Unit transformation factor\n
      !!     eta(double(:)): Intensity absorptivity\n
      !!     eps(double(:)): Intensity emissivity
      subroutine sourI(Atom,omega,itran,itermu,iterml, &
                       iJu,iJl,iz,if0,if1,Norma,Dw,pE,eta,eps)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Prof_class), intent(in):: Norma
      integer, intent(in):: itran,itermu,iterml,iJu,iJl,iz
      integer, intent(in):: if0,if1
      double precision, intent(in):: Dw,pE
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps,eta

      ! Local

      integer:: ifreq,i

      double precision:: el,eu,al,au,aul,rho,at,iDw,Dfreqw,prof
      double precision:: feps,feta


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      !
      ! Lower

      ! Level index
      i = Atom%irho(iterml)%irho_ij(iJl)

      ! Population
      rho = Atom%popu(i,iz)

      ! Absorptibity factor
      feta = rho*1d3*IPI41*Atom%fst(itran)%Blu(iJl,iJu)*pE*iDw

      !
      ! Upper

      ! Level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! Population
      rho = Atom%popu(i,iz)

      ! Emissivity factor
      feps = rho*1d3*IPI41*Atom%fst(itran)%Aul(iJu,iJl)*iDw

      ! If stored
      if (Norma%VRAM) then

        ! Get profile
        eps = feps*Norma%p(if0:if1)
        eta = feta*Norma%p(if0:if1)

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        aul = Atom%ldamp(itran,iz)

        ! Level quantities

        ! Damping parameter
        au = Atom%damp(itermu,iz)

        ! Energy
        eu = Atom%FSfreq(iJu,itermu)

        ! Damping parameter
        al = Atom%damp(iterml,iz)

        ! Energy
        el = Atom%FSfreq(iJl,iterml)

        ! Intermediate quantities
        at = (au+al+aul)*iDw
        Dfreqw = (eu - el)*iDw
        feps = feps*Norma%Norm(1)
        feta = feta*Norma%Norm(1)

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*iDw,at,prof)

          ! Save contribution
          eps(ifreq) = feps*prof
          eta(ifreq) = feta*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine sourI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the second order emissivity of a given atomic line.
      !! This subroutine only computes the positive contribution of
      !! the coherent scattering and the negative contribution for the
      !! flat-spectrum, i.e., the result needs to be added to the
      !! product of emissI to get the actual emissivity\n
      !!         Atom(Atom_class): Structure with atomic data\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!               vx(double): Velocity vector along X\n
      !!               vy(double): Velocity vector along Y\n
      !!               vz(double): Velocity vector along Z\n
      !!            lvel(logical): If dynamic node\n
      !!         omega(double(:)): Frequency array\n
      !!          Fed(Reda_class): Structure with redistribution
      !!                           output frequency data\n
      !!          Red(Redb_class): Structure with redistribution input
      !!                           frequency data\n
      !!       RWarr(Redb2_class): Structure with redistribution
      !!                           function data\n
      !!        Norma(Prof_class): Normalization factors for Voigt
      !!                           profiles or Voigt profiles\n
      !!           njdir(integer): Number of output directions to
      !!                           compute\n
      !!           jtran(integer): Output transition index term wise\n
      !!          fjtran(integer): Output transition index level
      !!                           wise\n
      !!          itermu(integer): Upper term of output transition\n
      !!          itermf(integer): Lower term of output transition\n
      !!             iJu(integer): Upper level of output transition\n
      !!             iJf(integer): Lower level of output transition\n
      !!              iz(integer): Height index\n
      !!             if0(integer): First frequency index for this
      !!                           transition\n
      !!             if1(integer): Last frequency index for this
      !!                           transition\n
      !!            Mif0(integer): First frequency for this CPU\n
      !!            Mif1(integer): Last frequency for this CPU\n
      !!              DwT(double): Thermal part of Doppler width\n
      !!               Dw(double): Doppler width of the output
      !!                           transition\n
      !!              vmi(double): Microturbulent velocity\n
      !!    Stokes(double(:,:,:)): Intensity\n
      !!     JKQ(dcomplex(:,:,:)): Mean intensity integrated over the
      !!                           absorption profile\n
      !!    JKQC(dcomplex(:,:,:)): Mean intensity with frequency
      !!                           dependence\n
      !!          eps2(double(:)): Intensity emissivity\n
      !!         rpf(double(:,:)): Factor for lambda operator
      subroutine emissI2ord(Atom,Geom,vx,vy,vz,lvel,omega,Fed,Red, &
                            RWarr,Norma,njdir,jtran,fjtran,itermu, &
                            itermf,iJu,iJf,iz,if0,if1,Mif0,Mif1, &
                            DwT,Dw,vmi,Stokes,JKQ,JKQC,eps2,rpf)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Geometry_class), intent(in):: Geom
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      type(Redb2_class), intent(inout):: RWarr
      type(Prof_class), intent(in):: Norma
      logical, intent(in):: lvel
      integer, intent(in):: jtran,fjtran,itermu,itermf,iJu,iJf
      integer, intent(in):: njdir,iz,if0,if1,Mif0,Mif1
      double precision, intent(in):: DwT,Dw,vmi,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: JKQ
      double precision, dimension(Red%ggf0:Red%ggf1), &
                        intent(in), target:: JKQC
      double precision, dimension(:,:,:), intent(in):: Stokes
      double precision, dimension(njdir,if0:if1), intent(out):: eps2
      double precision, dimension(njdir,if0:if1), intent(out):: rpf

      ! Local

      logical:: PRDc,LIRAM

      integer:: i,j,iR,itran,fitran,iterml,ffjtran,ffitran,ffktran
      integer:: ish1,ith1,iph1,jdir,ifreq,iifreq,iran,nfs
      integer:: jjfreq,jjfreq0,kkfreq0,kkfreq0b,nmfreq,iJl,iti

      double precision:: auxb,norme2,rLl,rLu,rLf,S,rJl,rJu,rJf
      double precision:: el,eu,ef,al,au,af,auf,aul,Dw1
      double precision:: StokesM,Jrad,Dfreq,Dfreq1,iDw
      double precision:: PRDin,prof,hanleden,rhoc,rhou
      double precision:: atfw,atf,atl,omegai
      double precision:: cost,sint,cosc,sinc,vfac1
      double precision, dimension(if0:if1):: CRD
      double precision, dimension(njdir,if0:if1):: PRD
      double precision, dimension(:), allocatable:: Jin
      double precision, dimension(:), allocatable:: Stokesin
      double precision, dimension(:), allocatable, target:: JKQinMV
      double precision, dimension(:), allocatable, target:: Warr2

      ! Pointers

      type(Redc_class), pointer:: p_red
      type(Redc2_class), pointer:: p_rwarr
      type(Redc2_class), target:: p_dummy
      integer, pointer:: p_mfreq
      double precision, dimension(:), pointer:: p_JKQ
      double precision, dimension(:), pointer:: p_warr2


      ! Routine name
      urou = 'emissI2ord'

      ! Initialize pointers
      nullify(p_red)
      nullify(p_rwarr)
      nullify(p_mfreq)
      nullify(p_JKQ)
      nullify(p_warr2)

      ! Construct mean intensity if angle-averaged and dynamic
      if (AVI.and.lvel) &
        call getJMV(Geom,Red,omega,vx,vy,vz,DwT, &
                    Atom%ntran,Atom%tif0,Atom%tif1, &
                    Stokes,JKQinMV)

      !
      ! Get terms and transition quantities
      !

      ! Damping parameters
      au = Atom%damp(itermu,iz)
      af = Atom%damp(itermf,iz)
      auf = Atom%ldamp(jtran,iz)
      atf = au + af + auf

      ! Spin
      S = Atom%Sval(itermu)

      ! Orbital angular momentum
      rLu = Atom%rLval(itermu)
      rLf = Atom%rLval(itermf)

      ! Get FS global index of output transition
      ffjtran = Atom%ifst_ij(fjtran,jtran)

      ! Energies
      ef = Atom%FSfreq(iJf,itermf)
      eu = Atom%FSfreq(iJu,itermu)

      ! Angular momentums
      rJf = Atom%rJval(iJf,itermf)
      rJu = Atom%rJval(iJu,itermu)

      ! Transition frequency
      Dfreq = eu - ef

      ! Upper level index
      i = Atom%irho(itermu)%irho_ij(iJu)

      ! Upper level density
      rhou = Atom%popu(i,iz)

      ! Trano index
      ffktran = Atom%itrano(ffjtran)

      !
      ! Initializations
      !
      eps2 = 0d0

      ! Inverse Doppler width
      iDw = 1d0/Dw


      !
      ! Flat contribution. Implicit branching
      !

      ! If stored
      if (Norma%VRAM) then

        CRD = Norma%p(if0:if1)*1d5*sqrt(IPI)*iDw

      ! If not stored
      else

        ! Scale with Doppler width and get norm
        atfw = atf*iDw
        norme2 = Norma%Norm(1)*1d5*sqrt(IPI)*iDw

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile u-f
          call voigtI((Dfreq - omega(ifreq))*iDw,atfw,prof)

          ! Flat spectrum contribution
          CRD(ifreq) = prof*Norme2

        end do ! frequencies

      end if ! Storing Voigt


      !
      ! Calculation of 2nd order emissivity
      !

      ! For all the possible low->up transitions
      do iti=1,Atom%tranoI(ffktran)%nt

        ! Transition indexes
        ffitran = Atom%tranoI(ffktran)%indT(iti)
        itran = Atom%ifst(ffitran)
        fitran = Atom%ifstj(ffitran)

        ! Term and level indexes
        iterml = Atom%fst(itran)%iterml
        iJl = Atom%fst(itran)%ilevell(fitran)

        ! Damping parameter lower level input transition
        al = Atom%damp(iterml,iz)
        aul = Atom%ldamp(itran,iz)
        atl = au + al + aul

        ! Angular momentum input lower level
        rLl = Atom%rLval(iterml)

        ! Energy input lower level
        el = Atom%FSfreq(iJl,iterml)

        ! Angular momentum
        rJl = Atom%rJval(iJl,iterml)

        ! Point to input transition
        p_red => Red%trani(iti)

        ! If IRAM, point to the redistribution subblock
        if (IRAM) then

          p_rwarr => RWarr%trani(iti)
          LIRAM = IRAM.and.p_rwarr%RAM

        ! If not, nothing stored
        else

          p_rwarr => p_dummy
          LIRAM = .False.

        end if ! Storing redistribution

        ! Get frequency size
        nmfreq = sum(p_red%mfreq)

        ! If AA
        if (AVI) then

          ! If dynamic
          if (lvel) then

            ! Get input radiation field
            call getJin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega, &
                        Jin,JKQinMV)

          ! If static
          else

            ! Get input radiation field
            call getJin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega,Jin,JKQC)

          end if

        ! If AD
        else

          ! If Rayleigh scattering and there is coherent
          if (ffjtran.eq.ffitran.and.Geom%V_CScatt(1).ge.1d0) then

            ! Scale dimension
            nmfreq = nmfreq*(Geom%nScatt-1)
            nfs = 1

          ! Raman scattering
          else

            ! Scale dimension
            nmfreq = nmfreq*Geom%nScatt
            nfs = 0

          end if ! Rayleigh/Raman scattering

          ! Get interpolated intensity
          call getStkinI(Geom,p_red,Fed,Red,Mif0,Mif1,omega, &
                         vx,vy,vz,lvel,Stokesin,Stokes)

        end if

        ! Flat spectrum J00
        JRad = JKQ(ffitran)

        ! Input transition frequency
        Dfreq1 = eu - el

        ! Doppler width for the input transition
        Dw1 = Dfreq1*sqrt(DwT*DwT + vmi**2d0)

        ! Hanle factor
        ! TODO ATTENTION TO THIS
        hanleden = 2d0*(au+auf)*iDw

        ! If storing redistribution function
        if (LIRAM) then

          ! Check we need to initialize
          PRDc = p_rwarr%iIPRD

        ! Not storing
        else

          ! Always calculate
          PRDc = .True.

        end if ! Storing redistribution function

        ! If need to compute it and there are frequencies
        if (PRDc.and.nmfreq.gt.0) then

          ! Calculate redistribution function Warr2
          call get_WarrI(Geom,Fed,p_red,p_rwarr,LIRAM,Mif0,Mif1, &
                         ffitran,ffjtran,nmfreq,omega,Dw,Dw1, &
                         Dfreq,Dfreq1,atl,atf,Warr2)

        end if ! initialized

        ! Lower input level index
        j = Atom%irho(iterml)%irho_ij(iJl)

        ! Lower input level SEE index
        iR = Atom%irho(iterml)%Jrho(iJl,iJl)%kq(0,0)
        rhoc = sqrt(2d0*rJl+1d0)*dble(Atom%crho(iR,iz))


        !
        ! Integral over input frequencies
        !

        ! If storing Warr
        if (LIRAM) then

          ! If there are frequencies
          if (nmfreq.gt.0) then

            ! Copy from RAM
            allocate(p_warr2(nmfreq))
            p_warr2 = dble(p_rwarr%IWarr2)

          end if ! There are frequencies

        ! If not storing
        else

          ! Just point to the one calculated here
          if (allocated(Warr2)) p_warr2 => Warr2

        end if ! Storing Warr


        !
        ! Angle-average (Integral)
        !
        if (AVI) then

          ! If dynamic
          if (lvel) then

            p_JKQ(Red%ggf0:Red%ggf1) => JKQinMV

          ! If static
          else

            p_JKQ(Red%ggf0:Red%ggf1) => JKQC

          end if

          ! Initialize frequency index
          jjfreq = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! Manage MPI
              if (iifreq.lt.Mif0) cycle
              if (iifreq.gt.Mif1) exit

              ! Point to dimension
              p_mfreq => p_red%mfreq(iifreq)

              ! If coherent wing
              if (p_mfreq.lt.1) then

                ! Fully coherent contribution
                PRD(1,ifreq) = p_JKQ(ifreq) - Jrad

                ! Skip rest
                cycle

              end if

              ! Initialize
              PRD(1,ifreq) = sum(p_warr2(jjfreq+1:jjfreq+p_mfreq)* &
                                 Jin(jjfreq+1:jjfreq+p_mfreq))

              ! Subtract the flat spectrum part due to just
              ! radiative excitation
              PRD(1,ifreq) = PRD(1,ifreq) - Jrad

              ! Update jjfreq
              jjfreq = jjfreq + p_mfreq

            end do ! output frequencies
          end do ! output frequencies ranges

          ! Clean radiation field
          if (allocated(Jin)) deallocate(Jin)

        !
        ! Angle-dependent (Integral)
        !
        else

          ! Initialize
          PRD = 0d0

          ! If axial symmetry
          if (axiali) then

            ! For each output direction
            do jdir=1,Geom%njdir

              ! Initialize frequency indexes
              jjfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fed%nran
                do ifreq=Fed%if0(iran),Fed%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Manage MPI
                  if (iifreq.lt.Mif0) cycle
                  if (iifreq.gt.Mif1) exit

                  ! Initialize
                  PRDin = 0d0

                  ! Point to dimension
                  p_mfreq => p_red%mfreq(iifreq)

                  ! For each polar direction (input)
                  do ith1=1,Geom%nTh

                    ! For each azimuthal direction (input)
                    do iph1=1,Geom%nPh2

                      ! Scattering index
                      ish1 = Geom%i_scatt(iph1,ith1,jdir)

                      ! Special treatment if forward for two terms
                      if ((ffjtran.eq.ffitran.and. &
                           Geom%V_CScatt(ish1).ge.1d0).or. &
                          (p_mfreq.lt.1)) then

                        ! If there are dynamics
                        if (lvel) then

                          ! Get cosine of direction
                          cost = Geom%V_mu(ith1)

                          ! Calculate Doppler shift factor
                          vfac1 = 1d0 - vz*cost

                          ! We will be using the inverse
                          vfac1 = 1d0/vfac1

                          ! Input frequency
                          omegai = omega(ifreq)*vfac1

                          ! Get Stokes
                          StokesM = getStkinInu(omega, &
                                                Stokes(:,1,ith1), &
                                                ifreq, &
                                                Atom%tif0(itran), &
                                                Atom%tif1(itran), &
                                                omegai)
                        ! Static
                        else

                          StokesM = Stokes(ifreq,1,ith1)

                        end if ! Dynamics

                        ! Add the directional weights
                        PRDin = PRDin + StokesM*Geom%W_mu(ith1)* &
                                        Geom%W_mux2(iph1)
                      else

                        !
                        ! Find initial index for kkfreq

                        ! Shift in indexes
                        kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                        ! Add contribution
                        PRDin = PRDin + &
                                Geom%W_mu(ith1)* &
                                Geom%W_mux2(iph1)* &
                        sum(Stokesin(jjfreq0+1:jjfreq0+p_mfreq)* &
                            p_warr2(kkfreq0b+1: &
                                    kkfreq0b+p_mfreq))

                      end if ! Forward scattering two terms

                    end do ! azimuthal nodes

                    ! Update jjfreq
                    jjfreq0 = jjfreq0 + p_mfreq

                  end do ! polar nodes

                  ! Update kkfreq0
                  kkfreq0 = kkfreq0 + p_mfreq*(Geom%nScatt-nfs)

                  ! Subtract the flat spectrum part due to
                  ! just radiative excitation
                  PRD(jdir,ifreq) = PRDin - Jrad

                end do ! output frequencies
              end do ! output frequencies ranges
            end do ! Output directions

          ! If non-axial symmetryc
          else

            ! For each output direction
            do jdir=1,Geom%njdir

              ! Initialize indexes
              jjfreq0 = 0
              kkfreq0 = 0

              ! For each output frequency
              iifreq = 0
              do iran=1,Fed%nran
                do ifreq=Fed%if0(iran),Fed%if1(iran)

                  ! Advance index
                  iifreq = iifreq + 1

                  ! Manage MPI
                  if (iifreq.lt.Mif0) cycle
                  if (iifreq.gt.Mif1) exit

                  ! Initialize
                  PRDin = 0d0

                  ! Point to dimension
                  p_mfreq => p_red%mfreq(iifreq)

                  ! For each polar direction
                  do ith1=1,Geom%nTh

                    ! For each azimuthal direction
                    do iph1=1,Geom%nPh2

                      ! Scattering index
                      ish1 = Geom%i_scatt(iph1,ith1,jdir)

                      ! Special treatment if forward for two-terms
                      if ((ffjtran.eq.ffitran.and. &
                           Geom%V_CScatt(ish1).ge.1d0).or. &
                          (p_mfreq.lt.1)) then

                        ! If there are dynamics
                        if (lvel) then

                          cost = Geom%V_mu(ith1)
                          sint = sqrt(1d0 - cost*cost)
                          cosc = Geom%v_mux(iph1)
                          sinc = Geom%v_muy(iph1)* &
                                 sqrt(1d0 - cosc*cosc)

                          ! Calculate Doppler shift factor
                          vfac1 = 1d0 - vx*sint*cosc - &
                                        vy*sint*sinc - &
                                        vz*cost

                          ! We will be using the inverse
                          vfac1 = 1d0/vfac1

                          ! Input frequency
                          omegai = omega(ifreq)*vfac1

                          StokesM = getStkinInu(omega, &
                                                Stokes(:,iph1,ith1), &
                                                ifreq, &
                                                Atom%tif0(itran), &
                                                Atom%tif1(itran), &
                                                omegai)
                        ! Static
                        else

                          StokesM = Stokes(ifreq,iph1,ith1)

                        end if ! Dynamics

                        PRDin = PRDin + &
                                StokesM*Geom%W_mu(ith1)* &
                                Geom%W_mux2(iph1)

                      ! No forward two-terms
                      else

                        ! Shift in indexes
                        kkfreq0b = kkfreq0 + (ish1-nfs-1)*p_mfreq

                        ! Add contribution
                        PRDin = PRDin + &
                                Geom%W_mu(ith1)* &
                                Geom%W_mux2(iph1)* &
                        sum(Stokesin(jjfreq0+1:jjfreq0+p_mfreq)* &
                            p_warr2(kkfreq0b+1:kkfreq0b+p_mfreq))

                      end if ! Forward scattering two terms

                      ! Update indexes
                      jjfreq0 = jjfreq0 + p_mfreq

                    end do ! azimuthal nodes
                  end do ! polar nodes

                  ! Update kkfreq0
                  kkfreq0 = kkfreq0 + p_mfreq*(Geom%nScatt-nfs)

                  ! Subtract the flat spectrum part due to
                  ! just radiative excitation
                  PRD(jdir,ifreq) = PRDin - Jrad

                end do ! output frequencies
              end do ! output frequencies ranges
            end do ! Output directions

          end if ! Axial symmetry

          ! Free Stokes
          if (allocated(Stokesin)) deallocate(Stokesin)

        end if ! AA/AD

        ! Clean warr2
        if (LIRAM.and.nmfreq.gt.0) deallocate(p_warr2)
        nullify(p_warr2)

        ! Apply hanle factor and Einstein coefficient
        auxb = rhoc*Atom%fst(itran)%Blu(iJl,iJu)/hanleden
        eps2 = eps2 + PRD*auxb

      end do ! Input transitions

      ! Calculate factor for ALI and eps2
      auxb = 1d-8*IPI2*iDw/c
      rpF = eps2*auxb
      auxb = 1d-2*IPI4*Atom%fst(jtran)%Aul(iJu,iJf)
      do jdir=1,njdir
        eps2(jdir,:) = rpF(jdir,:)*CRD*auxb
      end do
      rpF = rpF/rhou + 1d0

      ! Clean pointers/arrays
      if (associated(p_red)) nullify(p_red)
      if (associated(p_rwarr)) nullify(p_rwarr)
      if (associated(p_mfreq)) nullify(p_mfreq)
      if (associated(p_JKQ)) nullify(p_JKQ)
      if (associated(p_warr2)) nullify(p_warr2)
      if (allocated(JKQinMV)) deallocate(JKQinMV)

      return

      end subroutine emissI2ord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the redistribution function\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!       Fed(Reda_class): Structure with redistribution output
      !!                        frequency data\n
      !!     p_red(Redc_class): Structure with redistribution input
      !!                        frequency data for a given input
      !!                        transition\n
      !!  p_rwarr(Redc2_class): Structure to store redistribution
      !!                        functions\n
      !!        LIRAM(logical): If storing the redistribution
      !!                        function in RAM\n
      !!         Mif0(integer): First frequency for this CPU\n
      !!         Mif1(integer): Last frequency for this CPU\n
      !!      ffitran(integer): Input transition index\n
      !!      ffjtran(integer): Output transition index\n
      !!       nmfreq(integer): Total amount of input frequencies
      !!                        for this combination of output and
      !!                        input transitions in this CPU\n
      !!         omega(double): Frequency array\n
      !!            Dw(double): Doppler width of the output
      !!                        transition\n
      !!           Dw1(double): Doppler width of the input
      !!                        transition\n
      !!         Dfreq(double): Output transition energy\n
      !!        Dfreq1(double): Input transition energy\n
      !!           atl(double): Damping parameter of the input
      !!                        transition\n
      !!           atf(double): Damping parameter of the output
      !!                        transition\n
      !!      Warr2(double(:)): Redistribution function
      subroutine get_WarrI(Geom,Fed,p_red,p_rwarr,LIRAM, &
                           Mif0,Mif1,ffitran,ffjtran,nmfreq,omega, &
                           Dw,Dw1,Dfreq,Dfreq1,atl,atf,Warr2)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(Reda_class), intent(in):: Fed
      type(Redc_class), intent(in), pointer:: p_red
      type(Redc2_class), intent(inout):: p_rwarr
      logical, intent(in):: LIRAM
      integer, intent(in):: ffitran,ffjtran,nmfreq,Mif0,Mif1
      double precision, intent(in):: Dfreq,Dfreq1,atl,atf,Dw,Dw1
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), &
                        allocatable, intent(inout):: Warr2

      ! Local

      integer:: jjfreq0,kkfreq0,jjfreq,kkfreq,iifreq,iran,ifreq,jfreq
      integer:: ith1,ish1,stype

      double precision:: Norme2,omegao,omegai

      ! Pointers

      integer, pointer:: p_mfreq


      ! If Warr2 not allocated
      if (.not.allocated(Warr2)) then

        ! Allocate
        allocate(Warr2(nmfreq))

      ! If already allocated
      else

        ! If wrong size
        if (size(Warr2).ne.nmfreq) then

          ! Re-allocate
          deallocate(Warr2)
          allocate(Warr2(nmfreq))

        end if ! Wrong-size
      end if ! Allocated or not

      ! AA
      if (AVI) then

        ! Inititlize redistribution
        Warr2 = 0d0

        ! Initialize frequency index
        jjfreq = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fed%nran
          do ifreq=Fed%if0(iran),Fed%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Manage MPI
            if (iifreq.lt.Mif0) cycle
            if (iifreq.gt.Mif1) exit

            ! Point to dimension
            p_mfreq => p_red%mfreq(iifreq)

            ! Skip coherent
            if (p_mfreq.lt.1) cycle

            ! Get output frequency
            omegao = omega(ifreq) - Dfreq

            ! Initialize norm
            Norme2 = 0d0

            ! For each input frequency
            do kkfreq=jjfreq+1,jjfreq+p_mfreq

              ! Get input frequency
              omegai = p_red%omega(kkfreq) - Dfreq1

              ! For each direction in the integral AA quadrature
              do ith1=1,Geom%nThAA

                ! Add the contribution to the angular integral
                ! of the redistribution function
                Warr2(kkfreq) = Warr2(kkfreq) + &
                                Geom%W_muAA(ith1)* &
                        WfuncI(omegai,omegao,Dw,Dw1,atl,atf, &
                               Geom%V_muAA(ith1), &
                               Geom%V_siAA(ith1),0)*IPI42

              end do  ! direction nodes

              ! Store with weight
              Warr2(kkfreq) = Warr2(kkfreq)* &
                              p_red%w_freq(kkfreq)

              ! Add to norm
              Norme2 = Norme2 + Warr2(kkfreq)

            end do ! input frequencies

            ! Update jjfreq
            jjfreq = jjfreq + p_mfreq

            ! Apply normalization
            if (Norme2.gt.0d0) &
              Warr2(jjfreq-p_mfreq+1:jjfreq) = &
                              Warr2(jjfreq-p_mfreq+1:jjfreq)/ &
                              Norme2

          end do ! output frequencies
        end do ! output frequencies ranges

      ! AD
      else

        ! Initialize frequency indexes
        jjfreq0 = 0
        kkfreq0 = 0

        ! For each output frequency
        iifreq = 0
        do iran=1,Fed%nran
          do ifreq=Fed%if0(iran),Fed%if1(iran)

            ! Advance index
            iifreq = iifreq + 1

            ! Manage MPI
            if (iifreq.lt.Mif0) cycle
            if (iifreq.gt.Mif1) exit

            ! Point to dimension
            p_mfreq => p_red%mfreq(iifreq)

            ! Coherent wing
            if (p_mfreq.lt.1) cycle

            ! Get output frequency
            omegao = omega(ifreq) - Dfreq

            ! For each scattering angle
            do ish1=1,Geom%nScatt

              ! Check forward scattering two-terms
              if (ffitran.eq.ffjtran.and. &
                  Geom%V_CScatt(ish1).ge.1d0) cycle

              ! Check forward or backward
              if (Geom%V_SScatt(ish1).le.0d0) then
                stype = 1
              else
                stype = 0
              end if

              ! Inititlize Norm
              Norme2 = 0d0

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq
                kkfreq = kkfreq0 + jfreq

                ! Get input frequency
                omegai = p_red%omega(jjfreq) - Dfreq1

                ! Calculate redistribution function
                ! and apply weight
                Warr2(kkfreq) = WfuncI(omegai,omegao, &
                                       Dw,Dw1,atl,atf, &
                                       Geom%V_CScatt(ish1), &
                                       Geom%V_SScatt(ish1), &
                                       stype)* &
                                IPI42*p_red%w_freq(jjfreq)

                ! Add to norm
                Norme2 = Norme2 + Warr2(kkfreq)

              end do ! mfreq

              ! Normalize
              if (Norme2.gt.0d0) &
                Warr2(kkfreq-p_mfreq+1:kkfreq) = &
                           Warr2(kkfreq-p_mfreq+1:kkfreq)/Norme2

              ! Update kkfreq0
              kkfreq0 = kkfreq0 + p_mfreq

            end do ! Scattering angle

            ! Update jjfreq0
            jjfreq0 = jjfreq0 + p_mfreq

          end do ! Output frequencies
        end do  ! Output frequency ranges

      end if ! AA/AD

      ! If storing
      if (LIRAM) then

        ! Signal that this one does not need to be calculated
        ! anymore
        p_rwarr%iIPRD = .False.

        ! For the just calculated frequencies
        do jfreq=1,nmfreq

          ! Zero our small for single precision to avoid
          ! underflow
          if (Warr2(jfreq).gt.TINYWAR) then

            ! Get single precision
            p_rwarr%IWarr2(jfreq) = real(Warr2(jfreq))

          ! Too small
          else

            ! Make zero
            p_rwarr%IWarr2(jfreq) = 0e0

          end if ! Single precision value

        end do ! Calculated frequencies

      end if ! If storing

      end subroutine get_WarrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity and emissivity for the intensity
      !! for a given LTE line\n
      !!  line(LTEline_class): Structure with LTE line data\n
      !!     omega(double(:)): Frequency array\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!    aprof(Prof_class): Voigt profiles\n
      !!           Dw(double): Doppler width of the transition\n
      !!         vfac(double): Doppler shift factor\n
      !!           pE(double): Unit transformation factor\n
      !!       eta(double(:)): Intensity absorptivity\n
      !!       eps(double(:)): Intensity emissivity
      subroutine rt1ordILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,pE, &
                            eta,eps)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Prof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw,pE,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq

      double precision:: iDw,at,feta,feps,Dfreqw,vfacw,prof


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE*iDw

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul*iDw

      ! If stored in RAM
      if (aprof%VRAM) then

        ! Get stored profiles
        eta = feta*aprof%p
        eps = feps*aprof%p

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        at = line%damp(iz)*iDw

        ! Energy
        Dfreqw = (line%Eu - line%El)*iDw

        ! Shift
        vfacw = vfac*iDw

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          ! Save contribution
          eta(ifreq) = feta*prof
          eps(ifreq) = feps*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine rt1ordILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity for the intensity for a given LTE
      !! line\n
      !!  line(LTEline_class): Structure with LTE line data\n
      !!     omega(double(:)): Frequency array\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!    aprof(Prof_class): Voigt profiles\n
      !!           Dw(double): Doppler width of the transition\n
      !!         vfac(double): Doppler shift factor\n
      !!           pE(double): Unit transformation factor\n
      !!       eta(double(:)): Intensity absorptivity
      subroutine absorbILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,pE, &
                            eta)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Prof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw,pE,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq

      double precision:: iDw,at,feta,Dfreqw,vfacw,prof


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      ! Absorptibity factor
      feta = line%nl(iz)*1d3*IPI41*line%Blu*pE*iDw

      ! If stored in RAM
      if (aprof%VRAM) then

        ! Get profile
        eta = feta*aprof%p

      ! Not stored
      else

        ! Level quantities

        ! Damping parameter
        at = line%damp(iz)*iDw

        ! Energy
        Dfreqw = (line%Eu - line%El)*iDw

        ! Shift
        vfacw = vfac*iDw

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          ! Save contribution
          eta(ifreq) = feta*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine absorbILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity for the intensity for a given LTE
      !! line\n
      !!  line(LTEline_class): Structure with LTE line data\n
      !!     omega(double(:)): Frequency array\n
      !!          iz(integer): Height index\n
      !!         if0(integer): First frequency index for this
      !!                       transition\n
      !!         if1(integer): Last frequency index for this
      !!                       transition\n
      !!    aprof(Prof_class): Voigt profiles\n
      !!           Dw(double): Doppler width of the transition\n
      !!         vfac(double): Doppler shift factor\n
      !!       eps(double(:)): Intensity emissivity
      subroutine emissILTE(line,omega,iz,if0,if1,aprof,Dw,vfac,eps)

      ! I/O

      type(LTEline_class), intent(in):: line
      type(Prof_class), intent(in):: aprof
      integer, intent(in):: iz,if0,if1
      double precision, intent(in):: Dw,vfac
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eps

      ! Local

      integer:: ifreq

      double precision:: iDw,at,Dfreqw,vfacw,prof,feps


      ! Inverse Doppler width
      iDw = 1d0/Dw

      !
      ! Get population factor
      !

      ! Emissivity factor
      feps = line%nu(iz)*1d3*IPI41*line%Aul*iDw

      ! If stored in RAM
      if (aprof%VRAM) then

        ! Get profile
        eps = feps*aprof%p

      ! Not stored
      else

        ! Transition quantities

        ! Damping parameter
        at = line%damp(iz)*iDw

        ! Energy
        Dfreqw = (line%Eu - line%El)*iDw

        ! Shift
        vfacw = vfac*iDw

        ! For each frequency
        do ifreq=if0,if1

          ! Calculate profile
          call voigtI(Dfreqw - omega(ifreq)*vfacw,at,prof)

          ! Save contribution
          eps(ifreq) = feps*prof

        end do ! frequencies

      end if ! Type of profile calculation

      return

      end subroutine emissILTE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the absorptivity of a given photoionization
      !! transition\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!    itran(integer): Transition index\n
      !!  ilevell(integer): Lower level index\n
      !!       iz(integer): Height index\n
      !!      if0(integer): First frequency index for this
      !!                    transition\n
      !!      if1(integer): Last frequency index for this
      !!                    transition\n
      !!    eta(double(:)): Absorptivity
      subroutine photoabsI(Atom,itran,ilevell,iz,if0,if1,eta)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran,ilevell,iz,if0,if1
      double precision, dimension(if0:if1), intent(out):: eta

      ! Local

      integer:: ifreq,iterml,iJl,iR

      double precision:: rJl,rhol


      !
      ! Get indexes
      !

      ! Term index
      iterml = Atom%term(ilevell)

      ! J level index
      iJl = Atom%sublevel(ilevell)

      ! Angular momentum
      rJl = Atom%rJval(iJl,iterml)

      ! SEE index
      iR = Atom%irho(iterml)%Jrho(iJl,iJl)%kq(0,0)

      ! Population lower level
      rhol = sqrt(2d0*rJl+1d0)*dble(Atom%crho(iR,iz))

      ! For each frequency
      do ifreq=if0,if1

        ! Compute absorptivity
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*rhol

      end do ! frequencies

      return

      end subroutine photoabsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given recombination
      !! transition\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  omega(double(:)): Frequency array\n
      !!         T(double): Temperature\n
      !!        ne(double): Electron number density\n
      !!    itran(integer): Transition index\n
      !!  ilevelu(integer): Upper level index\n
      !!       iz(integer): Height index\n
      !!      if0(integer): First frequency index for this
      !!                    transition\n
      !!      if1(integer): Last frequency index for this
      !!                    transition\n
      !!    eps(double(:)): Emissivity\n
      !!    eta(double(:)): Stimulated emissivity
      !!      rhou(double): Factor for Lambda operator
      subroutine photoepsI(Atom,omega,T,ne,itran,ilevelu,iz,if0,if1, &
                           eps,eta,rhou)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran,ilevelu,iz,if0,if1
      double precision, intent(in):: T, ne
      double precision, intent(out):: rhou
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq, iJu, itermu, iR

      double precision:: c0,c1,exu,pE,Saha,rJu,tmp
      double precision, dimension(if0:if1):: omega3


      !
      ! Saha term constants
      !
      c0 = fktoJ/kb/T
      Saha = cSaha*ne*Atom%phot(itran)%glu* &
             exp(Atom%phot(itran)%edge*c0)/(T**(1.5d0))

      !
      ! Indexes
      !

      ! Get term index
      itermu = Atom%term(ilevelu)

      ! Get J level index
      iJu = Atom%sublevel(ilevelu)

      ! Get angular momentum
      rJu = Atom%rJval(iJu,itermu)

      ! Get SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! Get upper level rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Get upper level population
      tmp = sqrt(2d0*rJu+1d0)*rhou

      ! Apply Saha factor
      tmp = tmp*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! Compute numerator of rpf
      Saha = 1d0/Saha/sqrt(2d0*rJu+1d0)

      ! Get cubic frequency
      omega3 = omega(if0:if1)
      omega3 = omega3*omega3*omega3

      ! For each frequency
      do ifreq=if0,if1

        ! Exponential
        exu = c0*omega(ifreq)
        exu = diexp(exu)

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu*tmp

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies

      return

      end subroutine photoepsI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the emissivity of a given recombination transition
      !! with frequency quantities stored in RAM\n
      !!   Atom(Atom_class): Structure with the atomic data\n
      !!  omega3(double(:)): Frequency array to the power of 3\n
      !!     exu(double(:)): Exponential in the stimulated part of the
      !!                     photoionization emissivity\n
      !!          T(double): Temperature\n
      !!         ne(double): Electron number density\n
      !!     itran(integer): Transition index\n
      !!   ilevelu(integer): Upper level index\n
      !!        iz(integer): Height index\n
      !!       if0(integer): First frequency index for this
      !!                     transition\n
      !!       if1(integer): Last frequency index for this
      !!                     transition\n
      !!     eps(double(:)): Emissivity\n
      !!     eta(double(:)): Stimulated emissivity\n
      !!       rhou(double): Factor for Lambda operator
      subroutine photoepsIS(Atom,omega3,exu,T,ne,itran, &
                            ilevelu,iz,if0,if1,eps,eta,rhou)

      ! I/O

      type(Atom_class), intent(in):: Atom
      integer, intent(in):: itran,ilevelu,iz,if0,if1
      double precision, intent(in):: T, ne
      double precision, intent(out):: rhou
      double precision, dimension(if0:if1), intent(in):: omega3
      double precision, dimension(if0:if1), intent(in):: exu
      double precision, dimension(if0:if1), intent(out):: eta,eps

      ! Local

      integer:: ifreq,iJu,itermu,iR

      double precision:: c0,c1,pE,Saha,rJu,tmp


      !
      ! Saha term constants
      !
      c0 = fktoJ/kb/T
      Saha = cSaha*ne*Atom%phot(itran)%glu* &
             exp(Atom%phot(itran)%edge*c0)/(T**(1.5d0))

      !
      ! Indexes
      !

      ! Get term index
      itermu = Atom%term(ilevelu)

      ! Get J level index
      iJu = Atom%sublevel(ilevelu)

      ! Get angular momentum
      rJu = Atom%rJval(iJu,itermu)

      ! Get SEE index
      iR = Atom%irho(itermu)%Jrho(iJu,iJu)%kq(0,0)

      ! Get upper level rho00
      rhou = dble(Atom%crho(iR,iz))

      ! Get upper level population
      tmp = sqrt(2d0*rJu+1d0)*rhou

      ! Apply Saha factor
      tmp = tmp*Saha

      ! Compute exponential argument constant
      c0 = c2*1d4/T

      ! Compute energy constant part
      c1 = 2d21*c*convF

      ! Compute numerator of rpf
      Saha = 1d0/Saha/sqrt(2d0*rJu+1d0)

      ! For each frequency
      do ifreq=if0,if1

        ! Compute energy part
        pE = c1*omega3(ifreq)

        ! Stimulated part
        eta(ifreq) = Atom%phot(itran)%alpha(ifreq)*exu(ifreq)*tmp

        ! Emissivity
        eps(ifreq) = eta(ifreq)*pE

      end do ! frequencies

      return

      end subroutine photoepsIS

!#####################################################################
!#####################################################################
!#####################################################################

      !> Calculate the frequency dependent mean intensity for the
      !! angle-averaged second order emissivity in the presence of
      !! velocities in the comoving frame\n
      !!   Geom(Geometry_class): Structure with geometric data\n
      !!        Red(Redb_class): Structure with redistribution input
      !!                         frequency data\n
      !!       omega(double(:)): Frequency array\n
      !!             vx(double): Velocity vector along X\n
      !!             vy(double): Velocity vector along Y\n
      !!             vz(double): Velocity vector along Z\n
      !!            DwT(double): Thermal part of the Doppler width\n
      !!         ntran(integer): Number of transitions in the atom\n
      !!       tif0(integer(:)): Lower limit to search in Stokes
      !!                         interpolation\n
      !!       tif1(integer(:)): Upper limit to search in Stokes
      !!                         interpolation\n
      !!  Stokes(double(:,:,:)): Original intensity\n
      !!     JKQinMV(double(:)): Frequency dependent JKQ in the mean
      !!                         intensity in the comoving frame
      subroutine getJMV(Geom,Red,omega,vx,vy,vz,DwT,ntran,tif0,tif1, &
                        Stokes,JKQinMV)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(Redb_class), intent(in):: Red
      integer, intent(in):: ntran
      integer, dimension(:),intent(in):: tif0,tif1
      double precision, intent(in):: DwT,vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), &
                        allocatable, intent(out):: JKQinMV
      double precision, dimension(:,:,:), intent(in):: Stokes

      ! Local

      logical:: shift

      integer:: ifreq,ith1,iph1,if0,if1,jtran

      double precision:: vfac1,omegao,StokesM,cost,sint,cosc,sinc


      ! Allocate mean intensity
      allocate(JKQinMV(Red%ggf0:Red%ggf1))

      ! If velocity is below threshold
      shift = (vx*vx + vy*vy + vz*vz)*1d6*c.ge. &
               vrfrac*DwT
      vfac1 = 1d0


      ! Initialize
      JKQinMV = 0d0

      ! For each frequency, get mean intensity
      do ifreq=Red%ggf0,Red%ggf1

        ! Initialize
        if0 = 1
        if1 = nfreq

        ! Find transition for frequency limits
        do jtran=1,ntran

          ! If out of limits, skip
          if (ifreq.lt.tif0(jtran)) cycle
          if (ifreq.gt.tif1(jtran)) cycle

          ! Found
          if0 = tif0(jtran)
          if1 = tif1(jtran)
          exit

        end do

        ! For each polar direction
        do ith1=1,Geom%nTh

          ! Get director cosines
          if (shift) cost = Geom%V_mu(ith1)

          ! If axial symmetric
          if (axiali) then

            ! If there is Doppler shift
            if (shift) then

              ! Compute shift
              vfac1 = 1d0 + vz*cost

              ! Target frequency
              omegao = omega(ifreq)*vfac1

              ! Interpolate
              StokesM = getStkinInu(omega, &
                                    Stokes(:,1,ith1), &
                                    ifreq,if0,if1,omegao)

              ! Add to integral
              JKQinMV(ifreq) = JKQinMV(ifreq) + &
                               StokesM*Geom%W_mu(ith1)
            ! Static
            else

              ! Add to integral
              JKQinMV(ifreq) = JKQinMV(ifreq) + &
                               Stokes(ifreq,1,ith1)* &
                               Geom%W_mu(ith1)
            end if

          ! Non axial symmetric
          else

            ! Get direction cosines
            if (shift) sint = sqrt(1d0 - cost*cost)

            ! For each azimuth
            do iph1=1,Geom%nPh

              ! If there is Doppler shift
              if (shift) then

                ! Calculate Doppler shift (inverse)
                cosc = Geom%v_mux(iph1)
                sinc = Geom%v_muy(iph1)* &
                       sqrt(1d0 - cosc*cosc)
                vfac1 = 1d0 + vx*sint*cosc + &
                              vy*sint*sinc + &
                              vz*cost

                ! Get target frequency
                omegao = omega(ifreq)*vfac1

                ! Interpolate
                StokesM = getStkinInu(omega, &
                                      Stokes(:,iph1,ith1), &
                                      ifreq,if0,if1,omegao)

                ! Add to integral
                JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                 StokesM*Geom%W_mu(ith1)* &
                                 Geom%W_mux2(iph1)

              ! No shift
              else

                ! Add to integral
                JKQinMV(ifreq) = JKQinMV(ifreq) + &
                                 Stokes(ifreq,iph1,ith1)* &
                                 Geom%W_mu(ith1)* &
                                 Geom%W_mux2(iph1)

              end if ! Doppler shift

            end do ! Azimuth

          end if ! Axial symmetry

        end do ! Polar
      enddo ! Frequencies

      end subroutine getJMV

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the intensity into the requested frequency
      !!   omega(double(:)): Frequency array\n
      !!  Stokes(double(:)): Intensity\n
      !!     ifreq(integer): Frequency index of the output frequency
      !!                     associated to the requested input
      !!                     frequency\n
      !!       if0(integer): Lower limit to search in Stokes
      !!                     interpolation\n
      !!       if1(integer): Upper limit to search in Stokes
      !!                     interpolation\n
      !!          x(double): Frequency to interpolate into
      double precision function getStkinInu(omega,Stokes,ifreq, &
                                            if0,if1,x)

      ! I/O

      integer, intent(in):: ifreq,if0,if1
      double precision, intent(in):: x
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: Stokes

      ! Local

      integer:: jfreq

      double precision:: dx, dy


      ! Initialize as equals
      getStkinInu = Stokes(ifreq)

      ! If omegai > omega(ifreq)
      if (x.gt.omega(ifreq)) then

        ! If out of right boundary
        if (x.ge.omega(if1)-TINYO) then

          ! Get intensity from last point
          getStkinInu = Stokes(if1)
          return

        ! If within boundaries, look for where
        else

          ! Look after the perfect resonance
          do jfreq=ifreq,if1-1

            ! Skip before
            if (x.lt.omega(jfreq)-TINYO) cycle

            ! Between
            if (x.lt.omega(jfreq+1)) then

              ! Same
              if (abs(x - omega(jfreq)).lt.TINYO) then

                ! Get exact intensity
                getStkinInu = Stokes(jfreq)

              ! Actually between
              else

                ! Slope numerator
                dy = Stokes(jfreq+1) - Stokes(jfreq)

                ! Slope denominator
                dx = x - omega(jfreq)

                ! Linear interpolation
                getStkinInu = Stokes(jfreq) + &
                              dx*dy/(omega(jfreq+1) - omega(jfreq))

              end if ! Same or between

              ! Found, return
              return

            end if ! Found range

          end do ! Searching frequency

        end if ! Within boundaries

      ! If omegai < omegao
      else if (x.lt.omega(ifreq)) then

        ! If out of left boundary
        if (x.le.omega(if0)+TINYO) then

          ! Get intensity from first point
          getStkinInu = Stokes(if0)
          return

        ! If within boundaries, look for where
        else

          ! Look before the perfect resonance
          do jfreq=ifreq,if0+1,-1

            ! If larger, skip
            if (x.gt.omega(jfreq)+TINYO) cycle

            ! Between
            if (x.gt.omega(jfreq-1)) then

              ! Same
              if (abs(x - omega(jfreq)).lt.TINYO) then

                ! Get exact intensity
                getStkinInu = Stokes(jfreq)

              ! Actually between
              else

                ! Slope numerator
                dy = Stokes(jfreq) - Stokes(jfreq-1)

                ! Slope denominator
                dx = x - omega(jfreq-1)

                ! Linear interpolation
                getStkinInu = Stokes(jfreq-1) + &
                              dx*dy/(omega(jfreq) - omega(jfreq-1))

              end if ! Same or between

              ! Found, return
              return

            end if ! Found range

          end do ! Searching frequency

        end if ! Within boundaries
      end if ! omegai > omegao

      return

      end function getStkinInu

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the Stokes parameters into the input frequency
      !! axis\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!     p_red(Redc_class): Structure with redistribution input
      !!                        frequency data for a given input
      !!                        transition\n
      !!       Fed(Reda_class): Structure with redistribution
      !!                        output frequency data\n
      !!       Red(Redb_class): Structure with redistribution input
      !!                        frequency data\n
      !!         Mif0(integer): First frequency for this CPU\n
      !!         Mif1(integer): Last frequency for this CPU\n
      !!      omega(double(:)): Frequency array\n
      !!            vx(double): Velocity vector along X\n
      !!            vy(double): Velocity vector along Y\n
      !!            vz(double): Velocity vector along Z\n
      !!         lvel(logical): If dynamic node\n
      !!      Stkin(double(:)): Interpolated intensity\n
      !!    Stk(double(:,:,:)): Original intensity
      subroutine getStkinI(Geom,p_red,Fed,Red,Mif0,Mif1,omega, &
                           vx,vy,vz,lvel,Stkin,Stk)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(Redc_class), intent(in), pointer:: p_red
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      logical, intent(in):: lvel
      integer, intent(in):: Mif0,Mif1
      double precision, intent(in):: vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:,:), intent(in):: Stk
      double precision, dimension(:), allocatable, intent(out):: Stkin

      ! Local

      integer:: ith1,iph1,nmfreq,iran,ifreq,iifreq,lifreq,ibfreq
      integer:: jfreq,jjfreq0,jjfreq,kkfreq0,kkfreq,nblock

      double precision:: y0,dy,dx,vfac1,cost,sint,cosc,sinc

      ! Pointer

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! Get J size
      nmfreq = 0

      ! Run over all output frequencies
      iifreq = 0
      do iran=1,Fed%nran
        do ifreq=Fed%if0(iran),Fed%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! MPI
          if (iifreq.gt.Mif1) exit
          if (iifreq.lt.Mif0) cycle

          ! Input frequency number
          p_mfreq => p_red%mfreq(iifreq)

          ! Add size
          if (p_mfreq.gt.0) nmfreq = nmfreq + p_mfreq

        end do
      end do ! Output frequencies

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate intensity
      if (axiali) then
        allocate(Stkin(nmfreq*Geom%nTh))
      else
        allocate(Stkin(nmfreq*(Geom%nPh2*Geom%nTh)))
      end if

      ! If axial
      if (axiali) then

        ! If dynamic
        if (lvel) then

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! For each input direction
              do ith1=1,Geom%nth

                ! Get director cosines
                cost = Geom%V_mu(ith1)

                ! Calculate Doppler shift factor
                vfac1 = 1d0 - vz*cost

                ! We will be using the inverse
                vfac1 = 1d0/vfac1

                ! Reset the search frequency
                lifreq = Red%ggf0

                ! For each input frequency
                do jfreq=1,p_mfreq

                  ! Advance indexes
                  jjfreq = jjfreq0 + jfreq
                  kkfreq = kkfreq0 + jfreq

                  ! If out of range, take the value at the
                  ! boundary
                  if (p_red%omega(jjfreq)*vfac1.le. &
                      omega(Red%ggf0)+TINYO) then

                    ! We are still looking in the first one
                    lifreq = Red%ggf0

                    ! The index to take is 1
                    Stkin(kkfreq) = Stk(Red%ggf0,1,ith1)

                  ! If out of range, take the value at the
                  ! boundary
                  else if (p_red%omega(jjfreq)*vfac1.ge. &
                           (omega(Red%ggf1) - TINYO)) then

                    ! We are in the last frequency
                    lifreq = Red%ggf1

                    ! The index to take is nfreq
                    Stkin(kkfreq) = Stk(Red%ggf1,1,ith1)

                  ! If within the boundaries
                  else

                    ! Search between the last found frequency and
                    ! all but the boundary
                    do ibfreq=lifreq,nfreq-1

                      ! If small, skip
                      if (p_red%omega(jjfreq)*vfac1.lt. &
                          omega(ibfreq)-TINYO) cycle

                      ! Between
                      if (p_red%omega(jjfreq)*vfac1.lt. &
                          omega(ibfreq+1)) then

                        ! We are in the found frequency
                        lifreq = ibfreq

                        ! Same
                        if (abs(p_red%omega(jjfreq)*vfac1 - &
                            omega(ibfreq)).lt.TINYO) then

                          ! This frequency gives us the value
                          Stkin(kkfreq) = Stk(lifreq,1,ith1)

                        ! Actually between
                        else

                          ! The first index is the lower
                          y0 = Stk(lifreq,1,ith1)

                          ! The second index is the upper
                          dy = Stk(lifreq+1,1,ith1) - y0

                          ! Inverse of the distance
                          ! between the two outputs
                          dx = (p_red%omega(jjfreq)*vfac1 - &
                                omega(lifreq))/ &
                               (omega(lifreq+1) - omega(lifreq))

                          ! Interpolate
                          Stkin(kkfreq) = dx*dy + y0

                        end if ! Same or between

                        ! Found, leave loop
                        exit

                      end if ! Found range

                    end do ! Run output frequencies

                  end if ! Check if out of limits

                end do ! Run input frequencies

                ! Update indexes
                kkfreq0 = kkfreq0 + p_mfreq

              end do ! Input polar

              ! Update indexes
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! Output frequencies
          end do ! Output frequency ranges

        ! If static
        else

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! Reset the search frequency
              lifreq = Red%ggf0

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq

                ! If out of range, take the value at the
                ! boundary
                if (p_red%omega(jjfreq).le. &
                    omega(Red%ggf0)+TINYO) then

                  ! We are still looking in the first one
                  lifreq = Red%ggf0

                  ! Inclinations
                  do ith1=1,Geom%nTh

                    ! Output index
                    kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                    ! The index to take is 1
                    Stkin(kkfreq) = Stk(Red%ggf0,1,ith1)

                  end do

                ! If out of range, take the value at the
                ! boundary
                else if (p_red%omega(jjfreq).ge. &
                         (omega(Red%ggf1) - TINYO)) then

                  ! We are in the last frequency
                  lifreq = Red%ggf1

                  ! Inclinations
                  do ith1=1,Geom%nTh

                    ! Output index
                    kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                    ! The index to take is nfreq
                    Stkin(kkfreq) = Stk(Red%ggf1,1,ith1)

                  end do

                ! If within the boundaries
                else

                  ! Search between the last found frequency and
                  ! all but the boundary
                  do ibfreq=lifreq,nfreq-1

                    ! If small, skip
                    if (p_red%omega(jjfreq).lt. &
                        omega(ibfreq)-TINYO) cycle

                    ! Between
                    if (p_red%omega(jjfreq).lt.omega(ibfreq+1)) then

                      ! We are in the found frequency
                      lifreq = ibfreq

                      ! Same
                      if (abs(p_red%omega(jjfreq)- &
                              omega(ibfreq)).lt.TINYO) then

                        ! Inclinations
                        do ith1=1,Geom%nTh

                          ! Output index
                          kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                          ! The index to take is nfreq
                          Stkin(kkfreq) = Stk(lifreq,1,ith1)

                        end do

                      ! Actually between
                      else

                        ! Inverse of the distance
                        ! between the two outputs
                        dx = (p_red%omega(jjfreq) - &
                              omega(lifreq))/ &
                             (omega(lifreq+1) - omega(lifreq))

                        ! Inclinations
                        do ith1=1,Geom%nTh

                          ! Output index
                          kkfreq = kkfreq0 + p_mfreq*(ith1-1) + jfreq

                          ! The first index is the lower
                          y0 = Stk(lifreq,1,ith1)

                          ! The second index is the upper
                          dy = Stk(lifreq+1,1,ith1) - y0

                          ! Interpolate
                          Stkin(kkfreq) = dx*dy + y0

                        end do

                      end if ! Same or between

                      ! Found, leave loop
                      exit

                    end if ! Found range

                  end do ! Run output frequencies

                end if ! Check if out of limits

              end do ! Run input frequencies

              ! Update frequency index
              jjfreq0 = jjfreq0 + p_mfreq

              ! Update Stokes index
              kkfreq0 = kkfreq0 + p_mfreq*Geom%nTh

            end do ! Output frequencies
          end do ! Output frequency ranges

        end if ! Not dynamic

      ! Non-axially symmetric
      else

        ! If dynamic
        if (lvel) then

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! For each input direction
              do ith1=1,Geom%nth
                do iph1=1,Geom%nph

                  ! Get director cosines
                  cost = Geom%V_mu(ith1)
                  sint = sqrt(1d0 - cost*cost)
                  cosc = Geom%v_mux(iph1)
                  sinc = Geom%v_muy(iph1)*sqrt(1d0 - cosc*cosc)

                  ! Calculate Doppler shift factor
                  vfac1 = 1d0 - vx*sint*cosc - &
                                vy*sint*sinc - &
                                vz*cost

                  ! We will be using the inverse
                  vfac1 = 1d0/vfac1

                  ! Reset the search frequency
                  lifreq = Red%ggf0

                  ! For each input frequency
                  do jfreq=1,p_mfreq

                    ! Advance indexes
                    jjfreq = jjfreq0 + jfreq
                    kkfreq = kkfreq0 + jfreq

                    ! If out of range, take the value at the
                    ! boundary
                    if (p_red%omega(jjfreq)*vfac1.le. &
                        omega(Red%ggf0)+TINYO) then

                      ! We are still looking in the first one
                      lifreq = Red%ggf0

                      ! The index to take is 1
                      Stkin(kkfreq) = Stk(Red%ggf0,iph1,ith1)

                    ! If out of range, take the value at the
                    ! boundary
                    else if (p_red%omega(jjfreq)*vfac1.ge. &
                             (omega(Red%ggf1) - TINYO)) then

                      ! We are in the last frequency
                      lifreq = Red%ggf1

                      ! The index to take is nfreq
                      Stkin(kkfreq) = Stk(Red%ggf1,iph1,ith1)

                    ! If within the boundaries
                    else

                      ! Search between the last found frequency and
                      ! all but the boundary
                      do ibfreq=lifreq,nfreq-1

                        ! If small, skip
                        if (p_red%omega(jjfreq)*vfac1.lt. &
                            omega(ibfreq)-TINYO) cycle

                        ! Between
                        if (p_red%omega(jjfreq)*vfac1.lt. &
                            omega(ibfreq+1)) then

                          ! We are in the found frequency
                          lifreq = ibfreq

                          ! Same
                          if (abs(p_red%omega(jjfreq)*vfac1 - &
                                  omega(ibfreq)).lt.TINYO) then

                            ! This frequency gives us the value
                            Stkin(kkfreq) = Stk(lifreq,iph1,ith1)

                          ! Actually between
                          else

                            ! The first index is the lower
                            y0 = Stk(lifreq,iph1,ith1)

                            ! The second index is the upper
                            dy = Stk(lifreq+1,iph1,ith1) - y0

                            ! Inverse of the distance
                            ! between the two outputs
                            dx = (p_red%omega(jjfreq)*vfac1 - &
                                  omega(lifreq))/ &
                                 (omega(lifreq+1) - omega(lifreq))

                            ! Interpolate
                            Stkin(kkfreq) = dx*dy + y0

                          end if ! Same or between

                          ! Found, leave loop
                          exit

                        end if ! Found range

                      end do ! Run output frequencies

                    end if ! Check if out of limits

                  end do ! Run input frequencies

                  ! Update index
                  kkfreq0 = kkfreq0 + p_mfreq

                end do ! Input azimuth
              end do ! Input polar

              ! Update index
              jjfreq0 = jjfreq0 + p_mfreq

            end do ! Output frequencies
          end do ! Output frequency ranges

        ! If static
        else

          ! Initialize index
          jjfreq0 = 0
          kkfreq0 = 0

          ! For each output frequency
          iifreq = 0
          do iran=1,Fed%nran
            do ifreq=Fed%if0(iran),Fed%if1(iran)

              ! Advance index
              iifreq = iifreq + 1

              ! MPI
              if (iifreq.gt.Mif1) exit
              if (iifreq.lt.Mif0) cycle

              ! Input frequency number
              p_mfreq => p_red%mfreq(iifreq)

              ! Skip empty
              if (p_mfreq.lt.1) cycle

              ! Reset the search frequency
              lifreq = Red%ggf0

              ! For each input frequency
              do jfreq=1,p_mfreq

                ! Advance indexes
                jjfreq = jjfreq0 + jfreq

                ! If out of range, take the value at the
                ! boundary
                if (p_red%omega(jjfreq).le. &
                    omega(Red%ggf0)+TINYO) then

                  ! We are still looking in the first one
                  lifreq = Red%ggf0

                  ! Block counter
                  nblock = -1

                  ! Directions
                  do ith1=1,Geom%nTh
                    do iph1=1,Geom%nPh

                      ! Add block
                      nblock = nblock + 1

                      ! Output index
                      kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                      ! The index to take is 1
                      Stkin(kkfreq) = Stk(Red%ggf0,iph1,ith1)

                    end do
                  end do

                ! If out of range, take the value at the
                ! boundary
                else if (p_red%omega(jjfreq).ge. &
                         (omega(Red%ggf1) - TINYO)) then

                  ! We are in the last frequency
                  lifreq = Red%ggf1

                  ! Block counter
                  nblock = -1

                  ! Directions
                  do ith1=1,Geom%nTh
                    do iph1=1,Geom%nPh

                      ! Add block
                      nblock = nblock + 1

                      ! Output index
                      kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                      ! The index to take is nfreq
                      Stkin(kkfreq) = Stk(Red%ggf1,iph1,ith1)

                    end do
                  end do

                ! If within the boundaries
                else

                  ! Search between the last found frequency and
                  ! all but the boundary
                  do ibfreq=lifreq,nfreq-1

                    ! If small, skip
                    if (p_red%omega(jjfreq).lt. &
                        omega(ibfreq)-TINYO) cycle

                    ! Between
                    if (p_red%omega(jjfreq).lt.omega(ibfreq+1)) then

                      ! We are in the found frequency
                      lifreq = ibfreq

                      ! Same
                      if (abs(p_red%omega(jjfreq)- &
                              omega(ibfreq)).lt.TINYO) then

                        ! Block counter
                        nblock = -1

                        ! Directions
                        do ith1=1,Geom%nTh
                          do iph1=1,Geom%nPh

                            ! Add block
                            nblock = nblock + 1

                            ! Output index
                            kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                            ! This frequency gives us the value
                            Stkin(kkfreq) = Stk(lifreq,iph1,ith1)

                          end do
                        end do

                      ! Actually between
                      else

                        ! Inverse of the distance
                        ! between the two outputs
                        dx = (p_red%omega(jjfreq) - &
                              omega(lifreq))/ &
                             (omega(lifreq+1) - omega(lifreq))

                        ! Block counter
                        nblock = -1

                        ! Directions
                        do ith1=1,Geom%nTh
                          do iph1=1,Geom%nPh

                            ! Add block
                            nblock = nblock + 1

                            ! The first index is the lower
                            y0 = Stk(lifreq,iph1,ith1)

                            ! The second index is the upper
                            dy = Stk(lifreq+1,iph1,ith1) - y0

                            ! Output index
                            kkfreq = kkfreq0 + jfreq + p_mfreq*nblock

                            ! Interpolate
                            Stkin(kkfreq) = dx*dy + y0

                          end do
                        end do

                      end if ! Same or between

                      ! Found, leave loop
                      exit

                    end if ! Found range

                  end do ! Run output frequencies

                end if ! Check if out of limits

              end do ! Run input frequencies

              ! Update frequency index
              jjfreq0 = jjfreq0 + p_mfreq

              ! Update Stokes index
              kkfreq0 = kkfreq0 + p_mfreq*Geom%nTh*Geom%nPh

            end do ! Output frequencies
          end do ! Output frequency ranges

        end if ! Static
      end if ! Non-axial

      ! Nullify
      nullify(p_mfreq)

      end subroutine getStkinI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Interpolate the frequency dependent mean intensity into the
      !! input frequency axis\n
      !!  p_red(Redc_class): Structure with redistribution input
      !!                     frequency data for a given input
      !!                     transition\n
      !!    Fed(Reda_class): Structure with redistribution output
      !!                     frequency data\n
      !!    Red(Redb_class): Structure with redistribution input
      !!                     frequency data\n
      !!      Mif0(integer): First frequency for this CPU\n
      !!      Mif1(integer): Last frequency for this CPU\n
      !!    nmfreq(integer): Size of frequency space\n
      !!   omega(double(:)): Frequency array\n
      !!     Jin(double(:)): Interpolated frequency dependent mean
      !!                     intensity\n
      !!     J00(double(:)): Original frequency dependent mean
      !!                     intensity
      subroutine getJin(p_red,Fed,Red,Mif0,Mif1,nmfreq,omega,Jin,J00)

      ! I/O

      type(Redc_class), intent(in), pointer:: p_red
      type(Reda_class), intent(in):: Fed
      type(Redb_class), intent(in):: Red
      integer, intent(in):: nmfreq,Mif0,Mif1
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), allocatable, intent(out):: Jin
      double precision, dimension(Red%ggf0:Red%ggf1), intent(in):: J00

      ! Local

      integer:: lifreq,ibfreq,iran,ifreq,iifreq,jfreq,jjfreq0,jjfreq

      double precision:: y0,dy,dx

      ! Pointer

      integer, pointer:: p_mfreq


      ! Nullify
      nullify(p_mfreq)

      ! If no size, return
      if (nmfreq.lt.1) return

      ! Allocate J
      allocate(Jin(nmfreq))

      ! Initialize index
      jjfreq0 = 0

      ! For each output frequency
      iifreq = 0
      do iran=1,Fed%nran
        do ifreq=Fed%if0(iran),Fed%if1(iran)

          ! Advance index
          iifreq = iifreq + 1

          ! MPI
          if (iifreq.gt.Mif1) exit
          if (iifreq.lt.Mif0) cycle

          ! Input frequency number
          p_mfreq => p_red%mfreq(iifreq)

          ! Reset the search frequency
          lifreq = Red%ggf0

          ! For each input frequency
          do jfreq=1,p_mfreq

            ! Advance indexes
            jjfreq = jjfreq0 + jfreq

            ! If out of range, take the value at the
            ! boundary
            if (p_red%omega(jjfreq).le.omega(Red%ggf0)+TINYO) then

              ! We are still looking in the first one
              lifreq = Red%ggf0

              ! The index to take is 1
              Jin(jjfreq) = J00(Red%ggf0)

            ! If out of range, take the value at the
            ! boundary
            else if (p_red%omega(jjfreq).ge. &
                     (omega(Red%ggf1) - TINYO)) then

              ! We are in the last frequency
              lifreq = Red%ggf1

              ! The index to take is nfreq
              Jin(jjfreq) = J00(Red%ggf1)

            ! If within the boundaries
            else

              ! Search between the last found frequency and
              ! all but the boundary
              do ibfreq=lifreq,Red%ggf1-1

                ! Skip before
                if (p_red%omega(jjfreq).lt.omega(ibfreq)-TINYO) cycle

                ! Between
                if (p_red%omega(jjfreq).lt.omega(ibfreq+1)) then

                  ! We are in the found frequency
                  lifreq = ibfreq

                  ! Same
                  if (abs(p_red%omega(jjfreq) - &
                          omega(ibfreq)).lt.TINYO) then

                    ! This frequency gives us the value
                    Jin(jjfreq) = J00(lifreq)

                  ! Actually between
                  else

                    ! Inverse of the distance
                    ! between the two outputs
                    dx = (p_red%omega(jjfreq) - omega(lifreq))/ &
                         (omega(lifreq+1) - omega(lifreq))

                    ! The first index is the lower
                    y0 = J00(ibfreq)

                    ! Difference with next
                    dy = J00(ibfreq+1) - y0

                    ! Interpolate
                    Jin(jjfreq) = dy*dx + y0

                  end if ! Same or between

                  ! Found, leave loop
                  exit

                end if ! Found range

              end do ! Run output frequencies

            end if ! Check if out of limits

          end do ! Run input frequencies

          ! Update index in general
          jjfreq0 = jjfreq0 + p_mfreq

        end do ! Output frequencies
      end do ! Output frequency ranges

      ! Free pointers
      nullify(p_mfreq)

      end subroutine getJin

!#####################################################################
!#####################################################################
!#####################################################################

      end module rtcoeffiaux_mod
