      !> Get the radiation field tensors for the SEE
      module jcalccle_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     01/10/2022
!  Last version:
!     13/12/2024 V4.0.0
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     13/12/2024:    V4.0.0 - Updated for the new dimensions in the
!                             geometrical tensors (TdPA)
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
!    o Implement non-axial input spectra
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  getSEEJ
!    Compute the frequency integrated radiation field tensors for the
!  SEE
!
!  get_bottom_spect
!    Get boundary condition for the radiation field from the input
!  radiation
!
!  get_bottom_allen
!    Get boundary condition for the radiation field from Allen's
!  tabulation
!
!  getI
!    Get radiation intensity from Allen's tabulation
!
!  get_CLV
!    Get CLV coefficients for a given frequency from Allen's
!  tabulation
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use btens_mod
      use commons_mod
      use fieldb_mod
      use inter_mod
      use math_mod
      use parameters_mod , only: cZero , convF, IPI , PI , TINYB , &
                                 TINYO , cSaha , fktoJ , kb
      use planck_mod
      use profile_mod
      use stens_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the frequency integrated radiation field tensors for
      !! the SEE\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!              Trad(double): Radiation temperature of the
      !!                            illuminating disk\n
      !!        use_allen(logical): If using Allen's tabulation for
      !!                            input radiation\n
      !!      flat_cle_in(logical): If assuming flat spectrum input\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!        Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                            and J-symbols\n
      !!     Frec(Frequency_class): Structure with frequency data\n
      !!        spect(spect_class): Structure with the input spectra
      !!                            data\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!  GeomP(Coronapoint_class): Structure with geometric data for
      !!                            a CLE node\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!     JKQC(dcomplex(:,:,:)): Radiation field tensors with
      !!                            frequency dependence\n
      !!      JKQ(dcomplex(:,:,:)): Radiation field tensors
      !!                            integrated over the absorption
      !!                            profile\n
      !!        JPhot(double(:,:)): Intensity integrals in the
      !!                            photoionization rates
      subroutine getSEEJ(Atom,Atmo,Trad,use_allen,flat_cle_in, &
                         Bfield,Flgsg,Frec,spect,Geom,GeomP,MPID, &
                         JKQC,JKQ,Jphot)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Spect_class), intent(inout):: spect
      type(Fctsg_class), intent(inout):: Flgsg
      type(Frequency_class), intent(in):: Frec
      type(Coronapoint_class), intent(in):: GeomP
      type(Geometry_class), intent(inout):: Geom
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: use_allen,flat_cle_in
      double precision, intent(in):: Trad
      double precision, dimension(nxphot,2), intent(out):: Jphot
      complex(kind=8), dimension(-2:2,0:2,nxtran), intent(out):: JKQ
      complex(kind=8), dimension(-2:2,0:2,nfreq), intent(out):: JKQC

      ! Local

      logical:: first

      integer:: ia,ith,iph,iph1,ifreq,jfreq,cljfreq,ljfreq,iproc,istk
      integer:: if0,if1,jf0,jf1,lf0,lf1,K,iQ
      integer:: itermu,iterml,iJu,iJl,jdir
      integer:: jtran,ktran,fjtran,ffjtran,ffktran
      integer:: iaux,ntasks,itask,nbag,ios
      integer, dimension(0:nproc-1):: nf,pif0,pif1

      double precision:: CYp1h,CYm1h,theta,phi,vfac
      double precision:: ne,nh,sqT,Bs,Bt,Bp,ct,st,cc,sc
      double precision:: WA,W0,W1,DwT,Dw,gl
      double precision:: omegain,dx,feta,Tfeta
      double precision:: al,au,aul,at,Dfreq,u1,u2
      double precision:: prof,c0,c1,c3,Saha,arg,exu,I0
      double precision, dimension(:), allocatable:: dy
      double precision, dimension(:), allocatable:: b,c,d
      double precision, dimension(:,:,:), allocatable:: dyv
      double precision, dimension(:,:,:), allocatable:: buffer
      double precision, dimension(:,:,:,:), allocatable:: inspect

      complex(kind=8), dimension(:,:,:), allocatable:: integr

      ! Pointers

      double precision, pointer:: T,vx,vy,vz,vmi


      ! Initialize
      Jphot = 0d0
      JKQ = cZero
      JKQC = cZero

      ! Get frequency limits for this CPU
      if0 = MPID%iif0(pid)
      if1 = MPID%iif1(pid)

      ! Get thermal quantities
      T => Atmo%T(1)
      vx => Atmo%vx(1)
      vy => Atmo%vy(1)
      vz => Atmo%vz(1)
      vmi => Atmo%vmi(1)
      ne = Atmo%ne(1)
      nh = Atmo%nHa(1)
      sqT = sqrt(T)
      Bs = Bfield%Bstrength(1)
      Bt = Bfield%Btheta(1)
      Bp = Bfield%Bphi(1)

      ! Saha non-line-dependent part
      Saha = cSaha*ne/(T**(1.5d0))
      arg = fktoJ/kb/T
      c0 = 4d-8*pi/convF


      !
      ! If there is an input spectrum or we are doing PRD, set
      ! the quadrature for this height
      !
      if (spect%valid.or.PRD) then

        ! Half the CLV factor 1 (1 - CY)
        CYm1h = GeomP%CLV(1)*0.5d0
        CYp1h = (GeomP%CY + 1d0)*0.5d0

        ! Weights are given by half the cosine of the gamma angle
        Geom%W_mu = CYm1h*Geom%W_gauss

        ! Get the polar nodes for this particular height
        Geom%V_mu = CYm1h*Geom%V_gauss + CYp1h

        ! Get the polar nodes on the disk
        Geom%V_mu_disk = sqrt(1d0 - (1d0 - Geom%V_mu*Geom%V_mu)/ &
                                    (GeomP%SY*GeomP%SY))

        ! Allocate geometrical tensors
        allocate(Geom%TS(0:3,-2:2,0:2,Geom%nPh*Geom%nTh))
        allocate(Geom%TB(0:3,-2:2,0:2,Geom%nPh*Geom%nTh,1))

        ! Inititlize direction index
        jdir = 0

        ! For each polar direction
        do ith=1,Geom%nTh

          ! Get angle
          theta = acos(Geom%V_mu(ith))

          ! For each azimuth
          do iph=1,Geom%nPh

            ! Get angle
            phi = Geom%V_phi(iph)

            ! Advance direction index
            jdir = jdir + 1

            ! Get geometrical tensors in the vertical reference frame
            call Stens(theta,phi,GeomP%geom(3),Flgsg, &
                       Geom%TS(:,:,:,jdir))

            ! If there is magnetic field
            if (Bs.gt.TINYB) then

              ! Rotate to the magnetic field reference frame
              call Btens(Geom%TS(:,:,:,jdir), &
                         Geom%TB(:,:,:,jdir,1),Flgsg,Bt,Bp)

            ! No magnetic field
            else

              ! Copy vertical reference frame tensors
              Geom%TB(:,:,:,jdir,1) = Geom%TS(:,:,:,jdir)

            end if ! Magnetic field

          end do ! Phi angles
        end do ! Polar angles

      end if ! Input spectra or PRD

      !!!!!!!!!!!!!!!!!!!!!
      !                   !
      ! If valid spectrum !
      !                   !
      !!!!!!!!!!!!!!!!!!!!!
      if (spect%valid) then

        ! If something to allocate
        if (Frec%Mlif1(pid).ge.Frec%Mlif0(pid)) then

          ! If axial
          if (axial) then

            ! Allocate
            allocate(integr(0:0,0:2,Frec%Mlif0(pid):Frec%Mlif1(pid)))

          ! If non-axial
          else

            ! Allocate
            allocate(integr(0:2,0:2,Frec%Mlif0(pid):Frec%Mlif1(pid)))

          end if ! Axial

          ! Initialize to zero
          integr = cZero

        end if ! Something to allocate

        !
        ! Allocate intensity in grid
        !

        ! If dynamic or non-axial input
        if (dyn.or..not.spect%axial) then

          ! Allocate
          allocate(inspect(0:spect%nstk,if0:if1,Geom%nPh,Geom%nTh))

        ! Static and axial input
        else

          ! Allocate
          allocate(inspect(0:spect%nstk,if0:if1,1,Geom%nTh))

        end if ! Dynamic or non-axial input

        ! If MPI, decide sizes of jobs
        if (nproc.gt.1) then

          ! Axial
          if (spect%axial) then

            ! Number of tasks to perform
            ntasks = (spect%nstk+1)*spect%nfreq

            ! Angular size of package
            nbag = Geom%nTh

            ! Allocate
            allocate(buffer(1,Geom%nTh,ntasks))

          ! Not axial
          else

            ! Number of tasks to perform
            ntasks = (spect%nstk+1)*spect%nfreq

            ! Angular size of package
            nbag = Geom%nTh*Geom%nPh

            ! Allocate
            allocate(buffer(Geom%nPh,Geom%nTh,ntasks))

          end if ! Axial
        end if ! MPI

        !
        ! Allocate interpolation auxiliar
        !

        ! If velocities
        if (dyn) then

          ! Allocate
          allocate(dy(0:spect%nstk))

        ! If static
        else

          ! If axial input
          if (spect%axial) then

            ! Allocate
            allocate(dyv(0:spect%nstk,1,Geom%nTh))

          ! Not axial
          else

            ! Allocate
            allocate(dyv(0:spect%nstk,Geom%nPh,Geom%nTh))

          end if ! Axial
        end if ! dynamic

        ! Allocate spline coefficients
        allocate(b(spect%nmu),c(spect%nmu),d(spect%nmu))

        ! If MPI
        if (nproc.gt.1) then

          ! Let's split the work here
          iaux = ntasks/nproc

          ! We need at least one frequency per processor
          if(iaux.lt.1)then
            umsg = 'Too many processors for this frequency grid'
            call aborted
            return
          end if

          ! Determine number of tasks per CPU
          nf = iaux

          ! Complex splitting if they do not coincide in size
          if(iaux*nproc.ne.ntasks)then

            ! Add one to the first CPUs until filled
            iaux = ntasks - iaux*nproc
            do iproc=0,iaux-1
              nf(iproc) = nf(iproc) + 1
            end do

          end if ! No perfect split

          ! First CPU indexes
          pif0(0) = 1
          pif1(0) = nf(0)

          ! For each CPU
          do iproc=1,nproc-1

            ! Get indexes
            pif0(iproc) = pif1(iproc-1) + 1
            pif1(iproc) = pif0(iproc) + nf(iproc) - 1

          end do ! CPUs

          ! If input is axial
          if (spect%axial) then

            ! Initialize task index
            itask = 0

            ! For each polarization state
            do istk=0,spect%nstk

              ! For each frequency in the input (split)
              do ifreq=1,spect%nfreq

                ! Advance index
                itask = itask + 1

                ! If below limit, continue
                if (itask.lt.pif0(pid)) cycle
                ! If above limit, done
                if (itask.gt.pif1(pid)) exit

                ! Get splines
                call spline(spect%mu,spect%stokes(:,1,ifreq,istk), &
                            b,c,d,spect%nmu)

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! Get interpolation
                  buffer(1,ith,itask) = &
                              ispline(Geom%V_mu_disk(ith), &
                                      spect%mu, &
                                      spect%stokes(:,1,ifreq,istk), &
                                      b,c,d,spect%nmu)

                end do ! Polar directions
              end do ! Frequencies
            end do ! Stokes parameters

            ! Scale for message sending
            nf = nf*nbag
            do iproc=0,nproc-1
              pif1(iproc) = iproc*nf(iproc)
            end do

            ! Share info
            call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, &
                                buffer(1,1,1),nf,pif1, &
                                MPI_DOUBLE_PRECISION, &
                                MPI_COMM_RT,ierr)

            !
            ! Reorder
            !

            ! Initialize task index
            itask = 0

            ! For each polarization state
            do istk=0,spect%nstk

              ! For each frequency in the input (split)
              do ifreq=1,spect%nfreq

                ! Advance index
                itask = itask + 1

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! Get interpolation from buffer
                  spect%mustokes(istk,ifreq,1,ith) = &
                                                   buffer(1,ith,itask)

                end do ! Polar directions
              end do ! Frequencies
            end do ! Stokes parameters

          ! Input not axial
          else

            ! Not implemented!
            umsg = 'Non-axial input spectrum is not implemented'
            call aborted
            return

          end if ! Axial

        ! Serial
        else

          ! If input is axial
          if (spect%axial) then

            ! For each polarization state
            do istk=0,spect%nstk

              ! For each frequency in the input (split)
              do ifreq=1,spect%nfreq

                ! Get splines
                call spline(spect%mu,spect%stokes(:,1,ifreq,istk), &
                            b,c,d,spect%nmu)

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! Get interpolation
                  spect%mustokes(istk,ifreq,1,ith) = &
                              ispline(Geom%V_mu_disk(ith), &
                                      spect%mu, &
                                      spect%stokes(:,1,ifreq,istk), &
                                      b,c,d,spect%nmu)

                end do ! Polar directions
              end do ! Frequencies
            end do ! Stokes parameters

          ! Input not axial
          else

            ! Not implemented!
            umsg = 'Non-axial input spectrum is not implemented'
            call aborted
            return

          end if ! Axial
        end if ! MPI

        !
        ! Now integrate in relevant frequencies
        !

        ! If dynamic
        if (dyn) then

          ! For each polar direction
          do ith=1,Geom%nTh

            ! Cosine and sine
            ct = Geom%V_mu(ith)
            st = sqrt(1d0 - ct*ct)

            ! If axial, define vfac
            if (axial) vfac = 1d0 + vz*ct

            ! For each azimuthal direction
            do iph=1,Geom%nPh

              ! Index in input
              if (spect%axial) then
                iph1 = 1
              else
                iph1 = iph
              end if

              ! If not axial
              if (.not.axial) then

                ! Get director cosines
                cc = Geom%V_mux(iph)
                sc = Geom%V_muy(iph)*sqrt(1d0 - cc*cc)

                ! Get Doppler factor
                vfac = 1d0 + vz*ct + &
                             vx*st*cc + &
                             vy*st*sc

              end if

              ! Last searched frequency
              cljfreq = 1

              ! Indicate it is first
              first = .True.

              ! For each frequency
              do ifreq=if0,if1

                ! Update last jfreq
                ljfreq = cljfreq

                ! Get displaced source frequency
                omegain = Frec%omega(ifreq)*vfac

                ! If below lower limit
                if (omegain.le.spect%omega(1)) then

                  ! Constant extension
                  inspect(:,ifreq,iph,ith) = &
                                          spect%mustokes(:,1,iph1,ith)
                  cycle

                ! Or above limits
                else if (omegain.ge.spect%omega(spect%nfreq)) then

                  ! Constant extension
                  inspect(:,ifreq,iph,ith) = &
                                spect%mustokes(:,spect%nfreq,iph1,ith)
                  cycle

                end if

                ! Look for the surrounding frequencies
                do jfreq=ljfreq,spect%nfreq-1

                  ! If below current limit, skip
                  if (omegain.gt.spect%omega(jfreq).and. &
                      omegain.lt.spect%omega(jfreq+1)) then

                    ! If different than last time, calculate
                    ! linear interpolation coefficients
                    if (cljfreq.ne.ljfreq.or.first) then

                      ! Get interpolation data
                      dx = 1d0/ &
                         (spect%omega(jfreq+1) - spect%omega(jfreq))
                      dy = spect%mustokes(:,jfreq+1,iph1,ith) - &
                           spect%mustokes(:,jfreq,iph1,ith)

                      ! Flag it is not the first
                      if (first) first = .False.

                    end if ! Different than last time or first

                    ! Linear interpolation
                    inspect(:,ifreq,iph,ith) = &
                                  spect%mustokes(:,jfreq,iph1,ith) + &
                                  dy*dx*(omegain - spect%omega(jfreq))
                    ! Finish
                    exit

                  ! If beyond next frequency, continue
                  else if (omegain.gt.spect%omega(jfreq+1)) then

                    ! Update current last searched frequency
                    cljfreq = jfreq+1

                    ! Continue searching
                    cycle

                  ! Exactly the same than current frequency
                  else if (omegain.le.spect%omega(jfreq)) then

                    ! No interpolation needed
                    inspect(:,ifreq,iph,ith) = &
                                      spect%mustokes(:,jfreq,iph1,ith)

                    ! Finish
                    exit

                  ! Exactly the same than next frequency
                  else

                    ! Update current last searched frequency
                    cljfreq = jfreq + 1

                    ! No interpolation needed
                    inspect(:,ifreq,iph,ith) = &
                                    spect%mustokes(:,jfreq+1,iph1,ith)

                    ! Finish
                    exit

                  end if ! Searching condition

                end do ! input spectrum frequencies
              end do ! CPU frequencies
            end do ! Phi directions
          end do ! Polar directions

        ! Static
        else

          ! Last searched frequency
          cljfreq = 1

          ! Indicate it is first
          first = .True.

          ! For each frequency
          do ifreq=if0,if1

            ! Update last jfreq
            ljfreq = cljfreq

            ! Get displaced source frequency
            omegain = Frec%omega(ifreq)

            ! If below lower limit
            if (omegain.le.spect%omega(ljfreq)) then

              ! Constant extension
              inspect(:,ifreq,:,:) = spect%mustokes(:,1,:,:)
              cycle

            ! Or above limits
            else if (omegain.ge.spect%omega(spect%nfreq)) then

              ! Constant extension
              inspect(:,ifreq,:,:) = &
                                spect%mustokes(:,spect%nfreq,:,:)
              cycle

            end if

            ! Look for the surrounding frequencies
            do jfreq=ljfreq,spect%nfreq-1

              ! If below current limit, skip
              if (omegain.gt.spect%omega(jfreq).and. &
                  omegain.lt.spect%omega(jfreq+1)) then

                ! If different than last time, calculate
                ! linear interpolation coefficients
                if (cljfreq.ne.ljfreq.or.first) then

                  ! Interpolation data
                  dx = 1d0/ &
                     (spect%omega(jfreq+1) - spect%omega(jfreq))
                  dyv = spect%mustokes(:,jfreq+1,:,:) - &
                        spect%mustokes(:,jfreq,:,:)

                  ! Flag not the first anymore
                  if (first) first = .False.

                end if ! Different than last time or first

                ! Linear interpolation
                inspect(:,ifreq,:,:) = &
                              spect%mustokes(:,jfreq,:,:) + &
                              dyv*dx*(omegain - spect%omega(jfreq))
                ! Finish
                exit

              ! If beyond next frequency, continue
              else if (omegain.gt.spect%omega(jfreq+1)) then

                ! Update current last searched frequency
                cljfreq = jfreq+1

                ! Continue searching
                cycle

              ! Exactly the same than current frequency
              else if (omegain.le.spect%omega(jfreq)) then

                ! No interpolation needed
                inspect(:,ifreq,:,:) = spect%mustokes(:,jfreq,:,:)

                ! Finish
                exit

              ! Exactly the same than next frequency
              else

                ! Update current last searched frequency
                cljfreq = jfreq + 1

                ! No interpolation needed
                inspect(:,ifreq,:,:) = &
                                    spect%mustokes(:,jfreq+1,:,:)

                ! Finish
                exit

              end if ! Searching condition

            end do ! input spectrum frequencies
          end do ! CPU frequencies

        end if ! Dynamic

        !
        ! Time to integrate the input spectra in the radiation
        ! field tensors
        !

        ! If axial
        if (axial) then

          ! For each continuum frequency
          do ifreq=if0,if1

            ! Initialize direction index
            jdir = 0

            ! For each polar direction
            do ith=1,Geom%nTh

              ! Advance direction
              jdir = jdir + 1

              ! Weight
              WA = Geom%W_mu(ith)

              ! For each K
              do K=0,Krad

                ! Contribution
                JKQC(0,K,ifreq) = JKQC(0,K,ifreq) + &
                           WA*sum(inspect(0:spect%nstk,ifreq,1,ith)* &
                                  Geom%TS(0:spect%nstk,0,K,jdir))
              end do ! K

              ! Advance unused directions
              jdir = jdir - 1 + Geom%nPh

            end do ! Polar
          end do ! Frequency

          ! For each line frequency
          do ifreq=Frec%Mlif0(pid),Frec%Mlif1(pid)

            ! Initialize direction index
            jdir = 0

            ! For each polar direction
            do ith=1,Geom%nTh

              ! Advance direction
              jdir = jdir + 1

              ! Weight
              WA = Geom%W_mu(ith)

              ! For each K
              do K=0,Krad

                ! Get integr
                integr(0,K,ifreq) = integr(0,K,ifreq) + &
                           WA*sum(inspect(0:spect%nstk,ifreq,1,ith)* &
                                  Geom%TB(0:spect%nstk,0,K,jdir,1))
              end do ! K

              ! Advance unused directions
              jdir = jdir - 1 + Geom%nPh

            end do ! Polar
          end do ! Frequency

        ! Non-axial
        else

          ! For each continuum frequency
          do ifreq=if0,if1

            ! Initialize direction index
            jdir = 0

            ! For each polar direction
            do ith=1,Geom%nTh

              ! Azimuth
              do iph=1,Geom%nPh

                ! Advance direction
                jdir = jdir + 1

                ! Weight
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                ! For each K
                do K=0,Krad

                  ! For each Q >= 0
                  do iQ=0,K

                    ! Contribution
                    JKQC(iQ,K,ifreq) = JKQC(iQ,K,ifreq) + &
                         WA*sum(inspect(0:spect%nstk,ifreq,iph,ith)* &
                                Geom%TS(0:spect%nstk,iQ,K,jdir))
                  end do ! Q

                  ! For each Q<0
                  do iQ=1,K
                    JKQC(-iQ,K,ifreq) = Flgsg%sg(iQ)* &
                                        conjg(JKQC(iQ,K,ifreq))
                  end do

                end do ! K
              end do ! Azimuth
            end do ! Polar
          end do ! frequency

          ! For each line frequency
          do ifreq=Frec%Mlif0(pid),Frec%Mlif1(pid)

            ! Initialize direction index
            jdir = 0

            ! For each polar direction
            do ith=1,Geom%nTh

              ! Azimuth
              do iph=1,Geom%nPh

                ! Advance direction
                jdir = jdir + 1

                ! Weight
                WA = Geom%W_mu(ith)*Geom%W_mux(iph)

                ! For each K
                do K=0,Krad

                  ! For each Q >= 0
                  do iQ=0,K

                    ! Get integr
                    integr(iQ,K,ifreq) = integr(iQ,K,ifreq) + &
                         WA*sum(inspect(0:spect%nstk,ifreq,iph,ith)* &
                                Geom%TB(0:spect%nstk,iQ,K,jdir,1))
                  end do ! Q
                end do ! K
              end do ! Azimuth
            end do ! Polar
          end do ! frequency

        end if ! Axial

        ! For each atom
        do ia=1,nA

          ! Thermal part of the Doppler width
          DwT = Atom(ia)%cDopp*sqrt(T)

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! Get line indexes
            lf0 = Atom(ia)%ilf0(jtran)
            lf1 = Atom(ia)%ilf1(jtran)

            ! If absent line, skip
            if (lf0.gt.lf1) cycle

            ! Check if input radiation. If not, leave
            if (.not.Atom(ia)%bbspecin(jtran)) cycle

            ! True limits
            jf0 = max(if0,lf0)
            jf1 = min(if1,lf1)

            ! If not possible, skip
            if (jf0.gt.jf1) cycle

            ! Boundary weights
            W0 = Atom(ia)%W0(jtran)
            W1 = Atom(ia)%W1(jtran)

            ! Flattened line index
            ktran = jtran + Atom(ia)%tshift

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Damping
            al = Atom(ia)%damp(iterml,1)
            au = Atom(ia)%damp(itermu,1)
            aul = Atom(ia)%ldamp(jtran,1)

            ! Initialize sum of weights
            Tfeta = 0d0

            ! For each FS transition
            do fjtran=1,Atom(ia)%fst(jtran)%nt

              ! Idenfity involved levels
              iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
              iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

              ! Degeneration
              gl = 2d0*Atom(ia)%rJval(iJl,iterml) + 1d0

              ! Get the sequential index of this FS transition
              ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

              ! Get flattened line index
              ffktran = ffjtran + Atom(ia)%tfshift

              ! Get frequency of FS transition
              Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                      Atom(ia)%FSfreq(iJl,iterml)

              ! Add the microturbulence to Doppler width
              Dw = Dfreq*sqrt(DwT*DwT + vmi*vmi)

              ! Scale Dfreq to Dw
              Dfreq = Dfreq/Dw

              ! Absorptibity factor
              feta = Dfreq*Atom(ia)%fst(jtran)%Blu(iJl,iJu)*gl

              ! Add to total factor
              Tfeta = Tfeta + feta

              ! Add profile scale to feta
              feta = feta*1d-5*sqrt(IPI)/Dw

              ! Total damping
              at = (al + au + aul)/Dw

              ! Axial
              if (axial) then

                ! For each frequency
                do ifreq=jf0,jf1

                  ! If first point
                  if (ifreq.eq.lf0) then
                    WA = W0*feta
                  ! Last point
                  else if (ifreq.eq.lf1) then
                    WA = W1*feta
                  ! Any other
                  else
                    WA = Frec%W_freq(ifreq)*feta
                  end if

                  ! Get profile
                  call voigtI(Dfreq - Frec%omega(ifreq),at,prof)

                  ! Scale profile
                  WA = WA*prof

                  ! For each K
                  do K=0,Atom(ia)%Krad(jtran)

                    ! Contribution
                    JKQ(0,K,ktran) = JKQ(0,K,ktran) + &
                                      WA*integr(0,K,ifreq)
                  end do ! K
                end do ! Frequencies

              ! Not axial
              else

                ! For each frequency
                do ifreq=jf0,jf1

                  ! If first point
                  if (ifreq.eq.lf0) then
                    WA = W0*feta
                  ! Last point
                  else if (ifreq.eq.lf1) then
                    WA = W1*feta
                  ! Any other
                  else
                    WA = Frec%W_freq(ifreq)*feta
                  end if

                  ! Get profile
                  call voigtI(Dfreq - Frec%omega(ifreq),at,prof)

                  ! Scale profile
                  WA = WA*prof

                  ! For each K
                  do K=0,Atom(ia)%Krad(jtran)

                    ! For each Q>=0
                    do iQ=0,K

                      ! Contribution
                      JKQ(iQ,K,ktran) = JKQ(iQ,K,ktran) + &
                                        WA*integr(iQ,K,ifreq)
                    end do ! Q

                  end do ! K
                end do ! Frequencies

              end if ! Axial

            end do ! FS transitions

            ! Normalize
            Tfeta = 1d0/Tfeta
            JKQ(:,:,ktran) = JKQ(:,:,ktran)*Tfeta

            ! If not axial
            if (.not.axial) then

              ! For each K > 0
              do K=1,Atom(ia)%Krad(jtran)

                ! For each Q<0
                do iQ=1,K
                  JKQ(-iQ,K,ktran) = Flgsg%sg(iQ)* &
                                     conjg(JKQ(iQ,K,ktran))
                end do ! Q
              end do ! K

            end if ! Axial

          end do ! Transition

          ! If no photoionizations, skip
          if (Atom(ia)%nphot.le.0) cycle

          ! For each b-f transition
          do jtran=1,Atom(ia)%nphot

            ! Limits in frequency index of this transition
            lf0 = Atom(ia)%ipf0(jtran)
            lf1 = Atom(ia)%ipf1(jtran)

            ! If absent transition, leave
            if (lf0.gt.lf1) cycle

            ! Check if input radiation. If not, leave
            if (.not.Atom(ia)%bfspecin(jtran)) cycle

            ! True limits
            jf0 = max(if0,lf0)
            jf1 = min(if1,lf1)

            ! If not possible, skip
            if (jf0.gt.jf1) cycle

            ! Saha factor
            c3 = Saha*exp(Atom(ia)%phot(jtran)%edge*arg)* &
                 Atom(ia)%phot(jtran)%glu

            ! Boundary weights
            W0 = Atom(ia)%phot(jtran)%W0
            W1 = Atom(ia)%phot(jtran)%W1

            ! Flattened line index
            ktran = jtran + Atom(ia)%pshift

            ! Axial
            if (axial) then

              ! For each frequency
              do ifreq=jf0,jf1

                ! If first point
                if (ifreq.eq.lf0) then
                  c1 = W0
                ! Last point
                else if (ifreq.eq.lf1) then
                  c1 = W1
                ! Any other
                else
                  c1 = Frec%W_freq(ifreq)
                end if

                ! Weight with cross section and
                ! constants
                c1 = c0*c1*Atom(ia)%phot(jtran)%alpha(ifreq)/ &
                     Frec%omega(ifreq)

                ! Get exu
                exu = c3*diexp(c0*Frec%omega(ifreq))

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! Weight
                  WA = c1*Geom%W_mu(ith)*inspect(0,ifreq,iph,ith)

                  ! Add contribution
                  Jphot(ktran,1) = Jphot(ktran,1) + WA
                  Jphot(ktran,2) = Jphot(ktran,2) + WA*exu

                end do ! Polar
              end do ! Frequencies

            ! Not axial
            else

              ! For each frequency
              do ifreq=jf0,jf1

                ! If first point
                if (ifreq.eq.lf0) then

                  ! Get left weight
                  c1 = W0

                ! Last point
                else if (ifreq.eq.lf1) then

                  ! Get right weight
                  c1 = W1

                ! Any other
                else

                  ! Get weight
                  c1 = Frec%W_freq(ifreq)

                end if ! If extreme or not

                ! Weight with cross section and
                ! constants
                c1 = c0*c1*Atom(ia)%phot(jtran)%alpha(ifreq)/ &
                     Frec%omega(ifreq)

                ! Get exu
                exu = c3*diexp(c0*Frec%omega(ifreq))

                ! For each polar direction
                do ith=1,Geom%nTh

                  ! Azimuth
                  do iph=1,Geom%nPh

                    ! Weight
                    WA = c1*Geom%W_mu(ith)*Geom%W_mux(iph)* &
                         inspect(0,ifreq,iph,ith)

                    ! Add contribution
                    Jphot(ktran,1) = Jphot(ktran,1) + WA
                    Jphot(ktran,2) = Jphot(ktran,2) + WA*exu

                  end do ! Azimuth
                end do ! Polar
              end do ! Frequencies

            end if ! Axial

          end do ! photoionization
        end do ! For each atom

        ! Free
        if (allocated(integr)) deallocate(integr)
        if (allocated(inspect)) deallocate(inspect)
        if (allocated(buffer)) deallocate(buffer)
        if (allocated(dy)) deallocate(dy)
        if (allocated(dyv)) deallocate(dyv)
        if (allocated(b)) deallocate(b)
        if (allocated(c)) deallocate(c)
        if (allocated(d)) deallocate(d)

      end if ! If there is a valid spectrum

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                  !
      ! Finish off non-initialized lines !
      !                                  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! For each atom
      do ia=1,nA

        ! Thermal part of the Doppler width
        DwT = Atom(ia)%cDopp*sqrt(T)

        ! For each b-b transition
        do jtran=1,Atom(ia)%ntran

          ! Limits in frequency index of this transition
          lf0 = Atom(ia)%ilf0(jtran)
          lf1 = Atom(ia)%ilf1(jtran)

          ! If absent line, skip
          if (lf0.gt.lf1) cycle

          ! Check if input radiation. If yes, leave
          if (Atom(ia)%bbspecin(jtran)) cycle

          ! True limits
          jf0 = max(if0,lf0)
          jf1 = min(if1,lf1)

          ! If not possible, skip
          if (jf0.gt.jf1) cycle

          ! Boundary weights
          W0 = Atom(ia)%W0(jtran)
          W1 = Atom(ia)%W1(jtran)

          ! Flattened line index
          ktran = jtran + Atom(ia)%tshift

          !
          ! Flat
          !
          if (flat_cle_in) then

            ! Get frequency of transition
            Dfreq = Atom(ia)%Dfreq(jtran)

            ! Get CLV
            ljfreq = -1
            call getCLV(Dfreq,ljfreq,u1,u2)

            ! Get radiation from Allen
            if (use_allen) then

              ! Get Allen
              ljfreq = -1
              call getI(Dfreq,ljfreq,I0)

            ! Get radiation from Planck
            else

              ! Get Planck
              I0 = planck(Dfreq,Trad)

            end if

            ! Set value
            JKQ(0,0,ktran) = I0*0.5d0*(GeomP%CLV(1) + &
                                       GeomP%CLV(2)*u1 + &
                                       GeomP%CLV(3)*u2)
            JKQ(0,2,ktran) = I0*0.25d0*(GeomP%CLV(4) + &
                                        GeomP%CLV(5)*u1 + &
                                        GeomP%CLV(6)*u2)/sqrt(2d0)

          !
          ! No flat
          !
          else

            ! Identify involved terms
            itermu = Atom(ia)%fst(jtran)%itermu
            iterml = Atom(ia)%fst(jtran)%iterml

            ! Damping
            al = Atom(ia)%damp(iterml,1)
            au = Atom(ia)%damp(itermu,1)
            aul = Atom(ia)%ldamp(jtran,1)

            ! Initialize sum of weights
            Tfeta = 0d0

            ! For each FS transition
            do fjtran=1,Atom(ia)%fst(jtran)%nt

              ! Idenfity involved levels
              iJu = Atom(ia)%fst(jtran)%ilevelu(fjtran)
              iJl = Atom(ia)%fst(jtran)%ilevell(fjtran)

              ! Degeneration
              gl = 2d0*Atom(ia)%rJval(iJl,iterml) + 1d0

              ! Get the sequential index of this FS transition
              ffjtran = Atom(ia)%ifst_ij(fjtran,jtran)

              ! Get flattened line index
              ffktran = ffjtran + Atom(ia)%tfshift

              ! Get frequency of FS transition
              Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                      Atom(ia)%FSfreq(iJl,iterml)

              ! Add the microturbulence to Doppler width
              Dw = Dfreq*sqrt(DwT*DwT + vmi*vmi)

              ! Scale Dfreq to Dw
              Dfreq = Dfreq/Dw

              ! Absorptibity factor
              feta = Dfreq*Atom(ia)%fst(jtran)%Blu(iJl,iJu)*gl

              ! Add to total factor
              Tfeta = Tfeta + feta

              ! Add profile scale to feta
              feta = feta*1d-5*sqrt(IPI)/Dw

              ! Total damping
              at = (al + au + aul)/Dw

              ! Initialize CLV search index
              ljfreq = -1
              cljfreq = -1

              ! For each frequency
              do ifreq=jf0,jf1

                ! If first point
                if (ifreq.eq.lf0) then
                  WA = W0*feta
                ! Last point
                else if (ifreq.eq.lf1) then
                  WA = W1*feta
                ! Any other
                else
                  WA = Frec%W_freq(ifreq)*feta
                end if

                ! Get profile
                call voigtI(Dfreq - Frec%omega(ifreq),at,prof)

                ! Scale profile
                WA = WA*prof

                ! Get CLV
                call getCLV(Frec%omega(ifreq),ljfreq,u1,u2)

                ! Get radiation from Allen
                if (use_allen) then

                  ! Get Allen
                  call getI(Frec%omega(ifreq),cljfreq,I0)

                ! Get radiation from Planck
                else

                  ! Get Planck
                  I0 = planck(Frec%omega(ifreq),Trad)

                end if

                ! Contribution
                JKQ(0,0,ktran) = JKQ(0,0,ktran) + WA*I0* &
                                 0.5d0*(GeomP%CLV(1) + &
                                        GeomP%CLV(2)*u1 + &
                                        GeomP%CLV(3)*u2)
                JKQ(0,2,ktran) = JKQ(0,2,ktran) + WA*I0* &
                                0.25d0*(GeomP%CLV(4) + &
                                        GeomP%CLV(5)*u1 + &
                                        GeomP%CLV(6)*u2)/sqrt(2d0)
              end do ! Frequencies
            end do ! FS transitions

            ! Normalize
            Tfeta = 1d0/Tfeta
            JKQ(0,0:2,ktran) = JKQ(0,0:2,ktran)*Tfeta

          end if ! Flat input

        end do ! Transition

        ! If no photoionizations, skip
        if (Atom(ia)%nphot.le.0) cycle

        ! For each b-f transition
        do jtran=1,Atom(ia)%nphot

          ! Limits in frequency index of this transition
          lf0 = Atom(ia)%ipf0(jtran)
          lf1 = Atom(ia)%ipf1(jtran)

          ! If absent transition, leave
          if (Atom(ia)%phot(jtran)%absent) cycle

          ! Check if input radiation. If yes, leave
          if (Atom(ia)%bfspecin(jtran)) cycle

          ! True limits
          jf0 = max(if0,lf0)
          jf1 = min(if1,lf1)

          ! If not possible, skip
          if (jf0.gt.jf1) cycle

          ! Saha factor
          c3 = Saha*exp(Atom(ia)%phot(jtran)%edge*arg)* &
               Atom(ia)%phot(jtran)%glu

          ! Boundary weights
          W0 = Atom(ia)%phot(jtran)%W0
          W1 = Atom(ia)%phot(jtran)%W1

          ! Flattened line index
          ktran = jtran + Atom(ia)%pshift

          ! Initialize search frequency
          ljfreq = -1

          ! For each frequency
          do ifreq=jf0,jf1

            ! If first point
            if (ifreq.eq.lf0) then

              ! Get left weight
              c1 = W0

            ! Last point
            else if (ifreq.eq.lf1) then

              ! Get right weight
              c1 = W1

            ! Any other
            else

              ! Get weight
              c1 = Frec%W_freq(ifreq)

            end if ! If extreme or not

            ! Weight with cross section and
            ! constants
            c1 = c0*c1*Atom(ia)%phot(jtran)%alpha(ifreq)/ &
                 Frec%omega(ifreq)

            ! Get exu
            exu = c3*diexp(c0*Frec%omega(ifreq))

            ! Get radiation from Allen
            if (use_allen) then

              ! Get Allen
              call getI(Frec%omega(ifreq),cljfreq,I0)

            ! Get radiation from Planck
            else

              ! Get Planck
              I0 = planck(Frec%omega(ifreq),Trad)

            end if

            ! Geometry
            call getCLV(Frec%omega(ifreq),ljfreq,u1,u2)

            ! Weight
            WA = c1*Geom%W_mu(ith)*0.5d0*I0* &
                 (GeomP%CLV(1) + GeomP%CLV(2)*u1 + &
                  GeomP%CLV(3)*u2)

            ! Add contribution
            Jphot(ktran,1) = Jphot(ktran,1) + WA
            Jphot(ktran,2) = Jphot(ktran,2) + WA*exu

          end do ! Frequencies
        end do ! photoionization
      end do ! For each atom

      ! If doing MPI, we need to share the data
      if (nproc.gt.1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE,JKQ(-2,0,1),size(JKQ), &
                           MPI_DOUBLE_COMPLEX,MPI_SUM, &
                           MPI_COMM_RT,ios)
        call MPI_ALLREDUCE(MPI_IN_PLACE,JPhot(1,1),size(Jphot), &
                           MPI_DOUBLE_PRECISION,MPI_SUM, &
                           MPI_COMM_RT,ios)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !                                !
      ! Rotate the flat-spectrum lines !
      !                                !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! If magnetic field
      if (Bs.gt.TINYB) then

        ! For each atom
        do ia=1,nA

          ! For each b-b transition
          do jtran=1,Atom(ia)%ntran

            ! Check if input radiation. If yes, leave
            if (Atom(ia)%bbspecin(jtran)) cycle

            ! Flattened line index
            ktran = jtran + Atom(ia)%tshift

            ! Rotate
            call fieldB(JKQ(:,:,ktran),1,Flgsg,Bt,Bp,1)

          end do ! Transition
        end do ! For each atom

      end if ! There is magnetic field

      ! Free pointers
      nullify(T,vx,vy,vz,vmi)

      end subroutine getSEEJ

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get boundary condition for the radiation field from the input
      !! radiation\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          omega(double(:)): Frequency array\n
      !!                vx(double): Velocity vector along X\n
      !!                vy(double): Velocity vector along Y\n
      !!                vz(double): Velocity vector along Z\n
      !!  GeomP(Coronapoint_class): Structure with geometric data for
      !!                            a CLE node\n
      !!        spect(spect_class): Structure with the input spectra
      !!                            data\n
      !!       Stokes(double(:,:)): Stokes parameters
      subroutine get_bottom_spect(Atmo,omega,vx,vy,vz,GeomP, &
                                  spect,stokes)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(spect_class), intent(in):: spect
      type(Coronapoint_class), intent(in):: GeomP
      double precision, intent(in):: vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:), intent(out):: stokes

      ! Local

      logical:: first

      integer:: nf,nfs,nstk,istk
      integer:: ifreq,jfreq,cljfreq,ljfreq,if0,if1,jf0,jf1

      double precision:: mu_disk
      double precision:: ct,st,cc,sc,vfac,omegain,dx
      double precision, dimension(:), allocatable:: dy,b,c,d
      double precision, dimension(:,:), allocatable:: instokes

      ! If dynamic
      if (dyn) then

        ! Get angles
        ct = cos(GeomP%geom(1))
        st = sin(GeomP%geom(1))
        cc = cos(GeomP%geom(2))
        sc = sin(GeomP%geom(3))

        ! Get Doppler factor
        vfac = 1d0 + vz*ct + &
                       vx*st*cc + &
                       vy*st*sc
      ! Static
      else

        ! No Doppler factor
        vfac = 1d0

      end if

      ! Get dize of omega
      nf = size(omega)

      ! Size of spectra
      nstk = spect%nstk + 1

      ! Get the polar nodes on the disk
      mu_disk = sqrt(Atmo%ypos*Atmo%ypos + Atmo%zpos*Atmo%zpos)

      ! If overflow sine
      if (mu_disk.ge.1d0) then

        ! Cosine is 0
        mu_disk = 0d0

      ! If normal sine
      else

        ! Get cosine
        mu_disk = sqrt(1d0 - mu_disk*mu_disk)

      end if ! Control correct cosine

      ! Look for frequency limits
      jf0 = minloc(abs(omega(1)*vfac - spect%omega),1)
      jf1 = maxloc(abs(omega(nf)*vfac - spect%omega),1)

      ! Increasing order
      if (jf1.gt.jf0) then

        ! Keep order
        if0 = jf0
        if1 = jf1

      ! Decreasing order
      else

        ! Change order
        if0 = jf1
        if1 = jf0

      end if ! Order

      ! Size
      nfs = if1 - if0 + 1

      ! Allocate interpolation auxiliar
      allocate(dy(nstk))

      ! Allocate splines
      allocate(b(spect%nmu),c(spect%nmu),d(spect%nmu))

      ! Input Stokes for LOS
      allocate(instokes(nstk,if0:if1))

      ! If input is axial
      if (spect%axial) then

        ! For each polarization state
        do istk=1,nstk

          ! For each frequency in the input (split)
          do ifreq=if0,if1

            ! Get splines
            call spline(spect%mu,spect%stokes(:,1,ifreq,istk-1), &
                        b,c,d,spect%nmu)

            ! Get intensity
            instokes(istk,ifreq) = ispline(mu_disk,spect%mu, &
                                    spect%stokes(:,1,ifreq,istk-1), &
                                    b,c,d,spect%nmu)

          end do ! Frequencies
        end do ! Stokes parameters

      ! Input not axial
      else

        ! Not implemented!
        umsg = 'Non-axial input spectrum is not implemented'
        call aborted
        return

      end if ! Axial

      !
      ! Now interpolate in relevant frequencies
      !

      ! Last searched frequency
      cljfreq = if0

      ! Indicate it is first
      first = .True.

      ! For each frequency
      do ifreq=1,nf

        ! Update last jfreq
        ljfreq = cljfreq

        ! Get displaced source frequency
        omegain = omega(ifreq)*vfac

        ! If below lower limit
        if (omegain.le.spect%omega(if0)) then

          ! Constant extension
          stokes(1:nstk,ifreq) = instokes(:,if0)
          cycle

        ! Or above limits
        else if (omegain.ge.spect%omega(if1)) then

          ! Constant extension
          stokes(1:nstk,ifreq) = instokes(:,if1)
          cycle

        end if

        ! Look for the surrounding frequencies
        do jfreq=ljfreq,if1-1

          ! If below current limit, skip
          if (omegain.gt.spect%omega(jfreq).and. &
              omegain.lt.spect%omega(jfreq+1)) then

            ! If different than last time, calculate
            ! linear interpolation coefficients
            if (cljfreq.ne.ljfreq.or.first) then

              ! Get interpolation parameters
              dx = 1d0/(spect%omega(jfreq+1) - spect%omega(jfreq))
              dy = instokes(:,jfreq+1) - instokes(:,jfreq)

              ! Flag as not first
              if (first) first = .False.

            end if ! Different than last time or first

            ! Linear interpolation
            stokes(1:nstk,ifreq) = instokes(:,jfreq) + &
                                dy*dx*(omegain - spect%omega(jfreq))
            ! Finish
            exit

          ! If beyond next frequency, continue
          else if (omegain.gt.spect%omega(jfreq+1)) then

            ! Update current last searched frequency
            cljfreq = jfreq+1

            ! Continue searching
            cycle

          ! Exactly the same than current frequency
          else if (omegain.le.spect%omega(jfreq)) then

            ! No interpolation needed
            stokes(1:nstk,ifreq) = instokes(:,jfreq)

            ! Finish
            exit

          ! Exactly the same than next frequency
          else

            ! Update current last searched frequency
            cljfreq = jfreq + 1

            ! No interpolation needed
            stokes(1:nstk,ifreq) = instokes(:,jfreq+1)

            ! Finish
            exit

          end if ! Searching condition

        end do ! input spectrum frequencies
      end do ! CPU frequencies

      ! Free
      deallocate(dy,b,c,d,instokes)

      end subroutine get_bottom_spect

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get boundary condition for the radiation field from Allen's
      !! tabulation\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!          omega(double(:)): Frequency array\n
      !!                vx(double): Velocity vector along X\n
      !!                vy(double): Velocity vector along Y\n
      !!                vz(double): Velocity vector along Z\n
      !!  GeomP(Coronapoint_class): Structure with geometric data for
      !!                            a CLE node\n
      !!       Stokes(double(:,:)): Stokes parameters
      subroutine get_bottom_allen(Atmo,omega,vx,vy,vz,GeomP,stokes)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Coronapoint_class), intent(in):: GeomP
      double precision, intent(in):: vx,vy,vz
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:), intent(out):: stokes

      ! Local

      integer:: nf,ifreq,cljfreq,ljfreq

      double precision:: mu_disk,u1,u2,I0,ct,st,cc,sc,vfac,omegain


      ! Zero polarization
      stokes(2:4,:) = 0d0

      ! If dynamic
      if (dyn) then

        ! Get angles
        ct = cos(GeomP%geom(1))
        st = sin(GeomP%geom(1))
        cc = cos(GeomP%geom(2))
        sc = sin(GeomP%geom(3))

        ! Get Doppler factor
        vfac = 1d0 + vz*ct + &
                     vx*st*cc + &
                     vy*st*sc
      ! Static
      else

        ! No Doppler factor
        vfac = 1d0

      end if ! Dynamic

      ! Get dize of omega
      nf = size(omega)

      ! Get the polar nodes on the disk
      mu_disk = sqrt(Atmo%ypos*Atmo%ypos + Atmo%zpos*Atmo%zpos)

      ! If overflow sine
      if (mu_disk.ge.1d0) then

        ! Cosine is 0
        mu_disk = 0d0

      ! If normal sine
      else

        ! Get cosine
        mu_disk = sqrt(1d0 - mu_disk*mu_disk)

      end if ! Control correct cosine

      !
      ! Now interpolate in relevant frequencies
      !

      ! Last searched frequency
      ljfreq = -1
      cljfreq = -1

      ! For each frequency
      do ifreq=1,nf

        ! Get displaced source frequency
        omegain = omega(ifreq)*vfac

        ! Get Allen data
        call getCLV(omegain,ljfreq,u1,u2)
        call getI(omegain,cljfreq,I0)

        ! Save intensity
        stokes(1,ifreq) = I0*(1d0 - &
                              (1d0 - mu_disk)*u1 - &
                              (1d0 - mu_disk*mu_disk)*u2)

      end do ! Outputs frequencies

      end subroutine get_bottom_allen

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get radiation intensity from Allen's tabulation\n
      !!     freq(double): Frequency to interpolate into\n
      !!  ifreq0(integer): Frequency index to start search\n
      !!     Inte(double): Intensity
      subroutine getI(freq,ifreq0,Inte)

      ! I/O

      integer, intent(inout):: ifreq0
      double precision, intent(in):: freq
      double precision, intent(out):: Inte

      ! Parameters

      integer, parameter:: NN = 43

      ! Wavelength in microns
      double precision, dimension(NN), parameter:: &
        xi = (/ 0.2d0,0.22d0,0.24d0,0.26d0,0.28d0,0.3d0,0.32d0, &
                0.34d0,0.36d0,0.37d0,0.38d0,0.39d0,0.4d0,0.41d0, &
                0.42d0,0.43d0,0.44d0,0.45d0,0.46d0,0.48d0,0.5d0, &
                0.55d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.9d0,1d0, &
                1.1d0,1.2d0,1.4d0,1.6d0,1.8d0,2d0,2.5d0,3d0,4d0,5d0, &
                6d0,8d0,10d0,12d0 /) !mum

      ! Intensity in HanleRT units
      double precision, dimension(NN), parameter:: &
        yi = (/ 2.4d3,1.0164d4,1.6704d4,4.056d4,1.0192d5,2.205d5, &
                3.328d5,4.35812d5,5.35248d5,5.79087d5,6.68572d5, &
                7.52895d5,8.24d5,8.84206d5,9.31392d5,9.68876d5, &
                1.004784d6,1.03275d6,1.058d6,1.103616d6,1.1375d6, &
                1.21605d6,1.2672d6,1.29285d6,1.3181d6,1.2825d6, &
                1.2992d6,1.2717d6,1.26d6,1.2221d6,1.1664d6, &
                1.0388d6,9.216d5,7.7112d5,6.4d5,4.875d5,3.69d5, &
                2.272d5,1.55d5,1.152d5,6.08d4,3.5d4,2.592d4 /)

      ! Local

      integer:: ifreq

      double precision:: lamb,dxs,dx,dy


      ! Get lambda in microns
      lamb = 1d-1/freq

      ! Below lower limit
      if (lamb.le.xi(1)) then

        ifreq0 = 1
        Inte = yi(1)

      ! Above upper limit
      else if (lamb.ge.xi(NN)) then

        ifreq0 = NN
        Inte = yi(NN)

      ! In range
      else

        ! Negative index
        if (ifreq0.lt.0) then

          ! Search
          ifreq0 = minloc(abs(xi - lamb),1)

          ! If above
          if (xi(ifreq0).lt.lamb) ifreq0 = ifreq0 + 1

        end if

        ! For available frequencies
        do ifreq=ifreq0,2,-1

          ! If between two of them, interpolate
          if (lamb.gt.xi(ifreq-1).and.lamb.le.xi(ifreq)) then

            ! Indexing
            ifreq0 = ifreq

            ! Interpolation parameters
            dxs = lamb - xi(ifreq-1)
            dx = xi(ifreq) - xi(ifreq-1)
            dy = yi(ifreq) - yi(ifreq-1)

            ! Interpolation
            Inte = yi(ifreq-1) + dxs*dy/dx

            exit

          end if ! If interpolating

        end do ! Available frequencies

      end if ! In range?

      end subroutine getI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get CLV coefficients for a given frequency from Allen's
      !! tabulation\n
      !!     freq(double): Frequency to interpolate into\n
      !!  ifreq0(integer): Frequency index to start search\n
      !!       u1(double): First CLV coefficient\n
      !!       u2(double): Second CLV coefficient
      subroutine getCLV(freq,ifreq0,u1,u2)

      ! I/O

      integer, intent(inout):: ifreq0
      double precision, intent(in):: freq
      double precision, intent(out):: u1,u2

      ! Parameters

      ! Number of data nodes
      integer, parameter:: NN = 22

      ! Lambda in microns
      double precision, dimension(NN), parameter:: xi = &
        (/ 0.2d0,0.22d0,0.245d0,0.265d0,0.28d0,0.3d0,0.32d0,0.35d0, &
           0.37d0,0.38d0,0.4d0,0.45d0,0.5d0,0.55d0,0.6d0,0.8d0,1d0, &
           1.5d0,2d0,3d0,5d0,10d0 /)

      ! u1
      double precision, dimension(NN), parameter:: yi1 = &
        (/ 0.12d0,-1.30d0,-0.10d0,-0.10d0,0.38d0,0.74d0,0.88d0, &
           0.98d0,1.03d0,0.92d0,0.91d0,0.99d0,0.97d0,0.93d0,0.88d0, &
           0.73d0,0.64d0,0.57d0,0.48d0,0.35d0,0.22d0,0.15d0 /)

      ! u2
      double precision, dimension(NN), parameter:: yi2 = &
        (/ 0.33d0,1.6d0,0.85d0,0.90d0,0.57d0,0.20d0,0.03d0,-0.10d0, &
           -0.16d0,-0.05d0,-0.05d0,-0.17d0,-0.22d0,-0.23d0,-0.23d0, &
           -0.22d0,-0.20d0,-0.21d0,-0.18d0,-0.12d0,-0.07d0,-0.07d0 /)

      ! Local

      integer:: ifreq

      double precision:: lamb,dxs,dx,dy


      ! Get wavelength in microns
      lamb = 1d-1/freq

      ! Below lower limit
      if (lamb.le.xi(1)) then

        ! Extend as constant
        ifreq0 = 1
        u1 = yi1(1)
        u2 = yi2(1)

      ! Above upper limit
      else if (lamb.ge.xi(NN)) then

        ! Extend as constant
        ifreq0 = NN
        u1 = yi1(NN)
        u2 = yi2(NN)

      ! In range
      else

        ! Initialize if negative
        if (ifreq0.lt.0) ifreq0 = NN

        ! For available frequencies
        do ifreq=ifreq0,2,-1

          ! If between two of them, interpolate
          if (lamb.gt.xi(ifreq-1).and.lamb.le.xi(ifreq)) then

            ! Get last index
            ifreq0 = ifreq

            ! X axis interpolation parameters
            dxs = lamb - xi(ifreq-1)
            dx = xi(ifreq) - xi(ifreq-1)

            ! Interpolate u1
            dy = yi1(ifreq) - yi1(ifreq-1)
            u1 = yi1(ifreq-1) + dxs*dy/dx

            ! Interpolate u2
            dy = yi2(ifreq) - yi2(ifreq-1)
            u2 = yi2(ifreq-1) + dxs*dy/dx

            exit

          end if ! Interpolate

        end do ! Available frequencies

      end if ! If in range

      return

      end subroutine getCLV

!#####################################################################
!#####################################################################
!#####################################################################

      end module jcalccle_mod
