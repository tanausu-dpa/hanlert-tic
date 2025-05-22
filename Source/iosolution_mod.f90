      !> Solution reader and output writer
      module iosolution_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Contributors:
!     Hao Li (IAC/NSSCC)
!  Start:
!     20/04/2016
!  Last version:
!     15/05/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.3 - Generalized declarations of Atom to
!                             allow for empty arrays for any of
!                             them (TdPA)
!                           - Skip reading or writing atomic
!                             quantities if there are no active atoms
!                             in the run (TdPA)
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
!  readsol
!    Read a file with an existing self-consistent solution
!
!  getsol
!    Restore the self-consistent solution from RAM
!
!  setsol
!    Save the self-consistent solution in RAM
!
!  writesol
!    Save the self-consistent solution in a file
!
!  writesolI
!    Save the self-consistent solution for the intensity problem in a
!  file
!
!  writestk
!    Save the emergent Stokes parameters in a file
!
!  writestkI
!    Save the emergent intensity in a file
!
!  setstk
!    Save the emergent Stokes parameters in RAM
!
!  setstkI
!    Save the emergent intensity in RAM
!
!  write_CLEgeom
!    Write the position of the LOS in the CLE problem in a file
!
!  write_CLE
!    Write the emergent Stokes and optical depth in the CLE problem in
!  a file
!
!  writectr
!    Save the contribution function of a synthesis run in a file
!
!  writectr_inv
!    Save the contribution function of an inversion run in a file
!
!  writectrI
!    Save the intensity contribution function of a synthesis run in a
!  file
!
!  writectrI_inv
!    Save the intensity contribution function of an inversion run in a
!  file
!
!  setctr
!    Save the contribution function in RAM
!
!  setctrI
!    Save the intensity contribution function in RAM
!
!  writetau
!    Save the height of optical depth equal to one of a synthesis run
!  in a file
!
!  writetau_inv
!    Save the height of optical depth equal to one of an inversion run
!  in a file
!
!  settau
!    Save the height where the optical depth is equal to one in RAM
!
!  writecols
!    Save inelastic collisional rates of active atoms in a file
!
!  writedamp
!    Save damping parameter of active atoms in a file
!
!  writeqel
!    Save elastic collisional rates of active atoms in a file
!
!  writeback
!    Save background opacity, scattering coefficient, and emissivity
!  in a file
!
!  writeatmo
!    Write the model atmosphere in a file. The file is in ASCII for
!  a 1D synthesis and in binary for a 1.5D synthesis, with the same
!  format than the corresponding input for the latter
!
!  wAtmo
!    Write a 1D model atmosphere in a file in the 1D synthesis format
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use fieldb_mod
      use parameters_mod , only : pi , TINYB , TINYR , cZero , c , me
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read a file with an existing self-consistent solution\n
      !!     filename(character(:)): Name of the file to read\n
      !!      GeomI(Geometry_class): Structure with geometric data for
      !!                             the intensity problem\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       read_stokes(logical): If the frequency dependent
      !!                             radiation field could be read\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!   Stokes0(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates
      subroutine readsol(filename,GeomI,Geom,Flgsg,Bfield,Atom, &
                         read_stokes, &
                         Stokes,JKQ,JKQS,JKQC, &
                         Stokes0,J00,J00S,J00C,J00P)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      character(len=500), intent(in):: filename
      logical, intent(out):: read_stokes
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(out):: Stokes0
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(out):: Stokes
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(out):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(out):: J00P
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(out):: JKQC

      ! Local

      character(len=2):: label
      character(len=8):: scoord

      logical:: inaxial,instm,inAV

      integer:: ios,ia,iz,ifreq,i,ith,iph,idir,iS,itran,jtran
      integer:: it,iJ,iJ1,K,iQ,iR,iR0,iR1
      integer:: ia1,ia2,ia3,ia4,ia5
      integer:: ilabel,inaxial_int,instm_int,inAV_int
      integer:: rsize,jsize,psize,csize,ssize,nsize
      integer, dimension(:), allocatable:: ibuff

      double precision:: rJ,rJ1,da1,da2,rho0

      complex(kind=8):: integr
      complex(kind=8), dimension(-2:2,0:2,Rz0:Rz1):: JKQaux
      complex(kind=8), dimension(-nkx:nkx,Rz0:Rz1):: rhoKQaux


      ! Routine name
      urou = 'readsol'

      ! Master
      if (pid.eq.0) then

        !
        ! If 1D or inversion
        !
        if (run_mode.le.0) then

          ! Open file
          open (200,file=trim(filename), &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')
        !
        ! If 1.5D
        !
        else if (run_mode.eq.1) then

          ! Get LOS index
          write(scoord,'(I0.8)') icoords(3)

          ! Open file
          open (200,file=trim(filename)//'/Solution-'//scoord, &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='read', form='unformatted')

        end if ! 1D/1.5D
      end if ! Master


      !
      ! Dimensions, can be use for error handling
      !

      ! Master
      if (pid.eq.0) then

        ! Read label
        read (200,err=1100) label

        ! If polarization
        if (label.eq.'sp') then

          ! Flag as 0
          ilabel = 0

        ! If intensity
        else if (label.eq.'si') then

          ! Flag as 1
          ilabel = 1

        ! None
        else

          ! Flag as error
          ilabel = -1

        end if ! Label
      end if ! Master

      ! If MPI
      if (nproc.gt.1) then

        ! Control
        call control
        if (laborted) return

        ! Share label
        call MPI_BCAST(ilabel, 1, MPI_INTEGER, 0, MPI_COMM_RT, ierr)

      end if ! MPI

      ! Important to check the file label
      if(ilabel.ne.0.and.ilabel.ne.1) then

        ! Issue error
        umsg = 'The specified solution file does '// &
               'not have the correct ID of a solution file'
        call aborted
        return

      end if ! Wrong label

      ! Allocate J00 for photoionizations (common for both pathts)
      if (nA.gt.0) then
        allocate(J00P(nxphot,2,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00P)
      end if

      !
      ! Polarization read
      !
      if (ilabel.eq.0) then

        ! Initialize J00P
        if (nA.gt.0) J00P = 0d0

        !
        ! Allocations
        !

        ! If keeping Stokes parameters
        if (KSTK) then

          ! Allocate full height range
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1

        ! Not keeping Stokes parameters
        else

          ! Allocate two height nodes
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1

        end if ! Keep Stokes

        ! Count memory
        RRAMc = RRAMc + 1d-6*sizeof(Stokes)

        ! If atoms
        if (nA.gt.0) then

          ! JKQ for absorptivity
          allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(JKQ)

          ! JKQ for stimulated emission
          allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(JKQS)

        end if ! Atoms

        ! JKQ frequency dependent
        allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(JKQC)

        ! If the master
        if (pid.eq.0) then

          ! Read solution metadata
          read (200,err=1100) ia1
          read (200,err=1100) ia2
          read (200,err=1100) ia3
          read (200,err=1100) ia4
          read (200,err=1100) ia5
          read (200,err=1100) inaxial_int
          read (200,err=1100) instm_int
          read (200,err=1100) inAV_int

        end if ! Master

        ! If doing MPI
        if (nproc.gt.1) then

          ! Allocate transfer buffer
          allocate(ibuff(8))

          ! Control
          call control
          if (laborted) return

          ! If the master
          if (pid.eq.0) then

            ! Store metadata in buffer
            ibuff(1) = ia1
            ibuff(2) = ia2
            ibuff(3) = ia3
            ibuff(4) = ia4
            ibuff(5) = ia5
            ibuff(6) = inaxial_int
            ibuff(7) = instm_int
            ibuff(8) = inAV_int

          end if ! Master

          ! Share metadata
          call MPI_BCAST(ibuff(1), 8, MPI_INTEGER, 0, &
                         MPI_COMM_RT, ierr)

          ! Slave
          if (pid.ne.0) then

            ! Extract metadata
            ia1 = ibuff(1)
            ia2 = ibuff(2)
            ia3 = ibuff(3)
            ia4 = ibuff(4)
            ia5 = ibuff(5)
            inaxial_int = ibuff(6)
            instm_int = ibuff(7)
            inAV_int = ibuff(8)

          end if ! slave
        end if ! MPI

        ! Convert to logical
        inaxial = inaxial_int.eq.1

        ! Convert to logical
        instm = instm_int.eq.1

        ! Convert to logical
        inAV = inAV_int.eq.1

        ! Flag to read stokes parameter, initialize true
        read_stokes = .True.

        !
        ! Dimension checking
        !

        ! Check height nodes
        if (ia2.ne.nz) then

          ! Issue error
          umsg = 'Solution file with different number of '// &
                 'heights.'
          call aborted
          return

        end if ! Wrong height nodes

        ! Check number of atoms
        if (ia5.ne.nA) then

          ! Issue error
          umsg = 'Solution file with different number of '// &
                 'atoms.'
          call aborted
          return

        end if ! Wrong number of atoms

        ! Check number of frequencies, this one does not produce an
        ! abortion
        if (ia1.ne.nfreq) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of frequencies in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes and J^K_Q(nu).'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of frequencies

        ! Check polar nodes, this one does not produce an abortion
        if (ia3.ne.Geom%nTh.and.read_stokes.and..not.inAV) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of polar directions in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of polar nodes

        ! Check azimuthal nodes, this one does not produce an abortion
        if (ia4.ne.Geom%nPh.and.(ia4.ne.1.and.Geom%nPh.ne.1).and. &
            read_stokes.and..not.inAV) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'and is not axial, ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of azimuthal nodes

        ! Check azimuthal nodes for AD redistribution, this one does
        ! not produce an abortion
        if (ia4.ne.Geom%nPh2.and.(ia4.ne.1).and. &
            read_stokes.and..not.inAV) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'for eps^(2) and is not axial, ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of azimuthal nodes for AD PRD

        ! Warning when reading non-axial from AD for axial AD
        if (.not.inaxial.and.axial.and.(.not.(AV.and..not.dyn)).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then

          ! Issue warning
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read.'
          call verbose

        end if ! Reading non-axial for axial

        ! Warning when reading non-axial from AD for axial AV
        if (.not.inaxial.and.axial.and.(AV.and..not.dyn).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then

          ! Issue warning
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read to compute JKQC.'
          call verbose

        end if ! Reading non-axial AD for axial AV


        !
        ! Population and RhoKQ
        !

        ! For each atom
        do ia=1,nA

          ! Only the master reads
          if (pid.eq.0) then

            ! For each height
            do iz=1,nZ

              ! Read population at this height
              read (200,err=1100) da1

              ! If synthesis, get atomic population
              if (run_mode.ge.0) Atom(ia)%n(iz) = da1

            end do ! heights

            ! For each term
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! For each level
                do iJ1=1,Atom(ia)%nJ(it)

                  ! Get J'
                  rJ1 = Atom(ia)%rJval(iJ1,it)

                  ! For each K
                  do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                    ! For each Q
                    do iQ=-K,K

                      ! Get rho(J,J')KQ index
                      if (K.le.Atom(ia)%Kcut(it)) &
                        iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! For each height
                      do iz=1,nZ

                        ! Read real and imaginary parts
                        read (200,err=1100) da1,da2

                        ! If out of considered height range, skip
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! If not neglecting this multipole
                        if (K.le.Atom(ia)%Kcut(it)) then

                          ! Save
                          Atom(ia)%crho(iR,iz) = dcmplx(da1,da2)

                          ! Auxiliar variable for rotation
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end if ! Not neglecting this multipole

                      end do ! heights
                    end do ! Q

                    !
                    ! Rotate rhoKQ
                    !

                    ! Only if above the limit
                    if (K.le.Atom(ia)%Kcut(it)) then

                      ! For each height in the CPU domain
                      do iz=Rz0,Rz1

                        ! If there is magnetic field
                        if (Bfield%Bstrength(iz).gt.TINYB) then

                          ! Rotate the rhoKQ in the auxiliar variable
                          call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                    Bfield%Btheta(iz), &
                                    Bfield%Bphi(iz),1)

                          ! For each Q
                          do iQ=-K,K

                            ! Get index
                            iR = Atom(ia)%irho(it)% &
                                          Jrho(iJ1,iJ)%kq(iQ,K)

                            ! Store the rotated result in the rhoKQ
                            ! array
                            Atom(ia)%crho(iR,iz) = rhoKQaux(iQ,iz)

                          end do ! Q

                        end if ! There is B field

                      end do ! heights

                    end if ! K < Kcut

                  end do ! K
                end do ! J'
              end do ! J
            end do ! terms

          end if ! Master


          ! If there are slaves
          if (nproc.gt.1) then

            ! Control
            call control
            if (laborted) return

            ! Get buffer sizes
            nsize = nZ
            rsize = Atom(ia)%ndim*RnZ

            ! Send n
            call MPI_BCAST(Atom(ia)%n(1), nsize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Send rho
            call MPI_BCAST(Atom(ia)%crho(1,Rz0), rsize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

          end if ! MPI

        end do ! atoms

        !
        ! Flag the null rho(J,J')KQ
        !

        ! For each atom
        do ia=1,nA

          ! Initialize to non-null
          Atom(ia)%rhonull = .False.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)

              ! Get the rho00 indexes
              iR0 = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! For each level
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J'
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! Get the rho00 indexes
                iR1 = Atom(ia)%irho(it)%Jrho(iJ1,iJ1)%kq(0,0)

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(it))

                  ! For each Q
                  do iQ=-K,K

                    ! Get rho(J,J')KQ index
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! For each height in this domain
                    do iz=Rz0,Rz1

                      ! Get the inverse of rho00
                      rho0 = 1d0/sqrt(abs(Atom(ia)%crho(iR0,iz))* &
                                      abs(Atom(ia)%crho(iR1,iz)))

                      ! If rhoKQ/rho00 is lesser than double precision
                      ! flag null
                      if (abs(Atom(ia)%crho(iR,iz)*rho0).lt.TINYR) &
                        Atom(ia)%rhonull(iR,iz) = .True.

                    end do ! heights
                  end do ! Q
                end do ! K
              end do ! J'
            end do ! J
          end do ! terms
        end do ! atoms


        !
        ! JKQ
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Atomic shift
              jtran = itran + Atom(ia)%tshift

              ! For each K
              do K=0,2

                ! In K cut
                if (K.le.Atom(ia)%Krad(itran)) then

                  ! For each Q
                  do iQ=-K,K

                    ! For each height
                    do iz=1,nZ

                      ! Read real and imaginary parts of JKQ
                      read (200,err=1100) da1,da2

                      ! If out of considered height range, skip
                      if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                      ! Save in auxiliary
                      JKQaux(iQ,K,iz) = dcmplx(da1,da2)

                    end do ! heights
                  end do ! Q

                ! Out of cut
                else

                  ! Jump data
                  iQ = (2*K+1)*nZ*16
                  call fseek(200,iQ,1)

                  ! Set to zero
                  JKQaux(:,K,:) = cZero

                end if ! K cut
              end do ! K

              !
              ! Rotate JKQ
              !

              ! For each height
              do iz=Rz0,Rz1

                ! Register JKQ in the array
                JKQ(:,:,jtran,iz) = JKQaux(:,:,iz)

                ! If there is magnetic field, rotate
                if (Bfield%Bstrength(iz).gt.TINYB) &
                  call fieldB(JKQ(:,:,jtran,iz),1,Flgsg, &
                              Bfield%Btheta(iz),Bfield%Bphi(iz),1)

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (nproc.gt.1.and.nA.gt.0) then

          ! Control
          call control
          if (laborted) return

          ! Buffer size
          jsize = 15*nxtran*RnZ

          ! Share JKQ
          call MPI_BCAST(JKQ(-2,0,1,Rz0), jsize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

        end if ! MPI

        ! If there is stimulated emission in the input
        if(instm)then

          ! Only the master reads
          if (pid.eq.0) then

            ! For each atom
            do ia=1,nA

              ! For each transition
              do itran=1,Atom(ia)%ntran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tshift

                ! For each K
                do K=0,2

                  ! In K cut
                  if (K.le.Atom(ia)%Krad(itran)) then

                    ! For each Q
                    do iQ=-K,K

                      ! For each height
                      do iz=1,nZ

                        ! Read real and imaginary parts of JKQS
                        read (200,err=1100) da1,da2

                        ! If out of considered height range, skip
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! If stimulated emission, store
                        if (stm) JKQaux(iQ,K,iz) = dcmplx(da1,da2)

                      end do ! heights
                    end do ! Q

                  ! Out of cut
                  else

                    ! Jump
                    iQ = (2*K+1)*nZ*16
                    call fseek(200,iQ,1)

                    ! Set to zero
                    JKQaux(:,K,:) = cZero

                  end if ! K cut

                end do ! K

                !
                ! Rotate JKQS
                !

                ! If we are currently doing stimulated emission
                if(stm)then

                  ! For each height
                  do iz=Rz0,Rz1

                    ! Register JKQS in the array
                    JKQS(:,:,jtran,iz) = JKQaux(:,:,iz)

                    ! If there is magnetic field, rotate
                    if (Bfield%Bstrength(iz).gt.TINYB) &
                      call fieldB(JKQS(:,:,jtran,iz),1,Flgsg, &
                                  Bfield%Btheta(iz),Bfield%Bphi(iz),1)

                  end do ! heights

                end if ! stimulated emission

              end do ! transitions
            end do ! atoms

          end if ! Master

          ! If there are slaves
          if (nproc.gt.1.and.nA.gt.0) then

            ! Control
            call control
            if (laborted) return

            ! Buffer size
            jsize = 15*nxtran*RnZ

            ! Share
            call MPI_BCAST(JKQS(-2,0,1,Rz0), jsize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

          end if ! MPI

        ! No stimulated emission in the input
        else

          ! Just copy JKQ data then
          if (stm.and.nA.gt.0) JKQS = JKQ

        end if ! Stimulated emission in the input


        !
        ! Stokes or JKQC
        !

        ! If we can read Stokes
        if (read_stokes) then

          ! If input is only JKQ
          if (inAV) then

            ! Only the master reads
            if (pid.eq.0) then

              ! For each height
              do iz=1,nz

                ! For each frequency
                do ifreq=1,nfreq

                  ! For each K
                  do K=0,2

                    ! If within K cut
                    if (K.le.Krad) then

                      ! For each Q
                      do iQ=-K,K

                        ! Read real and imaginary parts of JKQS
                        read (200,err=1100) da1,da2

                        ! Skip if out of range
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! Store
                        JKQC(iQ,K,ifreq,iz) = dcmplx(da1,da2)

                      end do ! Q

                    ! Out of cut
                    else

                      ! Jump
                      iQ = (2*K+1)*16
                      call fseek(200,iQ,1)

                      ! Put to zero
                      JKQC(:,K,ifreq,iz) = cZero

                    end if ! K cut

                  end do ! K
                end do ! frequency
              end do ! heights

            end if ! Master

            ! If there are slaves
            if (nproc.gt.1) then

              ! Control
              call control
              if (laborted) return

              ! Buffer size
              csize = 15*nfreq*RnZ

              ! Share
              call MPI_BCAST(JKQC(-2,0,1,Rz0), csize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! MPI

            ! If we are doing angle dependent, the first iteration
            ! must be angle averaged if the input is angle average
            if (.not.AV) tbAD = .True.

          ! If input is AD
          else

            ! If Master
            if (pid.eq.0) then

              !
              ! Currently doing AV or not PRD, no need to read Stokes
              ! Therefore, no need to keep Stokes
              !
              if (.not.KSTK.or..not.(dyn.or..not.AV)) then

                ! Initialize
                JKQC = cZero

                ! For each height
                do iz=1,nZ

                  ! Initialize directions
                  idir = 0

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! Axial input
                    if (inaxial) then

                      ! Add one direction
                      idir = idir + 1

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Integral
                        integr = cZero

                        ! For each Stokes parameter
                        do iS=0,3

                          ! Read the data point
                          read (200,err=1100) da1

                          ! Skip out of limits
                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          ! Store
                          if (KSTK) Stokes(iS,ifreq,1,ith,iz) = da1

                          ! For each K
                          do K=0,Krad

                            ! Sum over Stokes parameters of Stk*TKQ
                            integr = da1*Geom%TS(iS,0,K,idir)

                            ! Add contribution to the JKQC integral
                            JKQC(0,K,ifreq,iz) = &
                                               JKQC(0,K,ifreq,iz) + &
                                               integr*Geom%W_mu(ith)

                          end do ! K
                        end do ! Stokes
                      end do ! frequencies

                      ! If storing, and in limits
                      if (KSTK.and.iz.ge.Rz0.and.iz.le.Rz1) then

                        ! Fill other azimuths
                        do iph=2,Geom%nph
                          Stokes(:,:,iph,ith,iz) = &
                                                  Stokes(:,:,1,ith,iz)
                        end do

                      end if ! Storing and in limits

                      ! Remove one direction and add azimuth size
                      idir = idir - 1 + Geom%nPh2

                    ! Not axial input
                    else

                      ! For each azimuthal direction
                      do iph=1,ia4

                        ! Advance direction
                        idir = idir + 1

                        ! For each frequency
                        do ifreq=1,nfreq

                          ! For each Stokes parameter
                          do iS=0,3

                            ! Read the data point
                            read (200,err=1100) da1

                            ! Skip out of limits
                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            ! Axial?
                            if (axial.and.iph.gt.1) cycle

                            ! Store
                            if (KSTK) &
                              Stokes(iS,ifreq,iph,ith,iz) = da1

                            ! If axial
                            if (axial) then

                              ! For each K
                              do K=0,Krad

                                ! Sum over Stk parameters of Stk*TKQ
                                integr = da1*Geom%TS(iS,0,K,idir)

                                ! Add contribution to the JKQC int.
                                JKQC(0,K,ifreq,iz) = &
                                                JKQC(0,K,ifreq,iz) + &
                                                integr*Geom%W_mu(ith)
                              end do

                            ! If not axial
                            else

                              ! For each K
                              do K=0,Krad

                                ! For each Q
                                do iQ=0,K

                                  ! Sum over Stk parameters of Stk*TKQ
                                  integr = da1* &
                                           Geom%TS(iS,iQ,K,idir)

                                  ! Add contribution to the JKQC
                                  JKQC(iQ,K,ifreq,iz) = &
                                               JKQC(iQ,K,ifreq,iz) + &
                                               integr* &
                                               Geom%W_mu(ith)* &
                                               Geom%W_mux(iph)

                                end do ! Q
                              end do ! K

                            end if ! Axial

                          end do ! Stokes
                        end do ! frequencies
                      end do ! azimuthal nodes

                      ! If axial, advance direction index
                      if (axial) &
                        idir = idir - 1 + Geom%nPh2

                    end if ! In axial

                  end do ! polar nodes
                end do ! heights

                ! Complete if not axial
                if (.not.axial) then

                  ! K
                  do K=1,Krad

                    ! Q
                    do iQ=1,K

                      ! Conjugation properties
                      JKQC(-iQ,K,:,Rz0:Rz1) = Flgsg%sg(iQ)* &
                                           conjg(JKQC(iQ,K,:,Rz0:Rz1))

                    end do ! Q
                  end do ! K

                end if ! Non-axial

              !
              ! Currently doing AD
              !
              else

                ! Initialize
                JKQC = cZero

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,Geom%nTh

                    ! If input axial
                    if (inaxial) then

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! For each Stokes parameter
                        do iS=0,3

                          ! Read the data
                          read (200,err=1100) da1

                          ! Skip out of limits
                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          ! Save
                          Stokes(iS,ifreq,1,ith,iz) = da1

                        end do ! Stokes parameters
                      end do ! frequencies

                      ! Skip out of limits
                      if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                      ! Rest of azimuths
                      do iph=2,Geom%nph

                        ! Fill
                        Stokes(:,:,iph,ith,iz) = Stokes(:,:,1,ith,iz)

                      end do ! azimuthal nodes

                    ! Not input axial
                    else

                      ! For each azimuthal direction
                      do iph=1,ia4

                        ! For each frequency
                        do ifreq=1,nfreq

                          ! For each Stokes parameter
                          do iS=0,3

                            ! Read the data
                            read (200,err=1100) da1

                            ! Skip out of limits
                            if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                            ! Axial?
                            if (axial.and.iph.gt.1) cycle

                            ! Save
                            Stokes(iS,ifreq,iph,ith,iz) = da1

                          end do ! Stokes parameters
                        end do ! frequencies
                      end do ! azimuthal nodes

                    end if ! Input axial

                  end do ! polar nodes

                  ! Skip out of limits
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! For each K
                    do K=0,Krad

                      ! For each Q
                      do iQ=0,K

                        ! Initialize direction
                        idir = 0

                        ! For each polar direction
                        do ith=1,Geom%nTh

                          ! For each azimuthal direction
                          do iph=1,Geom%nPh

                            ! Advance direction
                            idir = idir + 1

                            ! Sum over Stokes parameters of Stokes*TKQ
                            integr = sum(Stokes(:,ifreq,iph,ith,iz)* &
                                         Geom%TS(:,iQ,K,idir))

                            ! Add contribution to the JKQC integral
                            JKQC(iQ,K,ifreq,iz) =  &
                                               JKQC(iQ,K,ifreq,iz) + &
                                 integr*Geom%W_mu(ith)*Geom%W_mux(iph)

                          end do ! azimuthal nodes

                          ! If axial, advance more
                          if (axial) &
                            idir = idir - 1 + Geom%nPh2

                        end do ! polar nodes
                      end do ! Q
                    end do ! K
                  end do ! frequencies
                end do ! heights

                !
                ! Complete rest of multi-poles
                !

                ! K
                do K=1,Krad

                  ! Q
                  do iQ=1,K

                    ! Conjugation properties
                    JKQC(-iQ,K,:,Rz0:Rz1) = Flgsg%sg(iQ)* &
                                           conjg(JKQC(iQ,K,:,Rz0:Rz1))

                  end do ! Q
                end do ! K

                !
                ! Q!=0 JKQC
                !

                ! If the input was axial or we are doing axial
                if(Geom%axial.or.inaxial)then

                  ! Kill the Q!=0 components
                  JKQC(-2:-1,1:2,:,:) = cZero
                  JKQC(1:2,1:2,:,:) = cZero

                endif ! axiality
              end if ! Read AD or AV
            end if ! Master

            ! If there are slaves
            if (nproc.gt.1) then

              ! Control
              call control
              if (laborted) return

              ! Slaves require Stokes
              if (PRD.and.(dyn.or..not.AV)) then

                ! Buffer size
                ssize = 4*nfreq*Geom%nTh*Geom%nPh*Rnz

                ! Share
                call MPI_BCAST(Stokes(0,1,1,1,Rz0), ssize, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)

              end if ! Slaves require Stokes

              ! Buffer size
              csize = 15*nfreq*RnZ

              ! Share
              call MPI_BCAST(JKQC(-2,0,1,Rz0), csize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! There are slaves
          end if ! Input AV or AD
        end if ! Can read Stokes

      !
      ! INTENSITY READ
      !
      else if (ilabel.eq.1) then

        !
        ! Allocations
        !

        ! If keeping Stokes parameters
        if (KSTK) then

          ! Allocate whole height axis
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1

        ! Not keeping Stokes
        else

          ! Allocate only two height nodes
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1

        end if ! Keeping Stokes

        ! Memory count
        RRAMc = RRAMc + 1d-6*sizeof(Stokes0)

        ! Initialize
        Stokes0 = 0d0

        ! If atoms
        if (nA.gt.0) then

          ! J00 for absorptivity
          allocate(J00(nxt,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00)

          ! J00 for stimulated emission
          allocate(J00S(nxt,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00S)

        end if ! Atoms

        ! J00 frequency dependent
        allocate(J00C(nfreq,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00C)

        ! If the master
        if (pid.eq.0) then

          ! Read metadata
          read (200,err=1100) ia1
          read (200,err=1100) ia2
          read (200,err=1100) ia3
          read (200,err=1100) ia4
          read (200,err=1100) ia5
          read (200,err=1100) inaxial_int
          read (200,err=1100) instm_int
          read (200,err=1100) inAV_int

        end if ! Master

        ! If doing MPI
        if (nproc.gt.1) then

          ! Allocate buffer
          allocate(ibuff(8))

          ! Control
          call control
          if (laborted) return

          ! Master
          if (pid.eq.0) then

            ! Save metadata in buffer
            ibuff(1) = ia1
            ibuff(2) = ia2
            ibuff(3) = ia3
            ibuff(4) = ia4
            ibuff(5) = ia5
            ibuff(6) = inaxial_int
            ibuff(7) = instm_int
            ibuff(8) = inAV_int

          end if ! Master

          ! Share
          call MPI_BCAST(ibuff(1), 8, MPI_INTEGER, 0, &
                         MPI_COMM_RT, ierr)

          ! If a slave
          if (pid.ne.0) then

            ! Get metadata from buffer
            ia1 = ibuff(1)
            ia2 = ibuff(2)
            ia3 = ibuff(3)
            ia4 = ibuff(4)
            ia5 = ibuff(5)
            inaxial_int = ibuff(6)
            instm_int = ibuff(7)
            inAV_int = ibuff(8)

          end if ! slave
        end if ! MPI

        ! Convert to logical
        inaxial = inaxial_int.eq.1

        ! Convert to logical
        instm = instm_int.eq.1

        ! Convert to logical
        inAV = inAV_int.eq.1

        ! Flag to read stokes parameter, initialize true
        read_stokes = .True.

        !
        ! Dimension checking
        !

        ! Check height nodes
        if (ia2.ne.nz) then

          ! Issue error
          umsg = 'Solution file with different number of '// &
                 'heights.'
          call aborted
          return

        end if ! Wrong height nodes

        ! Check number of atoms
        if (ia5.ne.nA) then

          ! Issue error
          umsg = 'Solution file with different number of '// &
                 'atoms.'
          call aborted
          return

        end if ! Wrong number of atoms

        ! Check number of frequencies, this one does not produce an
        ! abortion
        if (ia1.ne.nfreq) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of frequencies in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes and J^K_Q(nu).'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of frequencies

        ! Check polar nodes, this one does not produce an abortion
        if (ia3.ne.GeomI%nTh.and.read_stokes) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of polar directions in '// &
                   'solution file different than in system; '// &
                   'ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of polar nodes

        ! Check azimuthal nodes, this one does not produce an abortion
        if (ia4.ne.GeomI%nPh.and.(ia4.ne.1.and.GeomI%nPh.ne.1).and. &
            read_stokes.and..not.inAV) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'and is not axial, ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of azimuthal nodes

        ! Check azimuthal nodes for AD redistribution, this one does
        ! not produce an abortion
        if (ia4.ne.GeomI%nPh2.and.(ia4.ne.1).and. &
            read_stokes.and..not.inAV) then

          ! Master
          if (pid.eq.0) then

            ! Issue warning
            umsg = ' - Warning: Number of azimuths in '// &
                   'solution file different than in system '// &
                   'for eps^(2) and is not axial, ignoring Stokes.'
            call verbose

          end if ! Master

          ! Cannot read frequency dependent data
          read_stokes = .False.

        end if ! Wrong number of azimuthal nodes for AD PRD

        ! Warning when reading non axial from AD for axial AD
        if (.not.inaxial.and.axiali.and.(.not.(AVI.and..not.dyn)) &
            .and..not.inAV.and.read_stokes.and.pid.eq.0) then

          ! Issue warning
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read.'
          call verbose

        end if ! Reading non-axial for axial

        ! Warning when reading non axial from AD for axial AV
        if (.not.inaxial.and.axiali.and.(AVI.and..not.dyn).and. &
            .not.inAV.and.read_stokes.and.pid.eq.0) then

          ! Issue warning
          umsg = ' - Warning: The Solution contains '// &
                 'non-axial Stokes parameters, but we are ' // &
                 'assuming axial symmetry. Only the first '// &
                 'azimuth will be read to compute J00C.'
          call verbose

        end if ! Reading non-axial AD for axial AV

        !
        ! Population and Rho00
        !

        ! For each atom
        do ia=1,nA

          ! Only the master reads
          if (pid.eq.0) then

            ! Initialize pop
            Atom(ia)%crho = cZero

            ! For each height
            do iz=1,nZ

              ! Read population at this height
              read (200,err=1100) da1

              ! If synthesis, save in atomic population
              if (run_mode.ge.0) Atom(ia)%n(iz) = da1

            end do ! heights

            ! For each terms
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! Get component index
                iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                ! For each height
                do iz=1,nZ

                  ! Read real and imaginary parts
                  read (200,err=1100) da1

                  ! Skip if out of considered range
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                  ! Save
                  Atom(ia)%crho(iR,iz) = dcmplx(da1,0d0)

                end do ! heights
              end do ! J
            end do ! terms

          end if ! Master

          ! If there are slaves
          if (nproc.gt.1) then

            ! Control
            call control
            if (laborted) return

            ! Buffer sizes
            nsize = nZ
            rsize = Atom(ia)%ndim*RnZ

            ! Share n
            call MPI_BCAST(Atom(ia)%n(1), nsize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! Share rho00
            call MPI_BCAST(Atom(ia)%crho(1,Rz0), rsize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

          end if ! MPI

        end do ! atoms

        !
        ! Flag the null rho(J,J')KQ
        !

        ! For each atom
        do ia=1,nA

          ! Initialize to non-null
          Atom(ia)%rhonull = .True.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)
              rJ = sqrt(2d0*rJ + 1d0)

              ! Get level index
              i = Atom(ia)%irho(it)%irho_ij(iJ)

              ! Get component index
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              do iz=Rz0,Rz1

                ! If rhoKQ/rho00 is lesser than double precision
                ! flag null
                if (abs(Atom(ia)%crho(iR,iz)).gt.0d0) &
                  Atom(ia)%rhonull(iR,iz) = .False.

                ! Define population
                Atom(ia)%popu(i,iz) = dble(Atom(ia)%crho(iR,iz))*rJ

              end do ! heights
            end do ! J
          end do ! terms
        end do ! atoms


        !
        ! J00
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%nftran

              ! Apply atomic index
              jtran = itran + Atom(ia)%tfshift

              ! For each height
              do iz=1,nZ

                ! Read real and imaginary parts of JKQ
                read (200,err=1100) da1

                ! Skip if out of considered range
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                ! Save
                J00(jtran,iz) = da1

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (nproc.gt.1.and.nA.gt.0) then

          ! Control
          call control
          if (laborted) return

          ! Buffer size
          jsize = nxt*RnZ

          ! Share
          call MPI_BCAST(J00(1,Rz0), jsize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, &
                         ierr)
        end if ! MPI

        ! If there is stimulated emission in the input
        if (instm) then

          ! Only the master reads
          if (pid.eq.0) then

            ! For each atom
            do ia=1,nA

              ! For each transition
              do itran=1,Atom(ia)%nftran

                ! Apply atomic shift
                jtran = itran + Atom(ia)%tfshift

                ! For each height
                do iz=1,nZ

                  ! Read real and imaginary parts of JKQ
                  read (200,err=1100) da1

                  ! Skip if out of considered range
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                  ! If stimulated emission, save
                  if (stm) J00S(jtran,iz) = da1

                end do ! heights
              end do ! transitions
            end do ! atoms

          end if ! Master

          ! If there are slaves
          if (nproc.gt.1.and.nA.gt.0) then

            ! Control
            call control

            ! Buffer size
            jsize = nxt*RnZ

            ! Share
            call MPI_BCAST(J00S(1,Rz0), jsize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          end if ! MPI

        ! No stimulated emission in the input
        else

          ! Just copy J00 if necessary
          if (stm.and.nA.gt.0) J00S = J00

        end if ! Stimulated emission in the input


        !
        ! J00P
        !

        ! Only the master reads
        if (pid.eq.0) then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%nphot

              ! Apply atomic shift
              jtran = itran + Atom(ia)%pshift

              ! For each height
              do iz=1,nZ

                ! Read photoionization rates
                read (200,err=1100) da1
                read (200,err=1100) da2

                ! If out of considered range, skip
                if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                ! Save
                J00P(jtran,1,iz) = da1
                J00P(jtran,2,iz) = da2

              end do ! heights
            end do ! transitions
          end do ! atoms

        end if ! Master

        ! If there are slaves
        if (nproc.gt.1.and.nA.gt.0) then

          ! Control
          call control
          if (laborted) return

          ! Buffer size
          psize = nxphot*2*RnZ

          ! Share
          call MPI_BCAST(J00P(1,1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

        end if ! MPI


        !
        ! Stokes or J00C
        !

        ! If we can read Stokes
        if (read_stokes) then

          ! If input is AV
          if (inAV) then

            ! Only the master reads
            if (pid.eq.0) then

              ! For each height
              do iz=1,nz

                ! For each frequency
                do ifreq=1,nfreq

                  ! Read real and imaginary parts of JKQS
                  read (200,err=1100) da1

                  ! Skip if out of considered range
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                  ! Save
                  J00C(ifreq,iz) = da1

                end do ! frequency
              end do ! heights

            end if ! Master

            ! If there are slaves
            if (nproc.gt.1) then

              ! Control
              call control
              if (laborted) return

              ! Buffer size
              csize = nfreq*RnZ

              ! Share
              call MPI_BCAST(J00C(1,Rz0), csize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

            end if ! MPI

            ! If we are doing angle dependent, the first iteration
            ! must be angle averaged if the input is angle average
            if (.not.AV) tbAD = .True.

          ! If input is AD
          else

            ! Master
            if (pid.eq.0) then

              !
              ! Currently doing AV without velocities, no need to
              ! read Stokes
              !
              if (.not.KSTK.or..not.(dyn.or..not.AV)) then

                ! Initialize
                J00C = 0d0

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,GeomI%nTh

                    ! Axial input
                    if (inaxial) then

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Read the data point
                        read (200,err=1100) da1

                        ! Skip out of limits
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! Store
                        if (KSTK) Stokes0(ifreq,1,ith,iz) = da1

                        ! Integral
                        J00C(ifreq,iz) = J00C(ifreq,iz) + &
                                         da1*GeomI%W_mu(ith)

                      end do ! Frequency

                      ! Skip if out of limits
                      if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                      ! If storing
                      if (KSTK) then

                        ! For the rest of azimuths
                        do iph=2,GeomI%nph

                          ! Fill
                          Stokes0(:,iph,ith,iz) = Stokes0(:,1,ith,iz)

                        end do ! Rest of azimuths

                      end if ! Storing Stokes

                    ! Not axial input
                    else

                      ! Read for azimuthal directions
                      do iph=1,ia4

                        ! For each frequency
                        do ifreq=1,nfreq

                          ! Read the data point
                          read (200,err=1100) da1

                          ! Axial?
                          if (axiali.and.iph.gt.1) cycle

                          ! Skip out of limits
                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          ! Store
                          if (KSTK) Stokes0(ifreq,iph,ith,iz) = da1

                          ! Add contribution to the J00C integral
                          if (axiali) then

                            ! Integrate
                            J00C(ifreq,iz) = J00C(ifreq,iz) + &
                                             da1*GeomI%W_mu(ith)

                          ! Non axial
                          else

                            ! Integrate
                            J00C(ifreq,iz) = J00C(ifreq,iz) + da1* &
                                             GeomI%W_mu(ith)* &
                                             GeomI%W_mux(iph)
                          end if ! Axial

                        end do ! Frequencies
                      end do ! Input azimuths

                    end if ! Type of axial input

                  end do ! Polar
                end do ! Height

              !
              ! Doing AD or keeping Stokes
              !
              else

                ! Initialize
                J00C = 0d0

                ! For each height
                do iz=1,nZ

                  ! For each polar direction
                  do ith=1,GeomI%nTh

                    ! If input axial
                    if (inaxial) then

                      ! For each frequency
                      do ifreq=1,nfreq

                        ! Read the data
                        read (200,err=1100) da1

                        ! Skip out of limits
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! Save
                        Stokes0(ifreq,1,ith,iz) = da1

                      end do ! frequencies

                      ! Skip out of limits
                      if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                      ! Rest of azimuths
                      do iph=2,GeomI%nph

                        ! Fill
                        Stokes0(:,iph,ith,iz) = Stokes0(:,1,ith,iz)

                      end do ! azimuthal nodes

                    ! Not input axial
                    else

                      ! For each azimuthal direction
                      do iph=1,ia4

                        ! For each frequency
                        do ifreq=1,nfreq

                          ! Read the data
                          read (200,err=1100) da1

                          ! Skip out of limits
                          if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                          ! Axial?
                          if (axial.and.iph.gt.1) cycle

                          ! Save
                          Stokes0(ifreq,iph,ith,iz) = da1

                        end do ! frequencies
                      end do ! azimuthal nodes

                    end if ! Input axial

                  end do ! polar nodes

                  ! Skip out of limits
                  if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! For each polar direction
                    do ith=1,GeomI%nTh

                      ! For each azimuthal direction
                      do iph=1,GeomI%nPh

                        ! Add contribution to the JKQC integral
                        J00C(ifreq,iz) = J00C(ifreq,iz) + &
                                         Stokes0(ifreq,iph,ith,iz)* &
                                         GeomI%W_mu(ith)* &
                                         GeomI%W_mux(iph)

                      end do ! azimuthal nodes
                    end do ! polar nodes
                  end do ! frequencies
                end do ! heights

              end if ! Read AV or AD
            end if ! Master

            ! If MPI
            if (nproc.gt.1) then

              ! Control
              call control
              if (laborted) return

              ! Slaves require Stokes
              if (PRD.and.(dyn.or..not.AV)) then

                ! Buffer size
                ssize = nfreq*GeomI%nTh*GeomI%nPh*Rnz

                ! Share
                call MPI_BCAST(Stokes0(1,1,1,Rz0), ssize, &
                               MPI_DOUBLE_PRECISION, 0, &
                               MPI_COMM_RT, ierr)
              end if

              ! Buffer size
              csize = nfreq*RnZ

              ! Share
              call MPI_BCAST(J00C(1,Rz0), csize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)

            end if ! MPI
          end if ! Input AV or AD
        end if ! Can read Stokes
      end if ! Type of input

      ! Close unit
      if (pid.eq.0) close (200)

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file '//trim(filename)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error reading solution file '//trim(filename)
      close(100)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine readsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Restore the self-consistent solution from RAM\n
      !!     SolF(Solution_F_class): Structure with the solution of
      !!                             the self-consistent problem and
      !!                             the corresponding emergent
      !!                             profiles, contribution function,
      !!                             and height for optical depth
      !!                             equal to one\n
      !!      GeomI(Geometry_class): Structure with geometric data for
      !!                             the intensity problem\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!   Stokes0(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates\n
      !!         intensity(logical): Need to fetch intensity
      !!                             solution
      subroutine getsol(SolF,GeomI,Geom,Flgsg,Bfield,Atom, &
                        Stokes,JKQ,JKQS,JKQC, &
                        Stokes0,J00,J00S,J00C,J00P, &
                        intensity)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Geometry_class), intent(in):: Geom,GeomI
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Solution_F_class), intent(in):: SolF
      logical, intent(in):: intensity
      double precision, dimension(:,:,:,:), &
                        allocatable, intent(inout):: Stokes0
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(inout):: Stokes
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(inout):: J00P
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(inout):: JKQC

      ! Local

      integer:: ia,iz,i,itran,it,iJ,iJ1,K,iQ,iR,iR0,iR1,psize

      double precision:: rJ,rJ1,rho0

      complex(kind=8), dimension(-nkx:nkx,nz):: rhoKQaux


      ! Routine name
      urou = 'getsol'

      !
      ! If there is a polarization solution
      !
      if (intensity) then

        !
        ! Allocations
        !

        ! If keeping Stokes
        if (KSTK) then

          ! Allocate full height range
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1

        ! Not keeping Stokes
        else

          ! Allocate only two heights
          allocate(Stokes0(nfreq,GeomI%nPh,GeomI%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1

        end if

        ! Count memory
        RRAMc = RRAMc + 1d-6*sizeof(Stokes0)

        ! If there are atoms
        if (nA.gt.0) then

          ! J00 for absorptivity
          allocate(J00(nxt,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00)

          ! J00 for stimulated emission
          allocate(J00S(nxt,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00S)

          ! J00 for photoionizations
          allocate(J00P(nxphot,2,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00P)

        end if ! Atoms

        ! J00 frequency dependent
        allocate(J00C(nfreq,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(J00C)

        ! Master
        if (pid.eq.0) then

          ! If keeping Stokes
          if (KSTK) then

            ! Save
            Stokes0 = dble(SolF%i_StkI_b(:,:,:,Rz0:Rz1))

            ! If MPI
            if (nproc.gt.1) then

              ! Buffer size
              psize = nfreq*GeomI%nph*GeomI%nth*Rnz

              ! Share
              call MPI_BCAST(Stokes0(1,1,1,Rz0), psize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)
            end if ! MPI

          ! Otherwise
          else

            ! Initialize
            Stokes0 = 0d0

          end if

          ! There are atoms
          if (nA.gt.0) then

            ! J00 bar
            J00 = SolF%i_J00_b(:,Rz0:Rz1)

            ! MPI
            if (nproc.gt.1) then

              ! Buffer size
              psize = nxt*RnZ

              ! Share
              call MPI_BCAST(J00(1,Rz0), psize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)
            end if ! MPI

            ! J00 photo
            J00P = SolF%i_J00P_b(:,:,Rz0:Rz1)

            ! MPI
            if (nproc.gt.1) then

              ! Buffer size
              psize = nxphot*2*Rnz

              ! Share
              call MPI_BCAST(J00P(1,1,Rz0), psize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)
            end if ! MPI
          end if ! Atoms

          ! J00 freq. dependent
          J00C = SolF%i_J00C_b(:,Rz0:Rz1)

          ! MPI
          if (nproc.gt.1) then

            ! Buffer size
            psize = nfreq*Rnz

            ! Share
            call MPI_BCAST(J00C(1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
          end if ! MPI

        ! Slaves
        else

          ! If keeping Stokes
          if (KSTK) then

            ! Buffer size
            psize = nfreq*GeomI%nph*GeomI%nth*Rnz

            ! Get from Master
            call MPI_BCAST(Stokes0(1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Otherwise
          else

            ! Initialize
            Stokes0 = 0d0

          end if

          ! Atoms
          if (nA.gt.0) then

            ! J00
            psize = nxt*RnZ
            call MPI_BCAST(J00(1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

            ! J00P
            psize = nxphot*2*Rnz
            call MPI_BCAST(J00P(1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
          end if ! Atoms

          ! J00C
          psize = nfreq*Rnz
          call MPI_BCAST(J00C(1,Rz0), psize, &
                         MPI_DOUBLE_PRECISION, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Master/Slave

        ! Copy in J00S
        if (nA.gt.0) J00S = J00

        ! For each atom
        do ia=1,nA

          ! Initialize
          Atom(ia)%crho = cZero

          ! Master
          if (pid.eq.0) then

            ! Copy popu
            Atom(ia)%popu = SolF%i_rhoes_b(ia)%rho

          end if ! Master/slave

          !
          ! And share populations
          !

          ! MPI
          if (nproc.gt.1) then

            ! Buffer size
            psize = Rnz*Atom(ia)%nlevel

            ! Share
            call MPI_BCAST(Atom(ia)%popu(1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)
          end if ! MPI

          ! For each height
          do iz=Rz0,Rz1

            ! For each term
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Level and KQ
                i = Atom(ia)%irho(it)%irho_ij(iJ)
                iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)
                rJ = sqrt(2d0*rJ + 1d0)

                ! Set up rho
                Atom(ia)%crho(iR,iz) = &
                                     dcmplx(Atom(ia)%popu(i,iz), 0d0)

                ! Define population
                Atom(ia)%popu(i,iz) = Atom(ia)%popu(i,iz)*rJ

              end do ! Level
            end do ! Term
          end do ! Height


          !
          ! Flag the null rho(J,J')KQ
          !

          ! Initialize to non-null
          Atom(ia)%rhonull = .True.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)
              rJ = sqrt(2d0*rJ + 1d0)

              ! Get level index
              i = Atom(ia)%irho(it)%irho_ij(iJ)

              ! Get component index
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! For each height
              do iz=Rz0,Rz1

                ! If rhoKQ/rho00 is lesser than double precision
                ! flag null
                if (abs(Atom(ia)%crho(iR,iz)).gt.0d0) &
                  Atom(ia)%rhonull(iR,iz) = .False.

              end do ! heights
            end do ! J
          end do ! terms
        end do ! Atoms

      ! Polarization
      else

        ! If storing Stokes parameters
        if (KSTK) then

          ! Allocate whole height range
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz1))
          giz0 = Rz0
          giz1 = Rz1

        ! Not storing Stokes
        else

          ! Allocate just two nodes
          allocate(Stokes(0:3,nfreq,Geom%nPh,Geom%nTh,Rz0:Rz0+1))
          giz0 = Rz0
          giz1 = Rz0+1

        end if ! Storing Stokes

        ! Memory count
        RRAMc = RRAMc + 1d-6*sizeof(Stokes)

        ! If atoms
        if (nA.gt.0) then

          ! JKQ for absorptivity
          allocate(JKQ(-2:2,0:2,nxtran,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(JKQ)

          ! JKQ for stimulated emission
          allocate(JKQS(-2:2,0:2,nxtran,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(JKQS)

          ! J00 for photoionizations
          allocate(J00P(nxphot,2,Rz0:Rz1))
          RRAMc = RRAMc + 1d-6*sizeof(J00P)

        end if ! Atoms

        ! JKQ frequency dependent
        allocate(JKQC(-2:2,0:2,nfreq,Rz0:Rz1))
        RRAMc = RRAMc + 1d-6*sizeof(JKQC)

        ! Master
        if (pid.eq.0) then

          ! If keeping Stokes
          if (KSTK) then

            ! Save
            Stokes = SolF%i_Stk_b(:,:,:,:,Rz0:Rz1)

            ! MPI
            if (nproc.gt.1) then

              ! Buffer size
              psize = 4*nfreq*Geom%nph*Geom%nth*Rnz

              ! Share
              call MPI_BCAST(Stokes(0,1,1,1,Rz0), psize, &
                             MPI_DOUBLE_PRECISION, 0, &
                             MPI_COMM_RT, ierr)
            end if ! MPI

          ! Otherwise
          else

            ! Initialize
            Stokes = 0d0

          end if

          ! Atoms
          if (nA.gt.0) then

            ! JKQ bar
            JKQ = SolF%i_JKQ_b(:,:,:,Rz0:Rz1)

            !
            ! Rotate JKQ
            !

            ! For each height
            do iz=Rz0,Rz1

              ! No field, skip
              if (Bfield%Bstrength(iz).le.TINYB) cycle

              ! For each transition
              do itran=1,nxtran

                ! Rotate
                call fieldB(JKQ(:,:,itran,iz),1,Flgsg, &
                            Bfield%Btheta(iz),Bfield%Bphi(iz),1)

              end do ! Transitions
            end do ! heights

            ! MPI
            if (nproc.gt.1) then

              ! Buffer size
              psize = 15*nxtran*RnZ

              ! Share
              call MPI_BCAST(JKQ(-2,0,1,Rz0), psize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)
            end if ! MPI

            ! If stimulated emission
            if (stm) then

              ! JKQS bar
              JKQS = SolF%i_JKQS_b(:,:,:,Rz0:Rz1)

              !
              ! Rotate JKQ
              !

              ! For each height
              do iz=Rz0,Rz1

                ! No field, skip
                if (Bfield%Bstrength(iz).le.TINYB) cycle

                ! For each transition
                do itran=1,nxtran

                  ! Rotate
                  call fieldB(JKQS(:,:,itran,iz),1,Flgsg, &
                              Bfield%Btheta(iz),Bfield%Bphi(iz),1)

                end do ! Transitions
              end do ! heights

              ! MPI
              if (nproc.gt.1) then

                ! Buffer size
                psize = 15*nxtran*RnZ

                ! Share
                call MPI_BCAST(JKQS(-2,0,1,Rz0), psize, &
                               MPI_DOUBLE_COMPLEX, 0, &
                               MPI_COMM_RT, ierr)
              end if ! MPI
            end if ! Stimulated emission
          end if ! Atoms

          ! JKQ freq. dependent
          JKQC = SolF%i_JKQC_b(:,:,:,Rz0:Rz1)

          ! MPI
          if (nproc.gt.1) then

            ! Buffer size
            psize = 15*nfreq*Rnz

            ! Share
            call MPI_BCAST(JKQC(-2,0,1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)
          end if ! MPI

        ! Slaves
        else

          ! If keeping Stokes
          if (KSTK) then

            ! Buffer size
            psize = 4*nfreq*Geom%nph*Geom%nth*Rnz

            ! Receive
            call MPI_BCAST(Stokes(0,1,1,1,Rz0), psize, &
                           MPI_DOUBLE_PRECISION, 0, &
                           MPI_COMM_RT, ierr)

          ! Otherwise
          else

            ! Initialize
            Stokes = 0d0

          end if

          ! Atoms
          if (nA.gt.0) then

            ! JKQ
            psize = 15*nxtran*RnZ
            call MPI_BCAST(JKQ(-2,0,1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)

            ! If stimulatted
            if (stm) then

              ! JKQS
              psize = 15*nxtran*RnZ
              call MPI_BCAST(JKQS(-2,0,1,Rz0), psize, &
                             MPI_DOUBLE_COMPLEX, 0, &
                             MPI_COMM_RT, ierr)

            end if ! Stimulated emission
          end if ! Atoms

          ! JKQC
          psize = 15*nfreq*Rnz
          call MPI_BCAST(JKQC(-2,0,1,Rz0), psize, &
                         MPI_DOUBLE_COMPLEX, 0, &
                         MPI_COMM_RT, ierr)

        end if ! Master/slave

        ! For each atom
        do ia=1,nA

          ! Initialize
          Atom(ia)%crho = cZero

          ! Master
          if (pid.eq.0) then

            ! For each term
            do it=1,Atom(ia)%nMulti

              ! For each level
              do iJ=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J
                rJ = Atom(ia)%rJval(iJ,it)

                ! For each level
                do iJ1=1,Atom(ia)%nJ(it)

                  ! Get J'
                  rJ1 = Atom(ia)%rJval(iJ1,it)

                  ! For each K
                  do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                    ! For each Q
                    do iQ=-K,K

                      ! Get rho(J,J')KQ index
                      if (K.le.Atom(ia)%Kcut(it)) &
                        iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! For each height
                      do iz=1,nZ

                        ! Skip is not in range
                        if (iz.lt.Rz0.or.iz.gt.Rz1) cycle

                        ! If considered multipole
                        if (K.le.Atom(ia)%Kcut(it)) then

                          ! Save
                          Atom(ia)%crho(iR,iz) = &
                                        SolF%i_rhoes_b(ia)%crho(iR,iz)

                          ! Auxiliar variable for rotation
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end if ! Considered multipole

                      end do ! heights
                    end do ! Q

                    !
                    ! Rotate rhoKQ
                    !

                    ! Only if not above the limit
                    if (K.le.Atom(ia)%Kcut(it)) then

                      ! For each height in the CPU domain
                      do iz=Rz0,Rz1

                        ! Skip if no magnetic field
                        if (Bfield%Bstrength(iz).le.TINYB) cycle

                        ! Rotate the rhoKQ in the auxiliar variable
                        call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                  Bfield%Btheta(iz), &
                                  Bfield%Bphi(iz),1)

                        ! For every Q
                        do iQ=-K,K

                          ! Get index
                          iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                          ! Store the rotated result in the rhoKQ
                          ! array
                          Atom(ia)%crho(iR,iz) = rhoKQaux(iQ,iz)

                        end do ! Q
                      end do ! heights

                    end if ! K < Kcut

                  end do ! K
                end do ! J'
              end do ! J
            end do ! terms

          end if ! Master/slave

          ! MPI
          if (nproc.gt.1) then

            ! Buffer size
            psize = Atom(ia)%ndim*Rnz

            ! Share
            call MPI_BCAST(Atom(ia)%crho(1,Rz0), psize, &
                           MPI_DOUBLE_COMPLEX, 0, &
                           MPI_COMM_RT, ierr)
          end if ! MPI


          !
          ! Flag the null rho(J,J')KQ
          !

          ! Initialize to non-null
          Atom(ia)%rhonull = .False.

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level
            do iJ=1,Atom(ia)%nJ(it)

              ! Get J
              rJ = Atom(ia)%rJval(iJ,it)

              ! Get the rho00 indexes
              iR0 = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! For each level
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J'
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! Get the rho00 indexes
                iR1 = Atom(ia)%irho(it)%Jrho(iJ1,iJ1)%kq(0,0)

                ! For each K
                do K=nint(abs(rJ-rJ1)), &
                     min(nint(rJ+rJ1),Atom(ia)%Kcut(it))

                  ! For each Q
                  do iQ=-K,K

                    ! Get rho(J,J')KQ index
                    iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                    ! For each height in this domain
                    do iz=Rz0,Rz1

                      ! Get the inverse of rho00
                      rho0 = 1d0/sqrt(abs(Atom(ia)%crho(iR0,iz))* &
                                      abs(Atom(ia)%crho(iR1,iz)))

                      ! If rhoKQ/rho00 is lesser than double precision
                      ! flag null
                      if (abs(Atom(ia)%crho(iR,iz)*rho0).lt.TINYR) &
                        Atom(ia)%rhonull(iR,iz) = .True.

                    end do ! heights
                  end do ! Q
                end do ! K
              end do ! J'
            end do ! J
          end do ! terms
        end do ! atoms

      end if ! Intensity/polarization initialization

      end subroutine getsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the self-consistent solution in RAM\n
      !!     SolF(Solution_F_class): Structure with the solution of
      !!                             the self-consistent problem and
      !!                             the corresponding emergent
      !!                             profiles, contribution function,
      !!                             and height for optical depth
      !!                             equal to one\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence\n
      !!   Stokes0(double(:,:,:,:)): Intensity\n
      !!           J00(double(:,:)): Mean intensity integrated over
      !!                             the absorption profile\n
      !!          J00S(double(:,:)): Mean intensity integrated over
      !!                             the emission profile\n
      !!          J00C(double(:,:)): Mean intensity with frequency
      !!                             dependence\n
      !!        J00P(double(:,:,:)): Intensity integrals in the
      !!                             photoionization rates\n
      !!         intensity(logical): Need to store just the intensity
      !!                             solution
      subroutine setsol(SolF,Flgsg,Bfield,Atom, &
                        Stokes,JKQ,JKQS,JKQC, &
                        Stokes0,J00,J00S,J00C,J00P, &
                        intensity)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      type(Solution_F_class), intent(inout):: SolF
      logical, intent(in):: intensity
      double precision, dimension(:,:,:,:),  &
                        allocatable, intent(in):: Stokes0
      double precision, dimension(:,:,:,:,:), &
                        allocatable, intent(in):: Stokes
      double precision, dimension(:,:), &
                        allocatable, intent(in):: J00
      double precision, dimension(:,:), &
                        allocatable, intent(in):: J00S
      double precision, dimension(:,:), &
                        allocatable, intent(in):: J00C
      double precision, dimension(:,:,:), &
                        allocatable, intent(in):: J00P
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(in):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(in):: JKQS
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(in):: JKQC

      ! Local

      integer:: ia,it,iJ,iJ1,iR,ii,itran,iz,K,iQ

      double precision:: rJ,rJ1

      complex(kind=8), dimension(-nkx:nkx,nz):: rhoKQaux


      ! Only master, rest can leave
      if (pid.gt.0) return


      !
      ! Intensity
      !
      if (intensity) then

        !
        ! Populations
        !

        ! For each atom
        do ia=1,nA

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J within the term
            do iJ=1,Atom(ia)%nJ(it)!,1,-1

              ! Get the level and rho00 index
              ii = Atom(ia)%irho(it)%irho_ij(iJ)
              iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

              ! Save rho00
              SolF%i_rhoes(ia)%rho(ii,:) = dble(Atom(ia)%crho(iR,:))

            end do ! Sublevel
          end do ! Term
        end do ! Atom

        ! Atoms
        if (nA.gt.0) then

          !
          ! Mean radiation field tensors
          !
          SolF%i_J00(:,1:Rz0-1) = 0d0
          SolF%i_J00(:,Rz1+1:nz) = 0d0
          SolF%i_J00(:,Rz0:Rz1) = J00(:,Rz0:Rz1)

          !
          ! Photoionization
          !
          SolF%i_J00P(:,:,1:Rz0-1) = 0d0
          SolF%i_J00P(:,:,Rz1+1:nz) = 0d0
          SolF%i_J00P(:,:,Rz0:Rz1) = J00P(:,:,Rz0:Rz1)

        end if ! Atoms

        !
        ! Radiation field tensors
        !
        SolF%i_J00C(:,1:Rz0-1) = 0d0
        SolF%i_J00C(:,Rz1+1:nz) = 0d0
        SolF%i_J00C(:,Rz0:Rz1) = J00C(:,Rz0:Rz1)

        !
        ! Stokes
        !
        if (KSTK) then
          SolF%i_StkI(:,:,:,1:Rz0-1) = 0d0
          SolF%i_StkI(:,:,:,Rz1+1:nz) = 0d0
          SolF%i_StkI(:,:,:,Rz0:Rz1) = Stokes0(:,:,:,Rz0:Rz1)
        end if

      !
      ! Polarization
      !
      else

        !
        ! Density matrix
        !

        ! For each atom
        do ia=1,nA

          ! Initialize
          SolF%i_rhoes(ia)%crho = cZero

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! For each level J within the term
            do iJ=1,Atom(ia)%nJ(it)!,1,-1

              ! Get J value
              rJ = Atom(ia)%rJval(iJ,it)

              ! For each J value within the term
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J value
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! For each K
                do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                  !
                  ! Rotate rhoKQ nto the vertical reference frame
                  if (K.le.Atom(ia)%Kcut(it)) then

                    ! For each height
                    do iz=1,nz

                      ! If out of bounds
                      if (iz.lt.Rz0.or.iz.gt.Rz1) then

                        ! Set to zero
                        rhoKQaux(-K:K,iz) = cZero

                      ! In bounds
                      else

                        ! For each Q
                        do iQ=-K,K

                          ! Get the rhoKQ index
                          iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                          ! Store the corresponding rhoKQ into an
                          ! auxiliar variable
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end do ! Q

                        ! If there is non-zero magnetic field, rotate
                        if (Bfield%Bstrength(iz).gt.TINYB) &
                          call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                   -Bfield%Btheta(iz), &
                                   -Bfield%Bphi(iz),-1)

                      end if ! Height bounds

                    end do ! heights

                    ! For each Q
                    do iQ=-K,K

                      ! Get the index
                      iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! Save
                      SolF%i_rhoes(ia)%crho(iR,:) = &
                                                   Atom(ia)%crho(iR,:)

                    end do ! Q

                  end if ! K<=Kcut

                end do ! K
              end do ! J'
            end do ! J
          end do ! Terms
        end do ! Atoms

        ! Atoms
        if (nA.gt.0) then

          !
          ! Mean radiation field tensors
          !

          ! Heights
          do iz=1,nz

            ! Out of limits
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! Set to zero and skip
              SolF%i_JKQ(:,:,:,iz) = cZero
              if (stm) SolF%i_JKQS(:,:,:,iz) = cZero
              cycle

            end if ! Out of limits

            ! Store JKQ and JKQS
            SolF%i_JKQ(:,:,:,iz) = JKQ(:,:,:,iz)
            if (stm) SolF%i_JKQS(:,:,:,iz) = JKQS(:,:,:,iz)

            ! If there is a magnetic field
            if (Bfield%Bstrength(iz).gt.TINYB) then

              ! For each transition
              do itran=1,nxtran

                ! Rotate
                call fieldB(SolF%i_JKQ(:,:,itran,iz),1,Flgsg, &
                            -Bfield%Btheta(iz),-Bfield%Bphi(iz),-1)

                ! If stimulated emission, rotate
                if (stm) &
                call fieldB(SolF%i_JKQS(:,:,itran,iz),1,Flgsg, &
                            -Bfield%Btheta(iz),-Bfield%Bphi(iz),-1)

              end do ! Transitions

            end if ! Magnetic field

          end do ! Heights

        end if ! Atoms

        !
        ! JKQC radiation field tensors
        !
        SolF%i_JKQC(:,:,:,1:Rz0-1) = cZero
        SolF%i_JKQC(:,:,:,Rz1+1:nz) = cZero
        SolF%i_JKQC(:,:,:,Rz0:Rz1) = JKQC(:,:,:,Rz0:Rz1)

        !
        ! Stokes
        !
        if (KSTK) then
          SolF%i_Stk(:,:,:,:,1:Rz0-1) = 0d0
          SolF%i_Stk(:,:,:,:,Rz1+1:nz) = 0d0
          SolF%i_Stk(:,:,:,:,Rz0:Rz1) = Stokes(:,:,:,:,Rz0:Rz1)
        end if

      end if ! Polarization/intensity

      return

      ! Fake assing to deceive the compiler
      rJ = J00S(1,1)

      end subroutine setsol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the self-consistent solution in a file\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!         suff(character(:)): Suffix for the file names of the
      !!                             radiation field and density
      !!                             matrix files\n
      !!           omega(double(:)): Frequency array\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!               z(double(:)): Height array\n
      !!  Stokes(double(:,:,:,:,:)): Stokes parameters\n
      !!     JKQ(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the absorption
      !!                             profile\n
      !!    JKQS(dcomplex(:,:,:,:)): Radiation field tensors
      !!                             integrated over the emission
      !!                             profile\n
      !!    JKQC(dcomplex(:,:,:,:)): Radiation field tensors with
      !!                             frequency dependence
      subroutine writesol(Input,suff,omega,Geom, &
                          Flgsg,Bfield,Atom,z,Stokes,JKQ,JKQS,JKQC)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(inout):: Atom
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(in):: Geom
      type(Fctsg_class), intent(in):: Flgsg
      type(Bfield_class), intent(in):: Bfield
      character(len=4), intent(in):: suff
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(0:3,nfreq,Geom%nPh,Geom%nTh, &
                                     giz0:giz1), intent(in):: Stokes
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(in):: JKQ
      complex(kind=8), dimension(:,:,:,:), &
                       allocatable, intent(in):: JKQS
      complex(kind=8), &
             dimension(-2:2,0:2,nfreq,Rz0:Rz1), intent(in):: JKQC

      ! Local

      integer, parameter:: izero = 0
      integer, parameter:: ione = 1
      double precision, parameter:: dzero = 0d0

      character(len=8):: scoord
      character(len=500):: filename

      logical:: saveJKQnu,saverKQ,saveJKQ,saveP,saveD,saveS,saveSol
      logical:: laux

      integer:: ierr,ii,iab,iran,ios,ia,iz,ifreq,i,ith,iph,iS,itran
      integer:: jtran,it,iJ,iJ1,K,iQ,iR,axial_int,stm_int,AV_int
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: rJ,rJ1,loffset

      complex(kind=8), dimension(-2:2,0:2,nz):: JKQaux
      complex(kind=8), dimension(-nkx:nkx,nz):: rhoKQaux


      ! Slaves
      if (pid.gt.0) then

        ! Control
        call control
        return

      end if ! Slaves

      !
      ! Shorter name and compute bools
      !
      filename = Input%folder
      saveSol = Input%keep_sol
      saveP = Input%keep_pop.and.nA.gt.0.and. &
              (suff.eq.'NONE'.or.run_mode.eq.0)
      saveD = Input%keep_dep.and.nA.gt.0.and. &
              (suff.eq.'NONE'.or.run_mode.eq.0)
      saverKQ = Input%keep_rhoKQ.and.nA.gt.0.and. &
                (suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQ = Input%keep_JKQ.and.nA.gt.0.and. &
                (suff.eq.'NONE'.or.run_mode.eq.0)
      saveS = Input%keep_stokesQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQnu = Input%keep_jkqnu.and.(suff.eq.'NONE'.or. &
                                        run_mode.eq.0)

      ! If not writing anything
      if (.not.(saveSol.or.saveP.or.saveD.or.saverKQ.or.saveJKQ.or. &
                saveS.or.saveJKQnu)) then

        ! Leave
        call control
        return

      end if ! Not writing

      ! Routine name
      urou = 'writesol'

      !
      ! Open files
      !

      ! Get LOS index string if 1.5D synthesis
      if (run_mode.eq.1) write(scoord,'(I0.8)') icoords(3)

      ! If writing solution
      if (saveSol) then

        ! If 1D
        if (run_mode.eq.0) then

          ! Open solution file
          open (200,file=trim(filename)//'/Solution', &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', &
                form='unformatted')

        ! If 1.5D
        else if (run_mode.eq.1) then

          ! Open solution file
          open (200,file=trim(filename)// &
                '/Solution-folder/Solution-'// &
                scoord, status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', form='unformatted')

        end if ! 1D vs 1.5D

      end if ! Writing solution

      ! If saving JKQ
      if (saveJKQ) then

        ! To write the final JKQ
        if (suff.eq.'NONE') then

          ! If 1D
          if (run_mode.eq.0) then

            ! Open file
            open (300,file=trim(filename)//'/Jout', &
                  status='unknown', iostat=ios, err=1002, &
                  access='stream', action='write', &
                  form='unformatted')

          ! 1.5D
          else if (run_mode.eq.1) then

            ! Open file
            open (300,file=trim(filename)// &
                  '/Solution-folder/Jout-'//scoord,&
                  status='unknown', iostat=ios, err=1002, &
                  access='stream', action='write', &
                  form='unformatted')

          end if ! 1D vs 1.5D

        ! With suffix
        else

          ! Open for 1D
          open (300,file=trim(filename)//'/Jout_'//suff, &
                status='unknown', iostat=ios, err=1002, &
                access='stream', action='write', &
                form='unformatted')

        end if ! Final or intermediate
      end if ! Saving JKQ

      ! If saving rhoKQ
      if (saverKQ) then

        ! To write the final rhoKQ
        if (suff.eq.'NONE') then

          ! If 1D
          if (run_mode.eq.0) then

            ! Open file
            open (400,file=trim(filename)//'/Rhoout', &
                  status='unknown', iostat=ios, err=1003, &
                  access='stream', action='write', &
                  form='unformatted')

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! Open file
            open (400,file=trim(filename)// &
                  '/Solution-folder/Rhoout-'// &
                  scoord, status='unknown', iostat=ios, err=1003, &
                  access='stream', action='write', &
                  form='unformatted')

          end if ! 1D vs 1.5D

        ! Intermediate file
        else

          ! Open for 1D
          open (400,file=trim(filename)//'/Rhoout_'//suff, &
                status='unknown', iostat=ios, err=1003, &
                access='stream', action='write', &
                form='unformatted')

        end if ! Final or intermediate file
      end if ! Saving rhoKQ

      !
      ! Convert logicals to integers
      !

      ! Axial symmetry
      if(Geom%axial)then
        axial_int = 1
      else
        axial_int = 0
      end if

      ! Stimulated emission
      if(stm)then
        stm_int = 1
      else
        stm_int = 0
      end if

      ! Angle averaged redistribution function
      if (KSTK) then
        AV_int = 0
      else
        AV_int = 1
      end if

      !
      ! Write headers with dimensions and flags
      !

      ! If saving solution
      if (saveSol) then

        ! Solution file metadata
        write(200,err=1100) 'sp'
        write(200,err=1100) nfreq
        write(200,err=1100) nZ
        write(200,err=1100) Geom%nTh
        write(200,err=1100) Geom%nPh
        write(200,err=1100) nA
        write(200,err=1100) axial_int
        write(200,err=1100) stm_int
        write(200,err=1100) AV_int

      end if ! Saving solution

      ! If saving JKQ
      if (saveJKQ) then

        ! JKQ file metadata
        write(300,err=1102) 'bj'
        write(300,err=1102) stm
        write(300,err=1102) nZ
        write(300,err=1102) nA
        write(300,err=1102) nxtran
        write(300,err=1102) z

      end if ! Saving JKQ

      ! If saving rhoKQ
      if (saverKQ) then

        ! rhoKQ file metadata
        write(400,err=1103) 'br'
        write(400,err=1103) nZ
        write(400,err=1103) nA
        write(400,err=1103) z

      end if ! Saving rhoKQ

      ! If 1.5D
      if (run_mode.eq.1) then

        ! If saving population or departure, prepare buffers
        if (saveP.or.saveD) &
          allocate(buffer(maxval(Input%lim_pop%nbuff)/4))

      end if ! 1.5D synthesis


      !
      ! Write the data
      !

      ! Only if saving anything
      if (saveSol.or.saveP.or.saveD.or.saverKQ) then

        !
        ! Population and rhoKQ

        ! For each atom
        do ia=1,nA

          !
          ! If saving population or departure coeff.
          !
          if (saveP.or.saveD) then

            !
            ! If 1D
            !
            if (run_mode.eq.0) then

              ! If writing populations
              if (saveP) then

                ! To write the populations
                open (500,file=trim(filename)//'/'// &
                      trim(Atom(ia)%file_label)//'.pop', &
                      status='unknown',iostat=ios, err=1004, &
                      access='stream', action='write', &
                      form='unformatted')

                ! Metadata
                write(500,err=1104) 'bp'
                write(500,err=1104) nZ
                write(500,err=1104) Atom(ia)%nlevel

              end if ! Saving population

              ! If writing departure c.
              if (saveD) then

                ! To write departure c.
                open (600,file=trim(filename)//'/'// &
                      trim(Atom(ia)%file_label)//'.dep', &
                      status='unknown',iostat=ios, err=1005, &
                      access='stream', action='write', &
                      form='unformatted')

                ! Metadata
                write(600,err=1105) 'bb'
                write(600,err=1105) nZ
                write(600,err=1105) Atom(ia)%nlevel

              end if ! Saving depature coeff.

              ! For each height
              do iz=1,nZ

                ! If out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! For each level
                  do it=1,Atom(ia)%nlevel

                    ! Write populations if requested
                    if (saveP) &
                      write(500,err=1104) Atom(ia)%popu(it,iz)

                    ! Write departure if requested
                    if (saveD) &
                      write(600,err=1104) Atom(ia)%popu(it,iz)/ &
                                          Atom(ia)%populte(it,iz)

                  end do ! Levels

                ! In bounds
                else

                  ! For each term
                  do it=1,Atom(ia)%nMulti

                    ! For each level J within the term
                    do iJ=1,Atom(ia)%nJ(it)!,1,-1

                      ! Get J value and level index
                      rJ = Atom(ia)%rJval(iJ,it)
                      i = Atom(ia)%irho(it)%irho_ij(iJ)
                      iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                      ! Write populations if requested
                      if (saveP) &
                        write(500,err=1104) Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))

                      ! Write departure coeff if requested
                      if (saveD) &
                        write(600,err=1105) Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))/ &
                                         Atom(ia)%populte(i,iz)
                    end do ! Levels
                  end do ! Terms

                end if ! Height bounds

              end do ! Heights

              ! Close files if opened
              if (saveP) close(500)
              if (saveD) close(600)

            !
            ! If 1.5D
            !
            else if (run_mode.eq.1) then

              ! Populations
              if (saveP.and.Input%lim_pop%nbuff(ia).gt.0) then

                ! Open file to write the populations
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.pop', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
                if (ierr.ne.0) goto 1004

                !
                ! Column offset
                !

                ! Get offset
                loffset = dble(icoords(3)-1)* &
                          dble(Input%lim_pop%nbuff(ia)) + &
                          dble(Input%lim_pop%head_size)
                do while(loffset.gt.offlimit)
                  offset = int(offlimit)
                  call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                  if (ierr.ne.0) goto 1014
                  loffset = loffset - offlimit
                end do
                offset = int(loffset)
                call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                if (ierr.ne.0) goto 1014

                ! Initialize buffer
                ii = 0

                ! If specified ranges
                if (Input%lim_pop%nran.gt.0) then

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! This atom not included, skip
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz))

                      end do ! Entries

                    ! In bounds
                    else

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz)))

                      end do ! Ranges to print

                    end if ! Height bounds

                  end do ! Heights

                ! Everything
                else

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each level
                      do it=1,Atom(ia)%nlevel

                        ! Advance index
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(it,iz))

                      end do ! Levels

                    ! In bounds
                    else

                      ! For each term
                      do it=1,Atom(ia)%nMulti

                        ! For each level J within the term
                        do iJ=1,Atom(ia)%nJ(it)!,1,-1

                          ! Get J value and level index
                          rJ = Atom(ia)%rJval(iJ,it)
                          iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                          ! Advance index
                          ii = ii + 1

                          ! Save
                          buffer(ii) = real(Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz)))
                        end do ! Levels
                      end do ! Terms

                    end if ! Height bounds

                  end do ! Heights

                end if ! Specific or everything

                ! Write buffer
                call MPI_FILE_WRITE(funit,buffer(1), &
                                    Input%lim_pop%nbuff(ia)/4, &
                                    MPI_REAL,MPI_STATUS_IGNORE,ierr)
                if (ierr.ne.0) goto 1304

              end if ! Saving populations

              ! Departure coefficients
              if (saveD.and.Input%lim_pop%nbuff(ia).gt.0) then

                ! Open file to write the populations
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.dep', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
                if (ierr.ne.0) goto 1005

                !
                ! Column offset
                !

                ! Get offset
                loffset = dble(icoords(3)-1)* &
                          dble(Input%lim_pop%nbuff(ia)) + &
                          dble(Input%lim_pop%head_size)
                do while(loffset.gt.offlimit)
                  offset = int(offlimit)
                  call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                  if (ierr.ne.0) goto 1015
                  loffset = loffset - offlimit
                end do
                offset = int(loffset)
                call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                if (ierr.ne.0) goto 1015

                ! Initialize buffer
                ii = 0

                ! If specified
                if (Input%lim_pop%nran.gt.0) then

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Advance buffer
                        ii = ii +1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                          Atom(ia)%populte(i,iz))

                      end do ! Ranges to print

                    ! In bounds
                    else

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance buffer
                        ii = ii +1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))

                      end do ! Ranges to print

                    end if ! Height bounds

                  end do ! Height

                ! Everything
                else

                  ! Each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each term
                      do it=1,Atom(ia)%nlevel

                        ! Advance index
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(it,iz)/ &
                                          Atom(ia)%populte(i,iz))

                      end do ! Levels

                    ! In bounds
                    else

                      ! For each term
                      do it=1,Atom(ia)%nMulti

                        ! For each level J within the term
                        do iJ=1,Atom(ia)%nJ(it)!,1,-1

                          ! Get J value and level index
                          rJ = Atom(ia)%rJval(iJ,it)
                          i = Atom(ia)%irho(it)%irho_ij(iJ)
                          iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                          ! Advance index
                          ii = ii + 1

                          ! Save
                          buffer(ii) = real(Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))/ &
                                         Atom(ia)%populte(i,iz))

                        end do ! Levels
                      end do ! Terms

                    end if ! Height bounds

                  end do ! Heights

                end if ! Specific or everything

                ! Write buffer
                call MPI_FILE_WRITE(funit,buffer(1), &
                                    Input%lim_pop%nbuff(ia)/4, &
                                    MPI_REAL,MPI_STATUS_IGNORE,ierr)
                if (ierr.ne.0) goto 1305

              end if! Saving populations

            end if! 1D vs 1.5D
          end if ! Saving populations or departure coefficients

          ! Write the population of the atom into solution
          if (saveSol) write(200,err=1100) Atom(ia)%n

          ! Saving rhoKQ
          if (saverKQ) then

            ! Write the population of the atom into rhoKQ files
            write(400,err=1103) Atom(ia)%n

            ! Write number of terms to rhoKQ file
            write(400,err=1103) Atom(ia)%nMulti

          end if

          ! For each term
          do it=1,Atom(ia)%nMulti

            ! Write number of levels to rhoKQ file
            if (saverKQ) write(400,err=1103) Atom(ia)%nJ(it)

            ! For each level J within the term
            do iJ=1,Atom(ia)%nJ(it)!,1,-1

              ! Get J value
              rJ = Atom(ia)%rJval(iJ,it)

              ! For each J value within the term
              do iJ1=1,Atom(ia)%nJ(it)!,1,-1

                ! Get J value
                rJ1 = Atom(ia)%rJval(iJ1,it)

                ! Write J values in rhoKQ file
                if (saverKQ) then
                  write(400,err=1103) nint(2d0*rJ)
                  write(400,err=1103) nint(2d0*rJ1)
                end if

                ! For each K
                do K=nint(abs(rJ-rJ1)),nint(rJ+rJ1)

                  ! If included multipole
                  if (K.le.Atom(ia)%Kcut(it)) then

                    ! For each height
                    do iz=1,nz

                      ! If out of bounds
                      if (iz.lt.Rz0.or.iz.gt.Rz1) then

                        ! Set to zero
                        rhoKQaux(-K:K,iz) = cZero

                      ! In bounds
                      else

                        ! For each Q
                        do iQ=-K,K

                          ! Get the rhoKQ index
                          iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                          ! Store the corresponding rhoKQ into an
                          ! auxiliar variable
                          rhoKQaux(iQ,iz) = Atom(ia)%crho(iR,iz)

                        end do ! Q

                        ! If there is non-zero magnetic field, rotate
                        if (Bfield%Bstrength(iz).gt.TINYB) &
                          call rhoB(rhoKQaux(-K:K,iz),1,K,Flgsg, &
                                   -Bfield%Btheta(iz), &
                                   -Bfield%Bphi(iz),-1)

                      end if ! Height bounds

                    end do ! heights

                    ! For each Q
                    do iQ=-K,K

                      ! Get the index
                      iR = Atom(ia)%irho(it)%Jrho(iJ1,iJ)%kq(iQ,K)

                      ! For each height
                      do iz=1,nZ

                        ! If saving rhoKQ
                        if (saverKQ) then

                          ! Write rhoKQ into rhoKQ file, and the null
                          write(400,err=1103) dble(rhoKQaux(iQ,iz))
                          write(400,err=1103) dimag(rhoKQaux(iQ,iz))

                          ! If out of bounds
                          if (iz.lt.Rz0.or.iz.gt.Rz1) then

                            ! Write a one for null
                            write(400,err=1103) ione

                          ! In bounds
                          else

                            ! If null
                            if (Atom(ia)%rhonull(iR,iz)) then

                              ! Write one
                              write(400,err=1103) ione

                            ! If not null
                            else

                              ! Write a zero
                              write(400,err=1103) izero

                            end if ! Null value
                          end if ! Height bounds
                        end if ! Saving rhoKQ

                        ! Write rhoKQ into solution file
                        if (saveSol) &
                          write(200,err=1100) &
                                           dble(rhoKQaux(iQ,iz)), &
                                           dimag(rhoKQaux(iQ,iz))

                      end do ! heights
                    end do ! Q

                  ! Not included multipole
                  else

                    ! For each Q
                    do iQ=-K,K

                      ! For each height
                      do iz=1,nZ

                        ! If saving rhoKQ
                        if (saverKQ) then

                          ! Write rhoKQ into rhoKQ file, and the null
                          ! flag
                          write(400,err=1103) dzero
                          write(400,err=1103) dzero
                          write(400,err=1103) ione

                        end if ! Saving rhoKQ

                        ! Write rhoKQ into solution file
                        if (saveSol) &
                          write(200,err=1100) dzero,dzero

                      end do ! heights
                    end do ! Q

                  end if ! K<=Kcut

                end do ! K
              end do ! J'
            end do ! J
          end do ! Terms
        end do ! Atoms

      end if ! Saving rhoKQ or the solution file


      !
      ! JKQ

      ! Only if saving anything in this block
      if (saveSol.or.saveJKQ) then

        ! For each atom
        do ia=1,nA

          ! Write number of transitions in JKQ file
          if (saveJKQ) write(300,err=1102) Atom(ia)%ntran

          ! For each transition
          do itran=1,Atom(ia)%ntran

            ! Apply atomic shift
            jtran = itran + Atom(ia)%tshift

            !
            ! Rotate JKQ

            ! For each height
            do iz=1,nz

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Set to zero
                JKQaux(:,:,iz) = cZero

              ! In bounds
              else

                ! get the jkq into an auxiliar variable
                JKQaux(:,:,iz) = JKQ(:,:,jtran,iz)

                ! If there is a non-zero magnetic field, rotate
                ! into the vertical reference frame
                if (Bfield%Bstrength(iz).gt.TINYB) &
                  call fieldB(JKQaux(:,:,iz),1,Flgsg, &
                              -Bfield%Btheta(iz), &
                              -Bfield%Bphi(iz),-1)

              end if ! Height bounds

            end do ! heights

            ! For each K
            do K=0,2

              ! For each Q
              do iQ=-K,K

                ! For each height
                do iz=1,nZ

                  ! Write the JKQ into the JKQ file
                  if (saveJKQ) &
                    write(300,err=1102) dble(JKQaux(iQ,K,iz)), &
                                        dimag(JKQaux(iQ,K,iz))

                  ! Write the JKQ into the solution file
                  if (saveSol) &
                    write(200,err=1100) dble(JKQaux(iQ,K,iz)), &
                                        dimag(JKQaux(iQ,K,iz))

                end do ! heights
              end do ! Q
            end do ! K
          end do ! transitions
        end do ! Atoms

        ! If there is stimulated emission, write JKQS
        if(stm)then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tshift

              !
              ! Rotate JKQ

              ! For each height
              do iz=1,nz

                ! If out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! Set to zero
                  JKQaux(:,:,iz) = cZero

                ! In bounds
                else

                  ! Get the Jkq into an auxiliar variable
                  JKQaux(:,:,iz) = JKQS(:,:,jtran,iz)

                  ! If there is a non-zero magnetic field, rotate
                  ! into the vertical reference frame
                  if (Bfield%Bstrength(iz).gt.TINYB) &
                    call fieldB(JKQaux(:,:,iz),1,Flgsg, &
                                -Bfield%Btheta(iz), &
                                -Bfield%Bphi(iz),-1)

                end if ! Height bounds

              end do ! heights

              ! For each K
              do K=0,2

                ! For each Q
                do iQ=-K,K

                  ! For each height
                  do iz=1,nZ

                    ! Write the JKQS into the JKQ file
                    if (saveJKQ) &
                      write(300,err=1102) dble(JKQaux(iQ,K,iz)), &
                                           dimag(JKQaux(iQ,K,iz))

                    ! Write the JKQS into the solution file
                    if (saveSol) &
                      write(200,err=1100) dble(JKQaux(iQ,K,iz)), &
                                           dimag(JKQaux(iQ,K,iz))

                  end do ! heights
                end do ! Q
              end do ! K
            end do ! transitions
          end do ! atoms

        end if ! stimulated emission
      end if ! Saving anything

      !
      ! Radiation field solution
      !

      ! Saving solution
      if (saveSol) then

        !
        ! Keeping full Stokes
        if (KSTK) then

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each degree of freedom
              do ith=1,Geom%nTh
                do iph=1,Geom%nPh
                  do ifreq=1,nfreq
                    do iS=0,3

                      ! Write 0
                      write(200,err=1100) dzero

                    end do
                  end do
                end do
              end do

            ! In bounds
            else

              ! For each polar direction
              do ith=1,Geom%nTh

                ! For each azimuthal direction
                do iph=1,Geom%nPh

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! For each Stokes parameter
                    do iS=0,3

                      ! Write Stokes parameter into the solution file
                      write(200,err=1100) Stokes(iS,ifreq,iph,ith,iz)

                    end do ! Stokes parameters
                  end do ! Frequencies
                end do ! azimuthal directions
              end do ! polar directions

            end if ! Height bounds

          end do ! heights

        ! Not keeping Stokes
        else

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each degree of freedom
              do ifreq=1,nfreq
                do K=0,2
                  do iQ=-K,K

                    ! Write 0 0
                    write(200,err=1100) dzero,dzero

                  end do
                end do
              end do

            ! In bounds
            else

              ! For each frequency
              do ifreq=1,nfreq

                ! For each K
                do K=0,2

                  ! For each Q
                  do iQ=-K,K

                    ! Write the JKQ(k) into the solution file
                    write(200,err=1100) dble(JKQC(iQ,K,ifreq,iz)), &
                                         dimag(JKQC(iQ,K,ifreq,iz))

                  end do ! Q
                end do ! K
              end do ! frequencies

            end if ! Height bounds

          end do ! heights

        end if ! AV or AD

        !
        ! Close files
        !
        close (200)

      end if ! Saving solution

      ! Close files if opened
      if (saveJKQ) close (300)
      if (saverKQ) close (400)

      ! If not storing anything else
      if (.not.saveJKQnu.and..not.saveS) then

        ! Free
        if (allocated(buffer)) deallocate(buffer)

        ! Control
        call control
        return

      end if ! Not storing anything else


      !
      ! Store the stokes in the quadrature
      !

      ! If there is no suffix
      if (suff.eq.'NONE') then

        ! Saving Stokes
        if (saveS) then

          ! If 1D
          if (run_mode.eq.0) then

            ! Open file
            open (250,file=trim(filename)//'/Stokesout', &
                  status='unknown', iostat=ios, err=1001, &
                  access='stream', action='write', &
                  form='unformatted')

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! Open file
            open (250,file=trim(filename)// &
                  '/Solution-folder/Stokesout-'//scoord, &
                  status='unknown', iostat=ios, err=1001, &
                  access='stream', action='write', &
                  form='unformatted')

          end if ! 1D vs 1.5D
        end if ! Saving Stokes

        ! JKQnu file
        if (saveJKQnu) then

          ! If 1D
          if (run_mode.eq.0) then

            ! Open file
            open (350,file=trim(filename)//'/JKQnuout', &
                  status='unknown', iostat=ios, err=1006, &
                  access='stream', action='write', &
                  form='unformatted')

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! Open file
            open (350,file=trim(filename)// &
                  '/Solution-folder/JKQnuout-'//scoord, &
                  status='unknown', iostat=ios, err=1006, &
                  access='stream', action='write', &
                  form='unformatted')

          end if ! 1D vs 1.5D
        end if ! Saving JKQnu

      ! If there is suffix
      else

        ! Saving Stokes, open 1D file
        if (saveS) &
          open (250,file=trim(filename)//'/Stokesout_'//suff, &
                status='unknown', iostat=ios, err=1001, &
                access='stream', action='write', &
                form='unformatted')

        ! Saving JKQnu file, open 1D file
        if (saveJKQnu) &
          open (350,file=trim(filename)//'/JKQnuout_'//suff, &
                status='unknown', iostat=ios, err=1006, &
                access='stream', action='write', &
                form='unformatted')

      end if ! suffix

      ! Saving Stokes
      if (saveS) then

        ! Write flag and dimensions
        write(250,err=1101) 'bo'
        write(250,err=1101) Nfreq

      end if

      ! Saving JKQnu
      if (saveJKQnu) then

        ! Write flag and dimensions
        write(350,err=1106) 'ko'
        write(350,err=1106) nz
        write(350,err=1106) nfreq

      end if

      !
      ! Write data
      !

      ! Frequency axis
      if (saveS) write(250,err=1101) omega
      if (saveJKQnu) write(350,err=1106) omega

      ! Stokes out file
      if (saveS) then

        ! Number of directions
        write(250,err=1101) Geom%nTh/2,Geom%nPh

        ! For each polar direction
        do ith=1,Geom%nTh

          ! Ignore the ones going down
          if (Geom%V_mu(ith).lt.0) cycle

          ! For each azimuthal direction
          do iph=1,Geom%nPh

            ! Write the angles (DEG) of this quadrature direction
            write(250,err=1101) Geom%V_theta(ith)*180D0/pi, &
                                Geom%V_phi(iph)*180D0/pi

            ! Write the emergent Stokes parameters
            write(250,err=1101) transpose(Stokes(:,:,iph,ith,giz0))

          end do ! azimuthal directions
        end do ! polar directions

      end if ! Stokes out file

      ! JKQnu file
      if (saveJKQnu) then

        ! For each height
        do iz=1,nz

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each degree of freedom
            do ifreq=1,nfreq
              do K=0,2
                do iQ=-K,K

                  ! Write 0
                  write(350,err=1106) dzero,dzero

                end do
              end do
            end do

          ! In bounds
          else

            ! For each frequency
            do ifreq=1,nfreq

              ! For each K multipole
              do K=0,2

                ! For each Q multipole
                do iQ=-K,K

                  ! Write
                  write(350,err=1106) JKQC(iQ,K,ifreq,iz)

                end do ! Q
              end do ! K
            end do ! Frequency

          end if ! Height bounds

        end do ! Height

      end if ! JKQnu file

      ! Close opened files
      if (saveS) close(250)
      if (saveJKQnu) close(350)

      ! Free
      if (allocated(buffer)) deallocate(buffer)

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing solution file'
      close(200)
      inquire(unit=300, opened=laux)
      if (laux) close(300)
      inquire(unit=400, opened=laux)
      if (laux) close(400)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening Stokesout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing Stokesout file'
      close(250)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1002  umsg = 'Error opening Jout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1102  umsg = 'Error writing Jout file'
      close(300)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1003  umsg = 'Error opening Rout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1103  umsg = 'Error writing Rout file'
      close(400)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1004  umsg = 'Error opening Population file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1014  umsg = 'Error seeking Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1104  umsg = 'Error writing Population file'
      close(500)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1304  umsg = 'Error writing Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1005  umsg = 'Error opening Departure file'
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1015  umsg = 'Error seeking Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1105  umsg = 'Error writing Departure file'
      close(600)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1305  umsg = 'Error writing Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1006  umsg = 'Error opening JKQnuout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1106  umsg = 'Error writing JKQnuout file'
      close(350)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writesol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the self-consistent solution for the intensity problem
      !! in a file\n
      !!       Input(Input_class): Structure with configuration data\n
      !!       suff(character(:)): Suffix for the file names of the
      !!                           radiation field and density matrix
      !!                           files\n
      !!         omega(double(:)): Frequency array\n
      !!     Geom(Geometry_class): Structure with geometric data\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!             z(double(:)): Height array\n
      !!  Stokes(double(:,:,:,:)): Intensity\n
      !!         J00(double(:,:)): Mean intensity integrated over the
      !!                           absorption profile\n
      !!        J00S(double(:,:)): Mean intensity integrated over the
      !!                           emission profile\n
      !!        J00C(double(:,:)): Mean intensity with frequency
      !!                           dependence\n
      !!      J00P(double(:,:,:)): Intensity integrals in the
      !!                           photoionization rates\n
      !!            keep(logical): If the solution file should have a
      !!                           different name to avoid being
      !!                           overwriten by the polarized
      !!                           solution
      subroutine writesolI(Input,suff,omega,Geom,Atom,z, &
                           Stokes,J00,J00S,J00C,J00P,keep)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Input_class), intent(in):: Input
      type(Geometry_class), intent(in):: Geom
      character(len=4), intent(in):: suff
      logical, intent(in):: keep
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(nfreq,Geom%nPh,Geom%nTh,giz0:giz1),&
                        intent(in):: Stokes
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: J00, J00S
      double precision, dimension(nfreq,Rz0:Rz1), intent(in):: J00C
      double precision, dimension(nxphot,2,Rz0:Rz1), intent(in):: J00P

      ! Local

      integer, parameter:: izero = 0
      integer, parameter:: ione = 1

      double precision, parameter:: dzero = 0d0

      character(len=8):: scoord
      character(len=500):: filename

      logical:: saveJ00nu, saverKQ, saveJKQ
      logical:: saveP, saveD, saveS, saveSol
      logical:: laux

      integer:: ierr,ii,iab,K,iQ,iR,it,iJ,iran,ios,ia,iz,ifreq,i
      integer:: axial_int,stm_int,AV_int,ith,iph,itran,jtran
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: rJ,loffset


      ! Slaves
      if (pid.gt.0) then

        ! Control
        call control
        return

      end if ! Slaves

      !
      ! Shorter name and compute bools
      !
      filename = Input%folder
      saveSol = Input%keep_sol
      saveP = Input%keep_pop.and.nA.gt.0.and. &
              (suff.eq.'NONE'.or.run_mode.eq.0)
      saveD = Input%keep_dep.and.nA.gt.0.and. &
              (suff.eq.'NONE'.or.run_mode.eq.0)
      saverKQ = Input%keep_rhoKQ.and.nA.gt.0.and. &
                (suff.eq.'NONE'.or.run_mode.eq.0)
      saveJKQ = Input%keep_JKQ.and.nA.gt.0.and. &
                (suff.eq.'NONE'.or.run_mode.eq.0)
      saveS = Input%keep_stokesQ.and.(suff.eq.'NONE'.or.run_mode.eq.0)
      saveJ00nu = Input%keep_jkqnu.and. &
                  (suff.eq.'NONE'.or.run_mode.eq.0)

      ! If not writing anything
      if (.not.(saveSol.or.saveP.or.saveD.or.saverKQ.or.saveJKQ.or. &
                saveS.or.saveJ00nu)) then

        ! Leave
        call control
        return

      end if ! Not writing

      ! Routine name
      urou = 'writesolI'

      !
      ! Open files
      !

      ! Get LOS index string if 1.5D synthesis
      if (run_mode.eq.1) write(scoord,'(I0.8)') icoords(3)

      ! If writing solution
      if (saveSol) then

        ! If 1D
        if (run_mode.eq.0) then

          ! If different name
          if (keep) then

            ! Open solution file
            open (200,file=trim(filename)//'/SolutionI', &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')

          ! Usual name
          else

            ! Open solution file
            open (200,file=trim(filename)//'/Solution', &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')

          end if ! File name

        ! If 1.5D
        else if (run_mode.eq.1) then

          ! If different name
          if (keep) then

            ! Open solution file
            open (200,file=trim(filename)// &
                  '/Solution-folder/SolutionI-'//scoord, &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')

          ! Usual name
          else

            ! Open solution file
            open (200,file=trim(filename)// &
                  '/Solution-folder/Solution-'//scoord, &
                  status='unknown', iostat=ios, err=1000, &
                  access='stream', action='write', form='unformatted')

          end if ! File name

        end if ! 1D vs 1.5D
      end if ! Saving solution file

      ! If saving JKQ
      if (saveJKQ) then

        ! To write the final JKQ
        if (suff.eq.'NONE') then

          ! If 1D
          if (run_mode.eq.0) then

            ! If different name
            if (keep) then

              ! Open file
              open (300,file=trim(filename)//'/JoutI', &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')

            ! If 1.5D
            else

              ! Open file
              open (300,file=trim(filename)//'/Jout', &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! If different name
            if (keep) then

              ! Open file
              open (300,file=trim(filename)// &
                    '/Solution-folder/JoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (300,file=trim(filename)// &
                    '/Solution-folder/Jout-'//scoord, &
                    status='unknown', iostat=ios, err=1002, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name
          end if ! 1D vs 1.5D

        ! Not final JKQ
        else

          ! Open file
          open (300,file=trim(filename)//'/Jout_'//suff, &
                status='unknown', iostat=ios, err=1002, &
                access='stream', action='write', &
                form='unformatted')

        end if ! Final solution
      end if ! Saving JKQ

      ! If saving rhoKQ
      if (saverKQ) then

        ! To write the final rhoKQ
        if (suff.eq.'NONE') then

          ! If 1D
          if (run_mode.eq.0) then

            ! If different name
            if (keep) then

              ! Open file
              open (400,file=trim(filename)//'/RhooutI', &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (400,file=trim(filename)//'/Rhoout', &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! If different name
            if (keep) then

              ! Open file
              open (400,file=trim(filename)// &
                    '/Solution-folder/RhooutI-'//scoord, &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (400,file=trim(filename)// &
                    '/Solution-folder/Rhoout-'//scoord, &
                    status='unknown', iostat=ios, err=1003, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name
          end if ! 1D vs 1.5D

        ! Not final rhoKQ
        else

          ! Open file
          open (400,file=trim(filename)//'/Rhoout_'//suff, &
                status='unknown', iostat=ios, err=1003, &
                access='stream', action='write', &
                form='unformatted')

        end if ! Final solution
      end if ! Saving rhoKQ

      !
      ! Convert logicals to integers
      !

      ! Axial symmetry
      if(Geom%axial)then
          axial_int = 1
      else
          axial_int = 0
      end if

      ! Stimulated emission
      if(stm)then
          stm_int = 1
      else
          stm_int = 0
      end if

      ! Angle averaged redistribution function
      if (KSTK) then
          AV_int = 0
      else
          AV_int = 1
      end if

      !
      ! Write headers with dimensions and flags
      !

      ! If saving solution
      if (saveSol) then

        ! Solution file metadata
        write(200,err=1100) 'si'
        write(200,err=1100) nfreq
        write(200,err=1100) nZ
        write(200,err=1100) Geom%nTh
        write(200,err=1100) Geom%nPh
        write(200,err=1100) nA
        write(200,err=1100) axial_int
        write(200,err=1100) stm_int
        write(200,err=1100) AV_int

      end if ! Saving solution

      ! If saving JKQ
      if (saveJKQ) then

        ! JKQ file metadata
        write(300,err=1102) 'bj'
        write(300,err=1102) stm
        write(300,err=1102) nZ
        write(300,err=1102) nA
        write(300,err=1102) nxt
        write(300,err=1102) z

      end if ! Saving JKQ

      ! If saving rhoKQ
      if (saverKQ) then

        ! rhoKQ file metadata
        write(400,err=1103) 'br'
        write(400,err=1103) nZ
        write(400,err=1103) nA
        write(400,err=1103) z

      end if ! Saving rhoKQ

      ! If 1.5D
      if (run_mode.eq.1) then

        ! If saving population or departure, prepare buffers
        if (saveP.or.saveD) &
          allocate(buffer(maxval(Input%lim_pop%nbuff)/4))

      end if ! 1.5D synthesis

      !
      ! Write the data
      !

      ! Only if saving anything
      if (saveSol.or.saverKQ.or.saveP.or.saveD) then

        !
        ! Population and rhoKQ

        ! For each atom
        do ia=1,nA

          !
          ! If saving population or departure coeff.
          !
          if (saveP.or.saveD) then

            !
            ! If 1D
            !
            if (run_mode.eq.0) then

              ! If writing populations
              if (saveP) then

                ! To write the populations
                open (500,file=trim(filename)//'/'// &
                      trim(Atom(ia)%file_label)//'.pop', &
                      status='unknown',iostat=ios, err=1004, &
                      access='stream', action='write', &
                      form='unformatted')

                ! Metadata
                write(500,err=1104) 'bp'
                write(500,err=1104) nZ
                write(500,err=1104) Atom(ia)%nlevel

              end if ! Saving populations

              ! If writing departure c.
              if (saveD) then

                ! To write departure c.
                open (600,file=trim(filename)//'/'// &
                      trim(Atom(ia)%file_label)//'.dep', &
                      status='unknown',iostat=ios, err=1005, &
                      access='stream', action='write', &
                      form='unformatted')

                ! Metadata
                write(600,err=1105) 'bb'
                write(600,err=1105) nZ
                write(600,err=1105) Atom(ia)%nlevel

              end if ! Saving departure coeff.

              ! For each height
              do iz=1,nZ

                ! If out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! For each level
                  do i=1,Atom(ia)%nlevel

                    ! Write populations if requested
                    if (saveP) &
                      write(500,err=1104) Atom(ia)%popu(i,iz)

                    ! Write departure coeff if requested
                    if (saveD) &
                      write(600,err=1105) Atom(ia)%popu(i,iz)/ &
                                           Atom(ia)%populte(i,iz)
                  end do ! Levels

                ! If in bounds
                else

                  ! For each level
                  do i=1,Atom(ia)%nlevel

                    ! Get J and KQ index
                    rJ = Atom(ia)%rJval(Atom(ia)%sublevel(i), &
                         Atom(ia)%term(i))
                    iR = Atom(ia)%irho(Atom(ia)%term(i))% &
                                  Jrho(Atom(ia)%sublevel(i), &
                                       Atom(ia)%sublevel(i))%kq(0,0)

                    ! Write populations if requested
                    if (saveP) &
                      write(500,err=1104) Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz))

                    ! Write departure coeff if requested
                    if (saveD) &
                      write(600,err=1105) Atom(ia)%n(iz)* &
                                         sqrt(2d0*rJ+1d0)* &
                                         dble(Atom(ia)%crho(iR,iz))/ &
                                         Atom(ia)%populte(i,iz)
                  end do ! Levels

                end if ! Height bounds

              end do ! Heights

              ! Close files if opened
              if (saveP) close(500)
              if (saveD) close(600)

            !
            ! If 1.5D
            !
            else if (run_mode.eq.1) then

              ! Populations
              if (saveP.and.Input%lim_pop%nbuff(ia).gt.0) then

                ! Open file to write the populations
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.pop', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
                if (ierr.ne.0) goto 1004

                !
                ! Column offset
                !

                ! Get offset
                loffset = dble(icoords(3)-1)* &
                          dble(Input%lim_pop%nbuff(ia)) + &
                          dble(Input%lim_pop%head_size)
                do while(loffset.gt.offlimit)
                  offset = int(offlimit)
                  call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                  if (ierr.ne.0) goto 1014
                  loffset = loffset - offlimit
                end do
                offset = int(loffset)
                call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                if (ierr.ne.0) goto 1014

                ! Initialize buffer
                ii = 0

                ! If specified
                if (Input%lim_pop%nran.gt.0) then

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! This atom not included, skip
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz))

                      end do ! Ranges to print

                    ! In bounds
                    else

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz)))

                      end do ! Ranges to print

                    end if ! Height bounds

                  end do ! Heights

                ! Everything
                else

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each level
                      do i=1,Atom(ia)%nlevel

                        ! Write populations
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz))

                      end do ! Levels

                    ! In bounds
                    else

                      ! For each level
                      do i=1,Atom(ia)%nlevel

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance index
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                          sqrt(2d0*rJ+1d0)* &
                                          dble(Atom(ia)%crho(iR,iz)))
                      end do ! Levels

                    end if ! Height bounds

                  end do ! Heights

                end if ! Specific or everything

                ! Write buffer
                call MPI_FILE_WRITE(funit,buffer(1), &
                                    Input%lim_pop%nbuff(ia)/4, &
                                    MPI_REAL,MPI_STATUS_IGNORE,ierr)
                if (ierr.ne.0) goto 1304

              end if ! Saving populations

              ! Departure coefficients
              if (saveD.and.Input%lim_pop%nbuff(ia).gt.0) then

                ! Open file to write the populations
                call MPI_FILE_OPEN(MPI_COMM_SELF, &
                                 trim(filename)//'/'// &
                                 trim(Atom(ia)%file_label)//'.dep', &
                                 MPI_MODE_WRONLY, MPI_INFO_NULL, &
                                 funit, ierr)
                if (ierr.ne.0) goto 1005

                !
                ! Column offset
                !

                ! Get offset
                loffset = dble(icoords(3)-1)* &
                          dble(Input%lim_pop%nbuff(ia)) + &
                          dble(Input%lim_pop%head_size)
                do while(loffset.gt.offlimit)
                  offset = int(offlimit)
                  call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                  if (ierr.ne.0) goto 1015
                  loffset = loffset - offlimit
                end do
                offset = int(loffset)
                call MPI_FILE_SEEK(funit,offset,MPI_SEEK_CUR,ierr)
                if (ierr.ne.0) goto 1015

                ! Initialize buffer
                ii = 0

                ! If specified
                if (Input%lim_pop%nran.gt.0) then

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                          Atom(ia)%populte(i,iz))

                      end do ! Ranges to print

                    ! In bounds
                    else

                      ! For each entry to write
                      do iran=1,Input%lim_pop%nran

                        ! Atom
                        iab = Input%lim_pop%indx(1,iran)

                        ! Skip if not this atom
                        if (ia.ne.iab) cycle

                        ! Level
                        i = Input%lim_pop%indx(2,iran)

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance buffer
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(ii)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))

                      end do ! Ranges to print

                    end if ! Height bounds

                  end do ! Heights

                ! Everything
                else

                  ! For each height
                  do iz=1,nz

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! For each level
                      do i=1,Atom(ia)%nlevel

                        ! Advance index
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%popu(i,iz)/ &
                                          Atom(ia)%populte(i,iz))

                      end do ! Levels

                    ! In bounds
                    else

                      ! For each level
                      do i=1,Atom(ia)%nlevel

                        ! Get necessary data
                        it = Atom(ia)%term(i)
                        iJ = Atom(ia)%sublevel(i)
                        rJ = Atom(ia)%rJval(iJ,it)
                        iR = Atom(ia)%irho(it)%Jrho(iJ,iJ)%kq(0,0)

                        ! Advance index
                        ii = ii + 1

                        ! Save
                        buffer(ii) = real(Atom(ia)%n(iz)* &
                                        sqrt(2d0*rJ+1d0)* &
                                        dble(Atom(ia)%crho(iR,iz))/ &
                                        Atom(ia)%populte(i,iz))
                      end do ! Levels

                    end if ! Height bounds

                  end do ! Heights

                end if ! Specific or everything

                ! Write buffer
                call MPI_FILE_WRITE(funit,buffer(1), &
                                    Input%lim_pop%nbuff(ia)/4, &
                                    MPI_REAL,MPI_STATUS_IGNORE,ierr)
                if (ierr.ne.0) goto 1305

              end if ! Saving populations
            end if! 1D vs 1.5D
          end if ! Saving populations or departure coefficients

          ! Write the population of the atom into solution file
          if (saveSol) write(200,err=1100) Atom(ia)%n

          ! Saving rhoKQ
          if (saverKQ) then

            ! Write the population of the atom into rhoKQ file
            write(400,err=1100) Atom(ia)%n

            ! Write number of terms to rhoKQ file
            write(400,err=1103) Atom(ia)%nlevel

          end if

          ! For each term
          do i=1,Atom(ia)%nlevel

            ! Write number of levels to rhoKQ file
            if (saverKQ) write(400,err=1103) 1

            ! Get J value and level index
            rJ = Atom(ia)%rJval(Atom(ia)%sublevel(i),Atom(ia)%term(i))
            iR = Atom(ia)%irho(Atom(ia)%term(i))% &
                          Jrho(Atom(ia)%sublevel(i), &
                               Atom(ia)%sublevel(i))%kq(0,0)

            ! Write J values in rhoKQ file
            if (saverKQ) then
              write(400,err=1103) nint(2d0*rJ)
              write(400,err=1103) nint(2d0*rJ)
            end if

            ! For each K
            do K=0,nint(2d0*rJ)

              ! For each Q
              do iQ=-K,K

                ! If population
                if (K.lt.1) then

                  ! For each height
                  do iz=1,nZ

                    ! If out of bounds
                    if (iz.lt.Rz0.or.iz.gt.Rz1) then

                      ! Write zero rhoKQ into rhoKQ file, and the null
                      ! flag
                      if (saverKQ) then
                        write(400,err=1103) dzero,dzero
                        write(400,err=1103) ione
                      end if

                      ! Write zero rhoKQ into solution file
                      if (saveSol) &
                        write(200,err=1100) dzero

                    ! In bounds
                    else

                      ! If saving rhoKQ
                      if (saverKQ) then

                        ! Write rhoKQ into rhoKQ file, and the null
                        ! flag
                        write(400,err=1103) dble(Atom(ia)%crho(iR,iz))
                        write(400,err=1103) dzero

                        ! If null
                        if (Atom(ia)%rhonull(iR,iz)) then

                          ! Write one
                          write(400,err=1103) ione

                        ! Not null
                        else

                          ! Write zero
                          write(400,err=1103) izero

                        end if ! If rho00 is null
                      end if ! If writing rhoKQ

                      ! Write rhoKQ into solution file if requested
                      if (saveSol) &
                        write(200,err=1100) dble(Atom(ia)%crho(iR,iz))

                    end if ! Height bounds

                  end do ! heights

                ! K>0 multipole
                else

                  ! For each height
                  do iz=1,nZ

                    ! If saving rhoKQ
                    if (saverKQ) then

                      ! Write rhoKQ into rhoKQ file, and the null flag
                      write(400,err=1103) dzero
                      write(400,err=1103) dzero
                      write(400,err=1103) ione

                    end if ! Saving rhoKQ

                  end do ! heights

                end if ! Multipole

              end do ! Q
            end do ! K
          end do ! Levels
        end do ! Atoms

      end if ! Saving rho or Solution


      !
      ! J00

      ! Only if saving anything in this block
      if (saveSol.or.saveJKQ) then

        ! For each atom
        do ia=1,nA

          ! Write number of transitions in JKQ file
          if (saveJKQ) write(300,err=1102) Atom(ia)%nftran

          ! For each transition
          do itran=1,Atom(ia)%nftran

            ! Apply atomic shift
            jtran = itran + Atom(ia)%tfshift

            ! For each height
            do iz=1,nZ

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Write zero JKQ into the JKQ file
                if (saveJKQ) write(300,err=1102) dzero,dzero

                ! Write zero JKQ into the solution file
                if (saveSol) write(200,err=1100) dzero

              ! In bounds
              else

                ! Write the JKQ into the JKQ file
                if (saveJKQ) write(300,err=1102) J00(jtran,iz),dzero

                ! Write the JKQ into the solution file
                if (saveSol) write(200,err=1100) J00(jtran,iz)

              end if ! Height bounds

            end do ! heights

            ! For each K>0
            do K=1,2

              ! For each Q
              do iQ=-K,K

                ! For each height
                do iz=1,nZ

                  ! Write the JKQ into the JKQ file
                  if (saveJKQ) write(300,err=1102) dzero,dzero

                end do ! heights
              end do ! Q
            end do ! K
          end do ! transitions
        end do ! Atoms

        ! If there is stimulated emission, write JKQS
        if(stm)then

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%nftran

              ! Apply atomic shift
              jtran = itran + Atom(ia)%tfshift

              ! For each height
              do iz=1,nZ

                ! If out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! Write the JKQ into the JKQ file
                  if (saveJKQ) write(300,err=1102) dzero,dzero

                  ! Write the JKQ into the solution file
                  if (saveSol) write(200,err=1100) dzero

                ! In bounds
                else

                  ! Write the JKQ into the JKQ file
                  if (saveJKQ) write(300,err=1102) J00S(jtran,iz), &
                                                   dzero

                  ! Write the JKQ into the solution file
                  if (saveSol) write(200,err=1100) J00S(jtran,iz)

                end if ! Height bounds

              end do ! heights

              ! Saving JKQ
              if (saveJKQ) then

                ! For each K
                do K=1,2

                  ! For each Q
                  do iQ=-K,K

                    ! For each height
                    do iz=1,nZ

                      ! Write the JKQS into the JKQ file
                      write(300,err=1102) dzero,dzero

                    end do ! heights
                  end do ! Q
                end do ! K

              end if ! Saving JKQ

            end do ! transitions
          end do ! atoms

        end if ! stimulated emission
      end if ! Saving anything

      !
      ! J00P

      ! Saving solution
      if (saveSol) then

        ! For each atom
        do ia=1,nA

          ! For each transition
          do itran=1,Atom(ia)%nphot

            ! Apply atomic shift
            jtran = itran + Atom(ia)%pshift

            ! For each height
            do iz=1,nZ

              ! If out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Write zero
                write(200,err=1100) dzero,dzero

              ! In bounds
              else

                ! Write the J00P into the solution file
                write(200,err=1100) J00P(jtran,1,iz)
                write(200,err=1100) J00P(jtran,2,iz)

              end if ! Height bounds

            end do ! heights

          end do ! transitions
        end do ! Atoms


        !
        ! Radiation field solution
        !

        !
        ! If not angle-averaged
        if (KSTK) then

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each degree of freedom
              do ith=1,Geom%nTh
                do iph=1,Geom%nPh
                  do ifreq=1,nfreq

                    ! Write 0
                    write(200,err=1100) dzero

                  end do
                end do
              end do

            ! In bounds
            else

              ! For each polar direction
              do ith=1,Geom%nTh

                ! For each azimuthal direction
                do iph=1,Geom%nPh

                  ! For each frequency
                  do ifreq=1,nfreq

                    ! Write Stokes parameter into the solution file
                    write(200,err=1100) Stokes(ifreq,iph,ith,iz)

                  end do ! Frequencies
                end do ! azimuthal directions
              end do ! polar directions

            end if ! Height bounds

          end do ! heights

        !
        ! If Angle-averaged
        else

          ! For each height
          do iz=1,nZ

            ! If out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each frequency
              do ifreq=1,nfreq

                ! Write 0
                write(200,err=1100) dzero

              end do ! frequencies

            ! In bounds
            else

              ! For each frequency
              do ifreq=1,nfreq

                ! Write the JKQ(k) into the solution file
                write(200,err=1100) J00C(ifreq,iz)

              end do ! frequencies

            end if ! Height bounds

          end do ! heights

        end if ! AV or AD

        !
        ! Close files
        !
        close (200)

      end if ! Saving Solution

      ! Close files if opened
      if (saveJKQ) close (300)
      if (saverKQ) close (400)

      ! If not storing anything else
      if (.not.saveJ00nu.and..not.saveS) then

        ! Free
        if (allocated(buffer)) deallocate(buffer)

        ! Control
        call control
        return

      end if ! Not storing anything else


      !
      ! Store the stokes in the quadrature
      !

      ! If there is no suffix
      if (suff.eq.'NONE') then

        ! Saving Stokes
        if (saveS) then

          ! If 1D
          if (run_mode.eq.0) then

            ! Different name
            if (keep) then

              ! Open file
              open (250,file=trim(filename)//'/StokesoutI', &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (250,file=trim(filename)//'/Stokesout', &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! Different name
            if (keep) then

              ! Open file
              open (250,file=trim(filename)// &
                    '/Solution-folder/StokesoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (250,file=trim(filename)// &
                    '/Solution-folder/Stokesout-'//scoord, &
                    status='unknown', iostat=ios, err=1001, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name
          end if ! 1D vs 1.5D
        end if ! Saving Stokes

        ! JKQnu file
        if (saveJ00nu) then

          ! If 1D
          if (run_mode.eq.0) then

            ! Different name
            if (keep) then

              ! Open file
              open (350,file=trim(filename)//'/JKQnuoutI', &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (350,file=trim(filename)//'/JKQnuout', &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name

          ! If 1.5D
          else if (run_mode.eq.1) then

            ! Different name
            if (keep) then

              ! Open file
              open (350,file=trim(filename)// &
                    '/Solution-folder/JKQnuoutI-'//scoord, &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            ! Usual name
            else

              ! Open file
              open (350,file=trim(filename)// &
                    '/Solution-folder/JKQnuout-'//scoord, &
                    status='unknown', iostat=ios, err=1006, &
                    access='stream', action='write', &
                    form='unformatted')

            end if ! File name
          end if ! 1D vs 1.5D
        end if ! Saving J00nu

      ! If there is suffix
      else

        ! Saving Stokes, open 1D file
        if (saveS) &
          open (250,file=trim(filename)//'/Stokesout_'//suff, &
                status='unknown', iostat=ios, err=1001, &
                access='stream', action='write', &
                form='unformatted')

        ! Saving JKQnu file, open 1D file
        if (saveJ00nu) &
          open (350,file=trim(filename)//'/JKQnuout_'//suff, &
                status='unknown', iostat=ios, err=1006, &
                access='stream', action='write', &
                form='unformatted')

      end if ! suffix

      ! Saving Stokes
      if (saveS) then

        ! Write flag and dimensions
        write(250,err=1101) 'bo'
        write(250,err=1101) Nfreq

      end if

      ! Saving JKQnu
      if (saveJ00nu) then

        write(350,err=1106) 'ko'
        write(350,err=1106) nz
        write(350,err=1106) nfreq

      end if

      !
      ! Write data
      !

      ! Frequency axis
      if (saveS) write(250,err=1101) omega
      if (saveJ00nu) write(350,err=1106) omega

      ! Saving Stokes
      if (saveS) then

        ! Number of directions
        write(250,err=1101) Geom%nTh/2,Geom%nPh

        ! For each polar direction
        do ith=1,Geom%nTh

          ! Ignore the ones going down
          if (Geom%V_mu(ith).lt.0) cycle

          ! For each azimuthal direction
          do iph=1,Geom%nPh

            ! Write the angles (DEG) of this quadrature direction
            write(250,err=1101) Geom%V_theta(ith)*180D0/pi, &
                                Geom%V_phi(iph)*180D0/pi

            ! Write the emergent Stokes parameters
            write(250,err=1101) Stokes(:,iph,ith,giz0)

            ! Fill polarization with zeros
            do ifreq=1,nfreq*3
              write(250,err=1101) dzero
            end do

          end do ! azimuthal directions
        end do ! polar directions

      end if ! Saving Stokes

      ! Save J00nu
      if (saveJ00nu) then

        ! For each height
        do iz=1,nz

          ! If out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! For each frequency
            do ifreq=1,nfreq

              ! Fill with zero
              write(350,err=1106) dzero

            end do ! Frequencies

          ! In bounds
          else

            ! Write J00C
            write(350,err=1106) J00C(:,iz)

          end if ! Height bounds

        end do ! Heights

      end if

      ! Close the Stokesout file
      if (saveS) close(250)
      if (saveJ00nu) close(350)

      ! Free
      if (allocated(buffer)) deallocate(buffer)

      ! Control
      call control

      return

1000  umsg = 'Error opening solution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing solution file'
      close(200)
      inquire(unit=300, opened=laux)
      if (laux) close(300)
      inquire(unit=400, opened=laux)
      if (laux) close(400)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening Stokesout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing Stokesout file'
      close(250)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1002  umsg = 'Error opening Jout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1102  umsg = 'Error writing Jout file'
      close(300)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1003  umsg = 'Error opening Rout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1103  umsg = 'Error writing Rout file'
      close(400)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1004  umsg = 'Error opening Population file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1014  umsg = 'Error seeking Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1104  umsg = 'Error writing Population file'
      close(500)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1304  umsg = 'Error writing Population file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1005  umsg = 'Error opening Departure file'
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1015  umsg = 'Error seeking Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1105  umsg = 'Error writing Departure file'
      close(600)
      if (saveSol) close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1305  umsg = 'Error writing Departure file'
      if (saveSol) close(200)
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1006  umsg = 'Error opening JKQnuout file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1106  umsg = 'Error writing JKQnuout file'
      close(350)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writesolI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the emergent Stokes parameters in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(double(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!       Stokes(double(:,:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writestk(filename,iph,ith,omega,Geom,Stokes,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:,:), intent(in):: Stokes

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,jj,i0,i1,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writestk'

      ! Convert the integers into appropriate length strings
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        ! Open file
        open(200,file=trim(filename)//'/Stokes_'//trim(cth)//'_'// &
             trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')

        ! Identification
        write(200,err=1100) 'be'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Stokes vector
        write(200,err=1100) transpose(Stokes)

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each Stokes parameter
          do jj=1,4

            ! For each entry to write
            do iran=1,buff%nran

              ! Range and size
              i0 = buff%indx(1,iran)
              i1 = buff%indx(2,iran)
              nn = i1-i0+1

              ! Fill buffer
              buffer(ii+1:ii+nn) = real(Stokes(jj,i0:i1))

              ! Advance index
              ii = ii + nn

            end do ! Ranges
          end do ! Stokes parameters

        ! All
        else

          ! Size
          nn = nfreq*4

          ! Fill buffer
          buffer = real(reshape(transpose(Stokes), (/ nn /)))

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      !
      ! If inversion
      !
      else if (run_mode.eq.-1) then

        ! Size
        nn = buff%buffer_size/4

        ! Allocate buffer
        allocate(buffer(nn))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes', &
                           MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Fill buffer
        buffer = real(reshape(transpose(Stokes), (/ nn /)))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing Stokes file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writestk

!#####################################################################
!#####################################################################
!#####################################################################

      !! Save the emergent intensity in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(double(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!        StokesI(double(:)): Intensity\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writestkI(filename,iph,ith,omega,Geom,StokesI,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: StokesI

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,i0,i1,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset
      double precision, dimension(nfreq,0:3):: Stokes


      ! Routine name
      urou = 'writestkI'

      ! Convert the integers into appropriate length strings
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        ! Convert intensity into polarization array
        Stokes(:,1:3) = 0d0
        Stokes(:,0) = StokesI

        ! Open file
        open(200,file=trim(filename)//'/StokesI_'// trim(cth)//'_'// &
             trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')

        ! Identification
        write(200,err=1100) 'be'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,Geom%L_theta(ith)*180d0/pi, &
                          Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Stokes vector
        write(200,err=1100) Stokes

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/StokesI_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size)/4 + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Range and size
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(StokesI(i0:i1))

            ! Advance index
            ii = ii + nn

          end do ! Ranges

        ! All
        else

          ! Fill buffer
          buffer(1:nfreq) = real(StokesI)

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/16, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      !
      ! If inversion
      !
      else if (run_mode.eq.-1) then

        ! Size
        nn = buff%buffer_size/16

        ! Allocate buffer
        allocate(buffer(nn))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Stokes', &
                           MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size)/4 + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Fill buffer
        buffer = real(StokesI)

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      return

1000  umsg = 'Error opening StokesI file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing StokesI file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing StokesI file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writestkI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the emergent Stokes parameters in RAM\n
      !!     e_Stokes(double(:,:)): RAM storage for Stokes
      !!                            parameters\n
      !!       Stokes(double(:,:)): Stokes parameters\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!           lrange(logical): If the data in buff needs to be
      !!                            considered
      subroutine setstk(e_Stokes,Stokes,buff,lrange)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      logical, intent(in):: lrange
      double precision, dimension(:,:), intent(in):: Stokes
      double precision, dimension(:,:), intent(out):: e_Stokes

      ! Local

      integer:: iran,ii,jj,i0,i1,nn


      ! Slaves leave
      if (pid.gt.0) return

      ! If specified and need to be considered
      if (buff%nran.gt.0.and.lrange) then

        ! For each Stokes parameter
        do jj=1,4

          ! Initialize buffer
          ii = 0

          ! For each entry to write
          do iran=1,buff%nran

            ! Range and size
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            e_Stokes(jj,ii+1:ii+nn) = Stokes(jj,i0:i1)

            ! Advance index
            ii = ii + nn

          end do ! Ranges
        end do ! Stokes parameters

      ! All
      else

        ! Copy
        e_Stokes = Stokes

      end if

      return

      end subroutine setstk

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the emergent intensity in RAM.\n
      !!     e_Stokes(double(:,:)): Output stokes parameters\n
      !!         Stokes(double(:)): Intensity\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!           lrange(logical): If the data in buff needs to be
      !!                            considered
      subroutine setstkI(e_Stokes,Stokes,buff,lrange)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      logical, intent(in):: lrange
      double precision, dimension(:), intent(in):: Stokes
      double precision, dimension(:,:), intent(out):: e_Stokes

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Slaves leave
      if (pid.gt.0) return

      ! If specified and need to be considered
      if (buff%nran.gt.0.and.lrange) then

        ! Initialize buffer
        ii = 0

        ! For each entry to write
        do iran=1,buff%nran

          ! Range and size
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_Stokes(1,ii+1:ii+nn) = Stokes(i0:i1)

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Copy
        e_Stokes(1,:) = Stokes

      end if

      return

      end subroutine setstkI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write the position of the LOS in the CLE problem in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!             wtau(logical): If writing tau
      subroutine write_CLEgeom(filename,Atmo,buff,wtau)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      logical, intent(in):: wtau

      ! Local

      integer:: ierr

      integer(kind=MPI_OFFSET_KIND):: offset

      double precision, dimension(2):: buffer


      ! Cartesian does not need to write, leave
      if (Atmo%mode.eq.0) return

      ! If no master
      if (pid.gt.0) then

        ! Control
        call control
        return

      end if ! Slaves

      ! Routine name
      urou = 'write_CLEgeom'

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Stokes', &
                         MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Get offset and data
      !

      ! If slab or non-cartesian
      if (Atmo%mode.eq.1.or.Atmo%mode.eq.2) then

        ! Compute offset
        offset = (icoords(3)-1)*16 + buff%head_size

        ! Slab
        if (Atmo%mode.eq.1) then

          ! Get height and theta
          buffer(1) = Atmo%z(1)
          buffer(2) = Atmo%ypos

        ! Non-cartesian
        else if (Atmo%mode.eq.2) then

          ! Get y and z position
          buffer(1) = Atmo%ypos
          buffer(2) = Atmo%zpos

        end if ! Type of model atmosphere
      end if ! Atmo mode

      ! Go to offset
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Write buffer
      call MPI_FILE_WRITE(funit, buffer(1), 2, &
                          MPI_DOUBLE_PRECISION, &
                          MPI_STATUS_IGNORE, ierr)
      if (ierr.ne.0) goto 1300

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)


      ! If writing tau
      if (wtau) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Tau', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1100


        ! Go to offset
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1110

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1), 2, &
                            MPI_DOUBLE_PRECISION,&
                            MPI_STATUS_IGNORE, ierr)
        if (ierr.ne.0) goto 1310

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! Tau

      ! Control
      call control

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error opening Tau file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1110  umsg = 'Error seeking Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1310  umsg = 'Error writing Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine write_CLEgeom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write the emergent Stokes and optical depth in the CLE
      !! problem in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              if0(integer): Lower limit index for frequency\n
      !!              if1(integer): Upper limit index for frequency\n
      !!               nf(integer): Total number of frequencies\n
      !!       Stokes(double(:,:)): Stokes parameters\n
      !!            tau(double(:)): Optical depth\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable\n
      !!             wtau(logical): If writing tau
      subroutine write_CLE(filename,if0,if1,nf,Stokes,tau, &
                           buff,wtau)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      logical, intent(in):: wtau
      integer, intent(in):: if0,if1,nf
      double precision, dimension(:), intent(in):: tau
      double precision, dimension(:,:), intent(in):: Stokes

      ! Local

      integer:: ierr,is,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'write_CLE'

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Stokes', &
                         MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size+buff%geom_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      !
      ! If splitting frequencies
      !
      if (nproc.gt.1) then

        ! Allocate buffer
        nn = nf
        allocate(buffer(nn))

        ! For each Stokes parameter
        do is=1,4

          ! Initial wavelength shift offset
          if (if0.gt.1) then
            loffset = dble((if0-1)*4)
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1010
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
          end if ! Not first wavelength

          ! Store in buffer
          buffer = real(Stokes(is,:))

          ! Write buffer
          call MPI_FILE_WRITE(funit,buffer(1),nn, &
                              MPI_REAL,MPI_STATUS_IGNORE,ierr)
          if (ierr.ne.0) goto 1300

          ! After V, skip rest
          if (is.eq.4) exit

          ! Final wavelength shift offset
          if (if1.lt.buff%nn) then
            loffset = dble((buff%nn-if1)*4)
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1010
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
          end if ! Not last wavelength

        end do ! Stokes parameter

      !
      ! Serial
      !
      else

        ! Allocate buffer
        nn = size(stokes)
        allocate(buffer(nn))

        ! Add to buffer
        buffer = real(reshape(transpose(Stokes),(/nn/)))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

      end if ! Serial/MPI (freqs)

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      ! Free
      deallocate(buffer)

      ! If writing tau
      if (wtau) then

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, trim(filename)//'/Tau', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1100

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size/4) + &
                  dble(buff%head_size+buff%geom_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1110
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1110

        ! Allocate buffer
        nn = nf

        ! Store in buffer
        buffer(1:nn) = real(tau)

        !
        ! If splitting frequencies
        !
        if (nproc.gt.1) then

          ! Initial wavelength shift offset
          if (if0.gt.1) then
            loffset = dble(if0-1)*4
            do while(loffset.gt.offlimit)
              offset = int(offlimit)
              call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
              if (ierr.ne.0) goto 1110
              loffset = loffset - offlimit
            end do
            offset = int(loffset)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1110
          end if ! Not first wavelength

        end if ! MPI (freqs)

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),nn, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1310

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

      end if ! Tau

      return

1000  umsg = 'Error opening Stokes file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing Stokes file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error opening Tau file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1110  umsg = 'Error seeking Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1310  umsg = 'Error writing Tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine write_CLE

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the contribution function of a synthesis run in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction
      !!          omega(double(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!              z(double(:)): Height array\n
      !!      Contr(double(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writectr(filename,iph,ith,omega,Geom,z,Contr,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(:,:,:), intent(in):: Contr

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,ifreq,iran,iz,jz,ii,jj,i0,i1,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectr'

      ! Convert the integers into appropriate length strings
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Contribution_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream',action='write',form='unformatted')

        ! Identification
        write(200,err=1100) 'bc'

        ! Number of frequencies, number of heights, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height axis
        write(200,err=1100) z

        ! Limited heights
        if (nz.ne.Rnz) then

          ! For each Stokes parameter
          do ios=1,4

            ! For each height less than lower limit, write zero
            do iz=1,Rz0-1

              ! For each frequency
              do ifreq=1,nfreq

                ! Write zero
                write(200,err=1100) 0d0

              end do ! Frequency
            end do ! Height

            ! Write contribution function
            write(200,err=1100) Contr(ios,:,:)

            ! For each height larger than upper limit, write zero
            do iz=Rz1+1,nz

              ! For each frequency
              do ifreq=1,nfreq

                ! Write zero
                write(200,err=1100) 0d0

              end do ! Frequency
            end do ! Height
          end do ! For each Stokes parameter

        ! All heights
        else

          ! For each Stokes parameter
          do ios=1,4

            ! Write contribution function, order: is, iz, ifreq
            write(200,err=1100) Contr(ios,:,:)

          end do ! Stokes parameters

        end if ! Limited heights

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Contribution_'// &
                           trim(cth)//'_'//trim(cph), &
                           MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each Stokes parameter
          do jj=1,4

            ! For each height
            do iz=1,nz

              ! Out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! For each entry to write
                do iran=1,buff%nran

                  ! Range and size
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = real(0)

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges

              ! In bounds
              else

                ! Shift height
                jz = iz - Rz0 + 1

                ! For each entry to write
                do iran=1,buff%nran

                  ! Range and size
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = real(Contr(jj,i0:i1,jz))

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges

              end if ! Height bounds

            end do ! Heights
          end do ! Stokes parameters

        ! All
        else

          ! Not full height range
          if (nz.ne.Rnz) then

            ! For each Stokes parameter
            do jj=1,4

              ! For each height
              do iz=1,nz

                ! Out of bounds
                if (iz.lt.Rz0.or.iz.gt.Rz1) then

                  ! Fill buffer
                  buffer(ii+1:ii+nfreq) = real(0)

                ! In bounds
                else

                  ! Shift height
                  jz = iz - Rz0 + 1

                  ! Fill buffer
                  buffer(ii+1:ii+nfreq) = real(Contr(jj,:,jz))

                end if ! Height bounds

                ! Advance index
                ii = ii + nfreq

              end do ! Heights
            end do ! Stokes parameters

          ! Full height range
          else

            ! Size
            nn = nfreq*nz

            ! For each Stokes parameter
            do jj=1,4

              ! Fill buffer
              buffer(ii+1:ii+nn) = real(reshape(Contr(jj,:,:), &
                                        (/ nn /)))

              ! Advance index
              ii = ii + nn

            end do ! Stokes parameters

          end if ! Full height range
        end if ! Full frequency range

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writectr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the contribution function of an inversion run in a
      !! file\n
      !!    filename(character(:)): Name of the file to write\n
      !!        Contr(real(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writectr_inv(filename,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:,:,:), intent(in):: Contr

      ! Local

      integer:: ierr,iz,jz,ii,jj,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectr_inv'

      ! Size
      nn = buff%buffer_size/4

      ! Allocate buffer
      allocate(buffer(nn))

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Contribution', &
                         MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Initialize offset buffer
      ii = 0

      ! Not full height range
      if (nz.ne.Rnz) then

        ! For each Stokes parameter
        do jj=1,4

          ! For each height
          do iz=1,nz

            ! Out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! Fill buffer
              buffer(ii+1:ii+buff%nn) = real(0)

            ! In bounds
            else

              ! Shift height
              jz = iz - Rz0 + 1

              ! Fill buffer
              buffer(ii+1:ii+buff%nn) = Contr(jj,:,jz)

            end if ! Height bounds

            ! Advance index
            ii = ii + buff%nn

          end do ! Heights
        end do ! Stokes parameters

      ! Full height range
      else

        ! For each Stokes parameter
        do jj=1,4

          ! Fill buffer
          buffer(ii+1:ii+buff%nn*nz) = &
                              reshape(Contr(jj,:,:), (/ buff%nn*nz /))

          ! Advance index
          ii = ii + buff%nn*nz

        end do ! Stokes parameters

      end if ! Full height range

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),nn, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      ! Free
      deallocate(buffer)

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writectr_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the intensity contribution function of a synthesis run
      !! in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(double(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!              z(double(:)): Height array\n
      !!        Contr(double(:,:)): Intensity contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writectrI(filename,iph,ith,omega,Geom,z,Contr,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega,z
      double precision, dimension(:,:), intent(in):: Contr

      ! Local

      character(len=4):: cph,cth

      integer:: ios,iz,jz,ifreq,ierr,iran,ii,i0,i1,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectrI'

      ! Convert the integers into appropriate length strings
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        ! Open file
        open(200,file=trim(filename)//'/Contribution_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream',action='write',form='unformatted')

        ! Identification
        write(200,err=1100) 'bc'

        ! Number of frequencies, number of heights, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height axis
        write(200,err=1100) z

        ! Before lower limit
        do iz=1,Rz0-1

          ! For each frequency
          do ifreq=1,nfreq

            ! Write zero
            write(200,err=1100) 0d0

          end do ! Frequencies
        end do ! Height below lower limit

        ! Write contribution function, order: is, iz, ifreq
        write(200,err=1100) Contr

        ! Above upper limit
        do iz=Rz1+1,nZ

          ! For each frequency
          do ifreq=1,nfreq

            ! Write zero
            write(200,err=1100) 0d0

          end do ! Frequencies
        end do ! Height above upper limit

        ! Fill the other Stokes
        do ios=1,3

          ! For each height
          do iz=1,nZ

            ! For each frequency
            do ifreq=1,nfreq

              ! Write zero
              write(200,err=1100) 0d0

            end do ! Frequencies
          end do ! Heights
        end do ! Stokes parameters

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))
        buffer = 0e0

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Contribution_'// &
                           trim(cth)//'_'//trim(cph), &
                           MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! Heights
          do iz=1,nz

            ! Out of bounds
            if (iz.lt.Rz0.or.iz.gt.Rz1) then

              ! For each entry to write
              do iran=1,buff%nran

                ! Range and size
                i0 = buff%indx(1,iran)
                i1 = buff%indx(2,iran)
                nn = i1-i0+1

                ! Fill buffer
                buffer(ii+1:ii+nn) = real(0)

                ! Advance index
                ii = ii + nn

              end do ! Ranges

            ! In bounds
            else

              ! Shift height
              jz = iz - Rz0 + 1

              ! For each entry to write
              do iran=1,buff%nran

                ! Range and size
                i0 = buff%indx(1,iran)
                i1 = buff%indx(2,iran)
                nn = i1-i0+1

                ! Fill buffer
                buffer(ii+1:ii+nn) = real(Contr(i0:i1,jz))

                ! Advance index
                ii = ii + nn

              end do ! Ranges

            end if ! Height bounds

          end do ! Heights

        ! All
        else

          ! Limited range
          if (nz.ne.Rnz) then

            ! For each height
            do iz=1,nz

              ! Out of bounds
              if (iz.lt.Rz0.or.iz.gt.Rz1) then

                ! Fill buffer
                buffer(ii+1:ii+nfreq) = real(0)

              ! In bounds
              else

                ! Shift height
                jz = iz - Rz0 + 1

                ! Fill buffer
                buffer(ii+1:ii+nfreq) = real(Contr(:,jz))

              end if ! Height bounds

              ! Advance index
              ii = ii + nfreq

            end do ! Heights

          ! Full range
          else

            ! Size
            nn = nfreq*nz

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(reshape(Contr, (/ nn /)))

          end if ! Height range
        end if ! Limited ranges

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writectrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the intensity contribution function of an inversion run
      !! in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!          Contr(real(:,:)): Intensity contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writectrI_inv(filename,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:,:), intent(in):: Contr

      ! Local

      integer:: iz,jz,ierr,ii,nn,nt
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writectrI_inv'

      ! Size
      nt = buff%buffer_size/4
      nn = nt/4

      ! Allocate buffer
      allocate(buffer(nt))
      buffer = 0e0

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Contribution', &
                         MPI_MODE_WRONLY,MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Limited range
      if (nz.ne.Rnz) then

        ! Initialize buffer
        ii = 0

        ! For each height
        do iz=1,nz

          ! Out of bounds
          if (iz.lt.Rz0.or.iz.gt.Rz1) then

            ! Fill buffer with zeros
            buffer(ii+1:ii+buff%nn) = real(0)

          ! In bounds
          else

            ! Shift height
            jz = iz - Rz0 + 1

            ! Fill buffer
            buffer(ii+1:ii+buff%nn) = Contr(:,jz)

          end if ! Height bounds

          ! Advance index
          ii = ii + buff%nn

        end do ! Heights

      ! Full range
      else

        ! Fill buffer
        buffer(1:buff%nn*nz) = reshape(Contr, (/ buff%nn*nz /))

      end if ! Height range

      ! Write buffer
      call MPI_FILE_WRITE(funit,buffer(1),nt, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      ! Free
      deallocate(buffer)

      return

1000  umsg = 'Error opening contribution file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing contribution file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writectrI_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the contribution function in RAM\n
      !!    e_Contr(double(:,:,:)): RAM storage for contribution
      !!                            function\n
      !!      Contr(double(:,:,:)): Contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine setctr(e_Contr,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:,:,:), intent(in):: Contr
      real, dimension(:,:,:), intent(out):: e_Contr

      ! Local

      integer:: iran,ii,jj,i0,i1,nn


      ! Slaves leave
      if (pid.gt.0) return

      ! Zero out of bounds
      e_Contr(:,:,1:Rz0-1) = 0.0
      e_Contr(:,:,Rz1+1:nz) = 0.0

      ! If specified
      if (buff%nran.gt.0) then

        ! For each Stokes parameter
        do jj=1,4

          ! Initialize buffer
          ii = 0

          ! For each entry to write
          do iran=1,buff%nran

            ! Range and size
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            e_Contr(jj,ii+1:ii+nn,Rz0:Rz1) = &
                                         real(Contr(jj,i0:i1,1:Rnz))

            ! Advance index
            ii = ii + nn

          end do ! Ranges
        end do ! Stokes parameters

      ! All
      else

        ! Copy
        e_Contr(:,:,Rz0:Rz1) = real(Contr(:,:,1:Rnz))

      end if ! Full frequency range

      return

      end subroutine setctr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the intensity contribution function in RAM\n
      !!    e_Contr(double(:,:,:)): RAM storage for intensity
      !!                            contribution function\n
      !!      Contr(double(:,:,:)): Intensity contribution function\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine setctrI(e_Contr,Contr,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:,:), intent(in):: Contr
      real, dimension(:,:,:), intent(out):: e_Contr

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Only master
      if (pid.gt.0) return

      ! Zero out of bounds
      e_Contr(:,:,1:Rz0-1) = 0.0
      e_Contr(:,:,Rz1+1:nz) = 0.0

      ! If specified
      if (buff%nran.gt.0) then

        ! Initialize buffer
        ii = 0

        ! For each entry to write
        do iran=1,buff%nran

          ! Range and size
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_Contr(1,ii+1:ii+nn,Rz0:Rz1) = real(Contr(i0:i1,1:Rnz))

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Copy
        e_Contr(1,:,Rz0:Rz1) = real(Contr(:,1:Rnz))

      end if ! Full frequency range

      return

      end subroutine setctrI

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the height of optical depth equal to one of a synthesis
      !! run in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              iph(integer): Index of the LOS azimuth
      !!                            direction\n
      !!              ith(integer): Index of the LOS polar direction\n
      !!          omega(double(:)): Frequency array\n
      !!      Geom(Geometry_class): Structure with geometric data\n
      !!            tau(double(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writetau(filename,iph,ith,omega,Geom,tau,buff)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      integer, intent(in):: iph,ith
      double precision, dimension(:), intent(in):: omega
      double precision, dimension(:), intent(in):: tau

      ! Local

      character(len=4):: cph,cth

      integer:: ios,ierr,iran,ii,i0,i1,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Routine name
      urou = 'writetau'

      ! Convert the integers into appropriate length strings
      if (ith.lt.1000.and.ith.ge.100) write(cth,'(i3)') ith
      if (ith.lt.100 .and.ith.ge.10 ) write(cth,'(i2)') ith
      if (ith.lt.10  .and.ith.ge.0  ) write(cth,'(i1)') ith
      if (iph.lt.1000.and.iph.ge.100) write(cph,'(i3)') iph
      if (iph.lt.100 .and.iph.ge.10 ) write(cph,'(i2)') iph
      if (iph.lt.10  .and.iph.ge.0  ) write(cph,'(i1)') iph

      !
      ! If 1D
      !
      if (run_mode.eq.0) then

        !
        ! Open file
        !
        open(200,file=trim(filename)//'/Tau_'//trim(cth)// &
             '_'//trim(cph), status='unknown', iostat=ios, err=1000, &
             access='stream', action='write', form='unformatted')

        ! Identification
        write(200,err=1100) 'bt'

        ! Number of frequencies, angles of LOS
        write(200,err=1100) nfreq,nZ,Geom%L_theta(ith)*180d0/pi, &
                            Geom%L_phi(iph)*180d0/pi

        ! Frequency axis
        write(200,err=1100) omega

        ! Height of tau=1
        write(200,err=1100) tau

        !
        ! Close file
        !
        close(200)

      !
      ! If 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(filename)//'/Tau_'//trim(cth)// &
                           '_'//trim(cph),MPI_MODE_WRONLY, &
                           MPI_INFO_NULL,funit,ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Range and size
            i0 = buff%indx(1,iran)
            i1 = buff%indx(2,iran)
            nn = i1-i0+1

            ! Fill buffer
            buffer(ii+1:ii+nn) = real(tau(i0:i1))

            ! Advance index
            ii = ii + nn

          end do ! Ranges

        ! All
        else

          ! Fill buffer
          buffer = real(tau)

        end if

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if

      return

1000  umsg = 'Error opening tau file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing tau file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writetau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the height of optical depth equal to one of an inversion
      !! run in a file\n
      !!    filename(character(:)): Name of the file to write\n
      !!              tau(real(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writetau_inv(filename,tau,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: filename
      real, dimension(:), intent(in):: tau

      ! Local

      integer:: ierr,nn
      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: loffset


      ! Routine name
      urou = 'writetau_inv'

      ! Size
      nn = buff%buffer_size/4

      ! Open file
      call MPI_FILE_OPEN(MPI_COMM_SELF, &
                         trim(filename)//'/Tau', &
                         MPI_MODE_WRONLY, &
                         MPI_INFO_NULL,funit,ierr)
      if (ierr.ne.0) goto 1000

      !
      ! Column offset
      !

      ! Get offset
      loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                dble(buff%head_size)
      do while(loffset.gt.offlimit)
        offset = int(offlimit)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010
        loffset = loffset - offlimit
      end do
      offset = int(loffset)
      call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
      if (ierr.ne.0) goto 1010

      ! Write buffer
      call MPI_FILE_WRITE(funit,tau(1),nn, &
                          MPI_REAL,MPI_STATUS_IGNORE,ierr)
      if (ierr.ne.0) goto 1100

      ! Close file
      call MPI_FILE_CLOSE(funit, ierr)

      return

1000  umsg = 'Error opening tau file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing tau file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writetau_inv

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save the height where the optical depth is equal to one in
      !! RAM\n
      !!          e_tau(double(:)): RAM storage for height where the
      !!                            optical depth is equal to one\n
      !!            tau(double(:)): Height where the optical depth is
      !!                            equal to one\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine settau(e_tau,tau,buff)

      ! I/O

      type(IO_helper_class), intent(in):: buff
      double precision, dimension(:), intent(in):: tau
      real, dimension(:), intent(out):: e_tau

      ! Local

      integer:: iran,ii,i0,i1,nn


      ! Slaves leave
      if (pid.gt.0) return

      ! Initialize buffer
      ii = 0

      ! If specified
      if (buff%nran.gt.0) then

        ! For each entry to write
        do iran=1,buff%nran

          ! Range and size
          i0 = buff%indx(1,iran)
          i1 = buff%indx(2,iran)
          nn = i1-i0+1

          ! Fill buffer
          e_tau(ii+1:ii+nn) = real(tau(i0:i1))

          ! Advance index
          ii = ii + nn

        end do ! Ranges

      ! All
      else

        ! Fill buffer
        e_tau = real(tau)

      end if

      return

      end subroutine settau

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save inelastic collisional rates of active atoms in a file\n
      !!      Atom(Atom_class(:)): Structures with atomic data\n
      !!     folder(character(:)): Path to the output folder\n
      !!  btt(IO_helper_class(:)): Information on what needs to be
      !!                           stored of term-term collisions\n
      !!  bll(IO_helper_class(:)): Information on what needs to be
      !!                           stored of level-level collisions
      subroutine writecols(Atom,folder,btt,bll)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(IO_helper_class), intent(in):: btt,bll
      character(len=500), intent(in):: folder

      ! Local

      integer:: ierr,ii,iz,it1,it
      integer:: ia,iterm,iterm1,ilevel,ilevel1,ios
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Slaves
      if (pid.gt.0.or.nA.eq.0) then

        ! Control
        call control
        return

      end if ! Slaves

      ! Routine name
      urou = 'writecols'

      !
      ! 1D
      !
      if (run_mode.eq.0) then

        ! Open term-term file
        open (200,file=trim(folder)//'/cols-TT', status='unknown', &
              iostat=ios, err=1000, access='stream', action='write', &
              form='unformatted')

        ! Open level-level file
        open (300,file=trim(folder)//'/cols-LL', status='unknown', &
              iostat=ios, err=1001, access='stream', action='write', &
              form='unformatted')

        ! Identification
        write(200,err=1100) 'ct'
        write(300,err=1100) 'cl'

        ! Write number of atoms
        write(200,err=1100) NA
        write(300,err=1101) NA

        ! Write number of heights
        write(200,err=1100) NZ
        write(300,err=1101) NZ

        ! For each atom
        do ia=1,nA

          ! Write number of terms/levels
          write(200,err=1100) Atom(ia)%nmulti
          write(300,err=1101) Atom(ia)%nlevel

          ! For each term pair
          do iterm=1,Atom(ia)%nMulti
            do iterm1=1,Atom(ia)%nMulti

              ! Write rate
              write(200,err=1100) Atom(ia)%Ccoeff(iterm1,iterm,:)

            end do ! Destiny
          end do ! Origin

          ! For each level pair
          do ilevel=1,Atom(ia)%nlevel
            do ilevel1=1,Atom(ia)%nlevel

              ! Write rate
              write(300,err=1101) Atom(ia)%CcoeffJ(ilevel1,ilevel,:)

            end do ! Destiny
          end do ! Origin
        end do ! Atoms

        ! Close files
        close(200)
        close(300)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer with the largest size
        if (btt%buffer_size.gt.bll%buffer_size) then
          allocate(buffer(btt%buffer_size/4))
        else
          allocate(buffer(bll%buffer_size/4))
        end if

        ! Open file term-term
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/cols-TT', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(btt%buffer_size) + &
                  dble(btt%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (btt%nran.gt.0) then

          ! For each entry to write
          do ia=1,btt%nran

            ! For each height
            do iz=1,nz

              ! Advance
              ii = ii + 1

              ! Write a->b rate in buffer
              buffer(ii) = real(Atom(btt%indx(1,ia))% &
                                Ccoeff(btt%indx(2,ia), &
                                       btt%indx(3,ia),iz))

            end do ! Height

            ! For each height
            do iz=1,nz

              ! Advance
              ii = ii + 1

              ! Write b->a rate in buffer
              buffer(ii) = real(Atom(btt%indx(1,ia))% &
                                Ccoeff(btt%indx(3,ia), &
                                       btt%indx(2,ia),iz))

            end do ! Height
          end do ! Atom

        ! All
        else

          ! For each atom
          do ia=1,nA

            ! For each origin term
            do it=1,Atom(ia)%nMulti

              ! For each destiny term
              do it1=1,Atom(ia)%nMulti

                ! For each height
                do iz=1,nz

                  ! Advance
                  ii = ii + 1

                  ! Save
                  buffer(ii) = real(Atom(ia)%Ccoeff(it1,it,iz))

                end do ! Height
              end do ! Term
            end do ! Term
          end do ! Atom

        end if ! Specify ranges

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),btt%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Open file level-level
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/cols-LL', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1001

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(bll%buffer_size) + &
                  dble(bll%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1011
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1011

        ! Initialize buffer
        ii = 0

        ! If specified
        if (bll%nran.gt.0) then

          ! For each entry to write
          do ia=1,bll%nran

            ! For each height
            do iz=1,nz

              ! Advance
              ii = ii + 1

              ! Save a->b rate
              buffer(ii) = real(Atom(bll%indx(1,ia))% &
                                CcoeffJ(bll%indx(2,ia), &
                                        bll%indx(3,ia),iz))

            end do ! Height

            ! For each height
            do iz=1,nz

              ! Advance
              ii = ii + 1

              ! Save b->a rate
              buffer(ii) = real(Atom(bll%indx(1,ia))% &
                                CcoeffJ(bll%indx(3,ia), &
                                        bll%indx(2,ia),iz))
            end do ! Height
          end do ! Atom

        ! All
        else

          ! For each atom
          do ia=1,nA

            ! For each origin level
            do it=1,Atom(ia)%nlevel

              ! For each destiny level
              do it1=1,Atom(ia)%nlevel

                ! For each height
                do iz=1,nz

                  ! Advance
                  ii = ii + 1

                  ! Save rate
                  buffer(ii) = real(Atom(ia)%CcoeffJ(it1,it,iz))

                end do ! Height
              end do ! Destiny level
            end do ! Origin level
          end do ! Atom

        end if ! Specify range

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),bll%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1301

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run type

      ! Control
      call control
      return

1000  umsg = 'Error opening collisions TT file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking collisions TT file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing collisions TT file'
      close(200)
      close(300)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing collisions TT file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1001  umsg = 'Error opening collisions LL file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1011  umsg = 'Error seeking collisions LL file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1101  umsg = 'Error writing collisions LL file'
      close(200)
      close(300)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1301  umsg = 'Error writing collisions LL file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writecols

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save damping parameter of active atoms in a file\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      folder(character(:)): Path to the output folder\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writedamp(Atom,Atmo,folder,buff)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder

      ! Local

      integer:: ierr,ii
      integer:: iran,ia,iterml,itermu,itran,ios
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset
      double precision, dimension(nz):: Dw,DwT


      ! Slaves
      if (pid.gt.0.or.nA.eq.0) then

        ! Control
        call control
        return

      end if ! Slaves

      ! Routine name
      urou = 'writedamp'


      !
      ! 1D
      !
      if (run_mode.eq.0) then

        ! Open file to write into
        open (200,file=trim(folder)//'/damping', status='unknown', &
              iostat=ios, err=1000, access='stream', action='write', &
              form='unformatted')

        ! Identification
        write(200,err=1100) 'da'

        ! Write number of atoms
        write(200,err=1100) NA

        ! Write number of heights
        write(200,err=1100) NZ

        ! For each atom
        do ia=1,nA

          ! Write number of terms/levels and transitions
          write(200,err=1100) Atom(ia)%ntran

          ! Thermal width
          DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

          ! For each transition
          do itran=1,Atom(ia)%ntran

            ! Doppler width
            Dw = Atom(ia)%Dfreq(itran)*sqrt(DwT*DwT + &
                                            Atmo%vmi*Atmo%vmi)

            ! Inverrse
            Dw = 1d0/Dw

            ! Identify involved terms
            itermu = -1

            ! Get terms
            itermu = Atom(ia)%fst(itran)%itermu
            iterml = Atom(ia)%fst(itran)%iterml

            ! Write total damping
            write(200,err=1100) (Atom(ia)%ldamp(itran,:) + &
                                 Atom(ia)%damp(itermu,:) + &
                                 Atom(ia)%damp(iterml,:))*Dw

          end do ! Transitions
        end do ! Atoms

        ! Close file
        close(200)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/damping', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            ia = buff%indx(1,iran)
            itran = buff%indx(2,iran)

            ! Thermal width
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

            ! Doppler width and inverse
            Dw = Atom(ia)%Dfreq(itran)* &
                 sqrt(DwT*DwT + Atmo%vmi*Atmo%vmi)
            Dw = 1d0/Dw

            ! Get terms
            itermu = Atom(ia)%fst(itran)%itermu
            iterml = Atom(ia)%fst(itran)%iterml

            ! Save
            buffer(ii+1:ii+nz) = real((Atom(ia)%ldamp(itran,:) + &
                                       Atom(ia)%damp(itermu,:) + &
                                       Atom(ia)%damp(iterml,:))*Dw)

            ! Advance buffer
            ii = ii + nz

          end do ! Ranges

        ! All
        else

          ! For each atom
          do ia=1,nA

            ! Thermal width
            DwT = Atom(ia)%cDopp*sqrt(Atmo%T)

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Doppler width
              Dw = Atom(ia)%Dfreq(itran)*sqrt(DwT*DwT + &
                                              Atmo%vmi*Atmo%vmi)

              ! Inverse
              Dw = 1d0/Dw

              ! Get terms
              itermu = Atom(ia)%fst(itran)%itermu
              iterml = Atom(ia)%fst(itran)%iterml

              ! Save
              buffer(ii+1:ii+nz) = real((Atom(ia)%ldamp(itran,:) + &
                                         Atom(ia)%damp(itermu,:) + &
                                         Atom(ia)%damp(iterml,:))*Dw)

              ! Advance buffer
              ii = ii + nz

            end do ! Transitions
          end do ! Atoms

        end if ! Specified transitions

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

    end if ! Run type

      ! Control
      call control

      return

1000  umsg = 'Error opening damping file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking damping file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing damping file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing damping file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writedamp

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save elastic collisional rates of active atoms in a file\n
      !!       Atom(Atom_class(:)): Structures with atomic data\n
      !!      folder(character(:)): Path to the output folder\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writeqel(Atom,folder,buff)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder

      ! Local

      integer:: ierr,ii
      integer:: iran,ia,iterml,itermu,itran,ios
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset


      ! Slaves
      if (pid.gt.0.or.nA.eq.0) then

        ! Control
        call control
        return

      end if ! Slaves

      ! Routine name
      urou = 'writeqel'


      !
      ! 1D
      !
      if (run_mode.eq.0) then

        ! Open file to write into
        open (200,file=trim(folder)//'/qel', status='unknown', &
              iostat=ios, err=1000, access='stream', action='write', &
              form='unformatted')

        ! Identification
        write(200,err=1100) 'qe'

        ! Write number of atoms
        write(200,err=1100) NA

        ! Write number of heights
        write(200,err=1100) NZ

        ! For each atom
        do ia=1,nA

          ! Write number of terms/levels and transitions
          write(200,err=1100) Atom(ia)%ntran

          ! For each transition
          do itran=1,Atom(ia)%ntran

            ! Get terms
            itermu = Atom(ia)%fst(itran)%itermu
            iterml = Atom(ia)%fst(itran)%iterml

            ! Write rate
            write(200,err=1100) Atom(ia)%qel(itran,:)

          end do ! Transitions

        end do ! Atoms

        ! Close files
        close(200)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/qel', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize buffer
        ii = 0

        ! If specified
        if (buff%nran.gt.0) then

          ! For each entry to write
          do iran=1,buff%nran

            ! Atom and transition
            ia = buff%indx(1,iran)
            itran = buff%indx(2,iran)

            ! Get terms
            itermu = Atom(ia)%fst(itran)%itermu
            iterml = Atom(ia)%fst(itran)%iterml

            ! Save
            buffer(ii+1:ii+nz) = real(Atom(ia)%qel(itran,:))

            ! Advance buffer
            ii = ii + nz

          end do ! Ranges

        ! All
        else

          ! For each atom
          do ia=1,nA

            ! For each transition
            do itran=1,Atom(ia)%ntran

              ! Get terms
              itermu = Atom(ia)%fst(itran)%itermu
              iterml = Atom(ia)%fst(itran)%iterml

              ! Save
              buffer(ii+1:ii+nz) = real(Atom(ia)%qel(itran,:))

              ! Advance buffer
              ii = ii + nz

            end do ! Transitions
          end do ! Atoms

        end if ! Specified ranges

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      ! Control
      call control

      return

1000  umsg = 'Error opening qel file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking qel file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing qel file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing qel file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writeqel

!#####################################################################
!#####################################################################
!#####################################################################

      !> Save background opacity, scattering coefficient, and
      !! emissivity in a file\n
      !!     Cont(Continuum_class): Structure with background opacity
      !!                            data\n
      !!          omega(double(:)): Frequency array\n
      !!      folder(character(:)): Path to the output folder\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writeback(Cont,omega,folder,MPID,buff)

      ! I/O

      type(Continuum_class), intent(inout):: Cont
      type(MPI_class), intent(in):: MPID
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder
      double precision, dimension(:), intent(in):: omega

      ! Local

      integer:: ierr,i0,i1,nn,ii,jj,iran
      integer:: iz,nfl,ndl,nxdir,idir
      integer:: iproc,iproc0,bsize,psize,ios
      integer, dimension(:), allocatable:: ndir
      integer(kind=MPI_OFFSET_KIND):: offset

      real, dimension(:), allocatable:: buffer

      double precision:: loffset
      double precision, dimension(:), allocatable:: dbuffer


      ! Routine name
      urou = 'writeback'

      ! Master opens file and writes dimensions
      if (pid.eq.0) then

        !
        ! 1D
        !
        if (run_mode.eq.0) then

          ! Open file
          open (200,file=trim(folder)//'/background', &
                status='unknown', iostat=ios, err=1000, &
                access='stream', action='write', form='unformatted')

          ! Identification
          write(200,err=1100) 'ba'

          ! Write number of frequencies
          write(200,err=1100) nfreq

          ! Write frequencies
          write(200,err=1100) omega

          ! Write number of heights
          write(200,err=1100) NZ

        end if ! 1D mode
      end if ! Master

      ! If MPI
      if (nproc.gt.1) then

        ! Control
        call control

        ! Allocate number of directions for every cpu
        allocate(ndir(0:nproc-1))

        ! Master initialize
        if (pid.eq.0) Cont%ndir = 0

        ! Master gather number of directions
        call MPI_GATHER(Cont%ndir, 1, MPI_INTEGER, ndir(0), 1, &
                        MPI_INTEGER, 0, MPI_COMM_RT, ierr)

        ! Master
        if (pid.eq.0) then

          ! Look for the maximum number of directions
          nxdir = maxval(ndir)

          ! Last processor checked
          iproc0 = 1

          ! Maximum dimension to receive
          bsize = MPID%nxfreq*nz*3*nxdir

          ! Allocate buffer
          allocate(dbuffer(bsize))

          ! Allocate background
          allocate(Cont%c(nfreq,3,nxdir,nz))
          Cont%c = 0d0

          ! For each process
          do iproc=1,nproc-1

            ! Size of this process
            nfl = MPID%nf(iproc)
            ndl = ndir(iproc)
            psize = nz*nfl*ndl*3

            ! Receive buffer
            do while (.True.)
              call MPI_RECV(dbuffer, psize, MPI_DOUBLE_PRECISION, &
                            iproc, iproc, MPI_COMM_RT, &
                            MPI_STATUS_IGNORE, ierr)
              if (ierr.eq.0) exit
            end do

            ! If the processor had the maximum number of
            ! directions
            if (ndl.eq.nxdir) then

              Cont%c(MPID%if0(iproc):MPID%if1(iproc),:,:,:) = &
                       reshape(dbuffer(1:psize), (/ nfl,3,ndl,nz /))

            ! If the buffer is smaller than the maximum
            else

              ! We have to replicate the directions, so for each
              ! height
              do iz=1,nz

                ! For each direction
                do idir=1,nxdir

                  ! Replicate same buffer in the relevant range
                  Cont%c(MPID%if0(iproc):MPID%if1(iproc),:,idir,iz)= &
                        reshape(dbuffer((iz-1)*(nfl*3)+1:iz*nfl*3), &
                                                        (/ nfl,3 /))

                end do ! directions
              end do ! heights

            end if ! Maximum number of directions

          end do ! Processors

          ! Free
          deallocate(dbuffer)

        ! Slaves
        else

          ! Compute size of data to send
          psize = nz*MPID%nf(pid)*3*Cont%ndir

          ! Send data
          do while (.True.)
            call MPI_SEND(Cont%c(MPID%if0(pid),1,1,1), &
                          psize, MPI_DOUBLE_PRECISION, 0, pid, &
                          MPI_COMM_RT, ierr)
            if (ierr.eq.0) exit
          end do

        end if ! Master or slave

        ! Free
        deallocate(ndir)

      ! If serial
      else

        ! Save nxdir
        nxdir = Cont%ndir

      end if ! MPI

      ! If Master
      if (pid.eq.0) then

        !
        ! 1D
        !
        if (run_mode.eq.0) then

          ! Write number of directions
          write(200,err=1100) nxdir

          ! Write the full block
          write(200,err=1100) Cont%c

          ! If MPI, deallocate background
          if (MPID%mpi) deallocate(Cont%c)

          ! Close files
          close(200)

        !
        ! 1.5D
        !
        else if (run_mode.eq.1) then

          ! Allocate buffer
          allocate(buffer(buff%buffer_size/4))

          ! Open file
          call MPI_FILE_OPEN(MPI_COMM_SELF, &
                             trim(folder)//'/background', &
                             MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                             ierr)
          if (ierr.ne.0) goto 1000

          !
          ! Column offset
          !

          ! Get offset
          loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                    dble(buff%head_size)
          do while(loffset.gt.offlimit)
            offset = int(offlimit)
            call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
            if (ierr.ne.0) goto 1010
            loffset = loffset - offlimit
          end do
          offset = int(loffset)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010

          ! Initialize buffer
          ii = 0

          ! If specified
          if (buff%nran.gt.0) then

            ! For each height
            do iz=1,nz

              ! For each variable
              do jj=1,3

                ! For each entry to write
                do iran=1,buff%nran

                  ! Range and size
                  i0 = buff%indx(1,iran)
                  i1 = buff%indx(2,iran)
                  nn = i1-i0+1

                  ! Fill buffer
                  buffer(ii+1:ii+nn) = real(Cont%c(i0:i1,jj,1,iz))

                  ! Advance index
                  ii = ii + nn

                end do ! Ranges
              end do ! Variables
            end do ! Heights

          ! All
          else

            ! Size
            nn = nfreq*nz*3

            ! Fill buffer
            buffer = real(reshape(Cont%c(:,:,1,:), (/ nn /)))

          end if

          ! Write buffer
          call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                              MPI_REAL,MPI_STATUS_IGNORE,ierr)
          if (ierr.ne.0) goto 1300

          ! Close file
          call MPI_FILE_CLOSE(funit, ierr)

          ! Free
          deallocate(buffer)

        end if ! Run mode
      end if ! Master

      ! Control
      call control

      return

1000  umsg = 'Error opening background file'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking background file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing background file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing background file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writeback

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write the model atmosphere in a file. The file is in ASCII
      !! for a 1D synthesis and in binary for a 1.5D synthesis, with
      !! the same format than the corresponding input for the latter\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!      Bfield(Bfield_class): Structure with magnetic field
      !!                            data\n
      !!      folder(character(:)): Path to the output folder\n
      !!  buff(IO_helper_class(:)): Information on what needs to be
      !!                            stored of the relevant variable
      subroutine writeatmo(Atmo,Bfield,folder,buff)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(IO_helper_class), intent(in):: buff
      character(len=500), intent(in):: folder

      ! Local

      integer:: iz, ios
      integer:: ierr,ii
      integer(kind=MPI_OFFSET_KIND):: offset

      double precision:: loffset
      double precision, dimension(:), allocatable:: buffer


      ! Routine name
      urou = 'writeatmo'

      ! Slaves
      if (pid.gt.0) then

        ! Control
        call control
        return

      end if ! Slaves

      !
      ! 1D or inversion
      !
      if (run_mode.le.0) then

        ! Master opens file and writes dimensions
        open (200,file=trim(folder)//'/atmos.dat', &
              status='unknown', iostat=ios, err=1000, action='write')

        ! If tau scale
        if (ztau) then

          ! Write header
          write(200,'(A)',err=1100) '! Atmospheric data'
          write(200,'(A)',err=1100) '!             Height [km]'// &
                                     '                tau cont'//&
                                     '        Chi cont [cm^-1]'//&
                                     '         Temperature [K]'//&
                                     ' Gas pressure [dyn/cm^2]'//&
                                     '        Density [g/cm^3]'//&
                                     '   Magnetic strength [G]'//&
                                     '    Magnetic inclination'//&
                                     '        Magnetic azimuth'//&
                                     '       Velocity x [km/s]'//&
                                     '       Velocity y [km/s]'//&
                                     '       Velocity z [km/s]'//&
                                     '  Microturbulence [km/s]'//&
                                     ' e^- pressure [dyn/cm^2]'//&
                                     '     e^- density [cm^-3]'//&
                                     '       H density [cm^-3]'//&
                                     '  atom H density [cm^-3]'//&
                                     '      H- density [cm^-3]'//&
                                     '    HI_0 density [cm^-3]'//&
                                     '    HI_1 density [cm^-3]'//&
                                     '    HI_2 density [cm^-3]'//&
                                     '    HI_3 density [cm^-3]'//&
                                     '    HI_4 density [cm^-3]'//&
                                     '     HII density [cm^-3]'

          ! For each height, output
          do iz=1,nz

            ! Write data
            write(200,'(1x,24(2x,es22.15))') &
                                         Atmo%zalt(iz), &
                                         Atmo%z(iz), &
                                         Atmo%chi500(iz), &
                                         Atmo%T(iz), &
                                         Atmo%Pg(iz), &
                                         Atmo%rho(iz), &
                                         Bfield%Bstrength(iz), &
                                         Bfield%Btheta(iz), &
                                         Bfield%Bphi(iz), &
                                         Atmo%vx(iz)*1d6*c, &
                                         Atmo%vy(iz)*1d6*c, &
                                         Atmo%vz(iz)*1d6*c, &
                                         Atmo%vmi(iz)*1d6*c, &
                                         Atmo%Pe(iz), &
                                         Atmo%ne(iz), &
                                         Atmo%nHT(iz), &
                                         Atmo%nHa(iz), &
                                         Atmo%nHm(iz), &
                                         Atmo%nh(iz,1), &
                                         Atmo%nh(iz,2), &
                                         Atmo%nh(iz,3), &
                                         Atmo%nh(iz,4), &
                                         Atmo%nh(iz,5), &
                                         Atmo%nh(iz,6)
          end do ! Height nodes

        ! If height scale
        else

          ! Write header
          write(200,'(A)',err=1100) '! Atmospheric data'
          write(200,'(A)',err=1100) '!             Height [km]'// &
                                     '                tau cont'//&
                                     '        Chi cont [cm^-1]'//&
                                     '         Temperature [K]'//&
                                     ' Gas pressure [dyn/cm^2]'//&
                                     '        Density [g/cm^3]'//&
                                     '   Magnetic strength [G]'//&
                                     '    Magnetic inclination'//&
                                     '        Magnetic azimuth'//&
                                     '       Velocity x [km/s]'//&
                                     '       Velocity y [km/s]'//&
                                     '       Velocity z [km/s]'//&
                                     '  Microturbulence [km/s]'//&
                                     ' e^- pressure [dyn/cm^2]'//&
                                     '     e^- density [cm^-3]'//&
                                     '       H density [cm^-3]'//&
                                     '  atom H density [cm^-3]'//&
                                     '      H- density [cm^-3]'//&
                                     '    HI_0 density [cm^-3]'//&
                                     '    HI_1 density [cm^-3]'//&
                                     '    HI_2 density [cm^-3]'//&
                                     '    HI_3 density [cm^-3]'//&
                                     '    HI_4 density [cm^-3]'//&
                                     '     HII density [cm^-3]'

          ! For each height
          do iz=1,nz

            ! Write data
            write(200,'(1x,24(2x,es22.15))') &
                                         Atmo%z(iz)*1d-5, &
                                         Atmo%zalt(iz), &
                                         Atmo%chi500(iz), &
                                         Atmo%T(iz), &
                                         Atmo%Pg(iz), &
                                         Atmo%rho(iz), &
                                         Bfield%Bstrength(iz), &
                                         Bfield%Btheta(iz), &
                                         Bfield%Bphi(iz), &
                                         Atmo%vx(iz)*1d6*c, &
                                         Atmo%vy(iz)*1d6*c, &
                                         Atmo%vz(iz)*1d6*c, &
                                         Atmo%vmi(iz)*1d6*c, &
                                         Atmo%Pe(iz), &
                                         Atmo%ne(iz), &
                                         Atmo%nHT(iz), &
                                         Atmo%nHa(iz), &
                                         Atmo%nHm(iz), &
                                         Atmo%nh(iz,1), &
                                         Atmo%nh(iz,2), &
                                         Atmo%nh(iz,3), &
                                         Atmo%nh(iz,4), &
                                         Atmo%nh(iz,5), &
                                         Atmo%nh(iz,6)
          end do ! Heights

        end if ! Height or tau scale

        ! Close files
        close(200)

      !
      ! 1.5D
      !
      else if (run_mode.eq.1) then

        ! Allocate buffer
        allocate(buffer(buff%buffer_size/4))

        ! Open file
        call MPI_FILE_OPEN(MPI_COMM_SELF, &
                           trim(folder)//'/atmo.hrt', &
                           MPI_MODE_WRONLY, MPI_INFO_NULL, funit, &
                           ierr)
        if (ierr.ne.0) goto 1000

        !
        ! Column offset
        !

        ! Get offset
        loffset = dble(icoords(3)-1)*dble(buff%buffer_size) + &
                  dble(buff%head_size)
        do while(loffset.gt.offlimit)
          offset = int(offlimit)
          call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
          if (ierr.ne.0) goto 1010
          loffset = loffset - offlimit
        end do
        offset = int(loffset)
        call MPI_FILE_SEEK(funit, offset, MPI_SEEK_CUR, ierr)
        if (ierr.ne.0) goto 1010

        ! Initialize
        ii = 0

        ! If tau scale
        if (ztau) then

          buffer(ii+1:ii+nz) = Atmo%zalt
          ii = ii + nz
          buffer(ii+1:ii+nz) = Atmo%z
          ii = ii + nz

        ! Height scale
        else

          buffer(ii+1:ii+nz) = Atmo%z*1d-5
          ii = ii + nz
          buffer(ii+1:ii+nz) = Atmo%zalt
          ii = ii + nz

        end if

        !
        ! Fill rest of buffer
        !

        ! Chi500
        buffer(ii+1:ii+nz) = Atmo%chi500
        ii = ii + nz
        ! T
        buffer(ii+1:ii+nz) = Atmo%T
        ii = ii + nz
        ! Pg
        buffer(ii+1:ii+nz) = Atmo%Pg
        ii = ii + nz
        ! rho
        buffer(ii+1:ii+nz) = Atmo%rho
        ii = ii + nz
        ! Bx
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             sin(Bfield%Btheta)* &
                             cos(Bfield%Bphi)
        ii = ii + nz
        ! By
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             sin(Bfield%Btheta)* &
                             sin(Bfield%Bphi)
        ii = ii + nz
        ! Bz
        buffer(ii+1:ii+nz) = Bfield%bstrength* &
                             cos(Bfield%Btheta)
        ii = ii + nz
        ! vx
        buffer(ii+1:ii+nz) = Atmo%vx*1d6*c
        ii = ii + nz
        ! vy
        buffer(ii+1:ii+nz) = Atmo%vy*1d6*c
        ii = ii + nz
        ! vz
        buffer(ii+1:ii+nz) = Atmo%vz*1d6*c
        ii = ii + nz
        ! vmi
        buffer(ii+1:ii+nz) = Atmo%vmi*1d6*c
        ii = ii + nz
        ! Pe
        buffer(ii+1:ii+nz) = Atmo%Pe
        ii = ii + nz
        ! ne
        buffer(ii+1:ii+nz) = Atmo%ne
        ii = ii + nz
        ! nHT
        buffer(ii+1:ii+nz) = Atmo%nHT
        ii = ii + nz
        ! nHa
        buffer(ii+1:ii+nz) = Atmo%nHa
        ii = ii + nz
        ! nHm
        buffer(ii+1:ii+nz) = Atmo%nHm
        ii = ii + nz
        ! nH
        buffer(ii+1:ii+nz*6) = reshape(Atmo%nH, (/ nz*6 /))

        ! Write buffer
        call MPI_FILE_WRITE(funit,buffer(1),buff%buffer_size/4, &
                            MPI_REAL,MPI_STATUS_IGNORE,ierr)
        if (ierr.ne.0) goto 1300

        ! Close file
        call MPI_FILE_CLOSE(funit, ierr)

        ! Free
        deallocate(buffer)

      end if ! Run mode

      ! Control
      call control

      return

1000  umsg = 'Error opening atmospheric file to write'
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1010  umsg = 'Error seeking atmospheric file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1100  umsg = 'Error writing atmospheric file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return
1300  umsg = 'Error writing atmospheric file'
      call MPI_FILE_CLOSE(funit, ierr)
      call abortedS(umsg,urou,.True.,.True.)
      call control
      return

      end subroutine writeatmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Write a 1D model atmosphere in a file in the 1D synthesis
      !! format\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!       outtypo(integer): What variables to write in the
      !!                         output\n
      !!   folder(character(:)): Path to the output folder\n
      !! filename(character(:)): Path to the original atmospheric
      !!                         file
      subroutine wAtmo(Atmo,outtypo,folder,filename)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      character(len=500), intent(in):: folder,filename
      integer, intent(in):: outtypo

      ! Local

      integer:: iz,ios,typo

      double precision:: zfac,vfac


      ! Routine name
      urou = 'wAtmo'

      ! Slaves
      if (pid.gt.0) then

        ! Control
        call control
        return

      end if ! Slaves

      ! If type of output larger than zero
      if (outtypo.gt.0) then

        ! Shifting one gives us the desired type
        typo = outtypo - 1

      ! Otherwise
      else

        ! Get from the model itself
        typo = Atmo%typo

      end if ! Type of output

      ! Open file
      open (200,file=trim(folder)//'/atmo.atmos', &
            status='unknown', iostat=ios, err=1000, action='write')

      ! Header
      write(200,'(A)',err=1100) '* Written by HANLERT as an '// &
                                'update of '//trim(filename)
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '  HANLERT-model'

      ! If tau scale
      if (ztau) then

        ! Get conversion factor and write
        zfac = 1d0
        write(200,'(A,1x,es23.16)',err=1100) &
          'TAU SCALE',1d2/Atmo%tfreq

      ! If height scale
      else

        ! Get conversion factor and write
        zfac = 1d-5
        write(200,'(A,1x,es23.16)',err=1100) &
          'HEIGHT SCALE',1d2/Atmo%tfreq

      end if ! Type of vertical scale

      ! Log g from original
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '* LG G'
      write(200,'(2x,es23.16)',err=1100) Atmo%logg

      ! Height number of nodes from original
      write(200,'(A)',err=1100) '*'
      write(200,'(A)',err=1100) '* NDEP'
      write(200,'(i6)',err=1100) Atmo%nz

      ! Velocity factor to recover km/s
      vfac = c*1d6

      ! First block
      write(200,'(A)',err=1100) '*'

      ! If Standard
      if (typo.eq.0) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'
        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Tau or height

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Hydrogen populations
        write(200,'(A)',err=1100) '*'
        write(200,'(A)',err=1100) '* HYDROGEN POPULATIONS (cm^-3)'
        write(200,'(A)',err=1100) '*     NH(1)      '// &
                                  '      NH(2)      '// &
                                  '      NH(3)      '// &
                                  '      NH(4)      '// &
                                  '      NH(5)      '// &
                                  '       NP        '

        ! For every height
        do iz=1,nZ

          ! Write hydrogen
          write(200,'(6(1x,es17.10))',err=1100) Atmo%nh(iz,:)

        end do ! Heights

        ! If there are helium populations
        if (Atmo%nhe(1,1).gt.-0.1) then

          ! Write header
          write(200,'(A)',err=1100) '*'
          write(200,'(A)',err=1100) '* HELIUM POPULATIONS (cm^-3)'

          ! For every height
          do iz=1,nZ

            ! Write populations
            write(200,'(4(1x,es17.10))',err=1100) Atmo%nhe(iz,:)

          end do ! Heights

        end if ! Helium populations present


      ! If electron number density
      else if (typo.eq.1) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)        Ne (cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Type of height scale

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Write tag
        write(200,'(A)',err=1100) 'ne'

      ! If electron pressure
      else if (typo.eq.2) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)     Pe (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)     Pe (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Type of height scale

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%Pe(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Write tag
        write(200,'(A)',err=1100) 'pe'

      ! If electron mass density
      else if (typo.eq.3) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)  dens_e (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)  dens_e (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Type of height scale

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%ne(iz)*me*1d3,Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Write tag
        write(200,'(A)',err=1100) 'rhoe'

      ! If gas pressure
      else if (typo.eq.4) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)     Pg (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)     Pg (dyn/cm^2)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Type of theight scale

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%Pg(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Write tag
        write(200,'(A)',err=1100) 'pg'

      ! If gas mass density
      else if (typo.eq.5) then

        ! Tau
        if (ztau) then

          ! Write
          write(200,'(A)',err=1100) '*    Optical depth'// &
                                    '         '// &
                                    'TEMP (K)    dens (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        ! Height
        else

          ! Write
          write(200,'(A)',err=1100) '*      HEIGHT (km)'// &
                                    '         '// &
                                    'TEMP (K)    dens (g*cm^-3)'// &
                                    '         vz (km/s)    '// &
                                    'Vmicro (km/s)'// &
                                    '         vx (km/s)'// &
                                    '         vy (km/s)'

        end if ! Type of height data

        ! For every height
        do iz=1,nZ

          ! Write data
          write(200,'(7(1x,es17.10))',err=1100) Atmo%z(iz)*zfac, &
            Atmo%T(iz),Atmo%rho(iz),Atmo%vz(iz)*vfac, &
            Atmo%vmi(iz)*vfac,Atmo%vx(iz)*vfac,Atmo%vy(iz)*vfac

        end do ! Heights

        ! Write tag
        write(200,'(A)',err=1100) 'rho'

      end if ! Type of output

      ! Close files
      close(200)

      ! Control
      call control

      return

1000  umsg = 'Error opening atmospheric model file to write'
      call abortedS(umsg,urou,.True.,.True.)
      call control
1100  umsg = 'Error writing atmospheric model file'
      close(200)
      call abortedS(umsg,urou,.True.,.True.)
      call control

      end subroutine wAtmo

!#####################################################################
!#####################################################################
!
! The routines beyond this point have debugging purposes and are not
! compiled unless specific flags are activated in the compilation.
! Therefore, they are not documented in the header and do not fully
! follow the source standard.
!
!#####################################################################
!#####################################################################

#ifdef DEBUGRHO00
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump rho00 solution into a file\n
      !!     Atom(Atom_class(:)): Structures with atomic data\n
      !!  folder(character(500)): Output folder\n
      !!           iter(integer): Iteration index
      subroutine dump_rho00(Atom,folder,iter)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      character(len=500), intent(in):: folder
      integer, intent(in):: iter

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: ia,iz,it,iJ
      double precision:: ff


      ! Sanity
      if (nA.lt.1) return

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_rho00'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_rho00_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.-2.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Write
      if (iter.eq.-2) then
        write(800,*) 'Initial rho00'
      else if (iter.eq.-1) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'Final rho00'
      else
        write(800,*) ''
        write(800,*) ''
        write(800,'("Iteration:",i3)') iter
      end if
      do ia=1,na
        write(800,'("  o Atom ",A)') Atom(ia)%element
        do it=1,Atom(ia)%nMulti
        do iJ=1,Atom(ia)%nJ(it)
          do iz=Rz0,Rz1
            ff = sqrt(2d0*Atom(ia)%rJval(iJ,it)+1d0)*Atom(ia)%n(iz)
            write(800,'(i4,1x,i2,1x,i3,3(1x,es15.8))') &
                  it,iJ,iz, &
                  dble(Atom(ia)%crho(Atom(ia)%irho(it)% &
                                           Jrho(iJ,iJ)%kq(0,0),iz)), &
               ff*dble(Atom(ia)%crho(Atom(ia)%irho(it)% &
                                           Jrho(iJ,iJ)%kq(0,0),iz)), &
                                           Atom(ia)%n(iz)
          end do
        end do
        end do
      end do

      ! Close
      close(800)

      end subroutine dump_rho00

#endif
#ifdef DEBUGJ00
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump integrated J00 into a file\n
      !!     Atom(Atom_class(:)): Structures with the atomic data\n
      !!    J00(double(:,:,:,:)): Mean intensity integrated over
      !!                          absorption profile\n
      !!   J00S(double(:,:,:,:)): Mean intensity integrated over
      !!                          emission profile\n
      !!     J00P(double(:,:,:)): Intensity integrals in the
      !!                          photoionization rates
      !!  folder(character(500)): Output folder\n
      !!           iter(integer): Iteration index
      subroutine dump_j00(Atom,J00,J00S,J00P,folder,iter)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      character(len=500), intent(in):: folder
      integer, intent(in):: iter
      double precision, dimension(nxt,Rz0:Rz1), intent(in):: J00, J00S
      double precision, dimension(nxphot,2,Rz0:Rz1), intent(in):: J00P

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: ia,iz,itran,jtran
      double precision, parameter:: ff=1d0/299792458d5


      ! Sanity
      if (nA.lt.1) return

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_j00'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_j00_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.-2.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Write
      if (iter.eq.-2) then
        write(800,*) 'Initial J00'
      else if (iter.eq.-1) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'Final J00'
      else
        write(800,*) ''
        write(800,*) ''
        write(800,'("Iteration:",i3)') iter
      end if
      do ia=1,na
        write(800,'("  o Atom ",A)') Atom(ia)%element
        write(800,'("    ",A)') 'J00'
        do itran=1,Atom(ia)%nftran
          jtran = itran + Atom(ia)%tfshift
          do iz=Rz0,Rz1
            write(800,'(i4,1x,i3,2(1x,es15.8))') &
                  itran,iz,J00(jtran,iz)*ff,J00S(jtran,iz)*ff
          end do
        end do
        write(800,'("    ",A)') 'J00P'
        do itran=1,Atom(ia)%nphot
          jtran = itran + Atom(ia)%pshift
          do iz=Rz0,Rz1
            write(800,'(i4,1x,i3,3(1x,es15.8))') &
                  itran,iz,&
                  J00P(jtran,1,iz)*1d8, &
                  J00P(jtran,2,iz)*1d8, &
                  Atom(ia)%phot(jtran)%TEI(iz)*1d8
          end do
        end do
      end do

      ! Close
      close(800)

      end subroutine dump_j00

#endif
#ifdef DEBUGLAMBDA
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump lambda operators into a file\n
      !!     Atom(Atom_class(:)): Structures with the atomic data\n
      !!       LamL(double(:,:)): Lambda operator for bound-bound
      !!                          transitions\n
      !!     LamP(double(:,:,:)): Lambda operator for bound-free
      !!                          transitions\n
      !!  folder(character(500)): Output folder\n
      !!           iter(integer): Iteration index
      subroutine dump_lambda(Atom,LamL,LamP,folder,iter)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      character(len=500), intent(in):: folder
      integer, intent(in):: iter
      double precision, dimension(nxb,nxt,Rz0:Rz1), intent(in):: LamL
      double precision, dimension(nxb,nxphot,2,Rz0:Rz1), &
                                                     intent(in):: LamP

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: ia,iz,itran,jtran
      double precision, parameter:: ff=1d0/299792458d5


      ! Sanity
      if (nA.lt.1) return

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_lambda'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_lambda_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.1.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Write
      if (iter.gt.1) then
        write(800,*) ''
        write(800,*) ''
      end if
      write(800,'("Iteration:",i3)') iter
      do ia=1,na
        write(800,'("  o Atom ",A)') Atom(ia)%element
        write(800,'("    ",A)') 'b-b'
        do itran=1,Atom(ia)%nftran
          jtran = itran + Atom(ia)%tfshift
          do iz=Rz0,Rz1
            write(800,'(i4,1x,i3,2(1x,es15.8))') &
                  itran,iz,LamL(1,jtran,iz)*ff
          end do
        end do
        write(800,'("    ",A)') 'b-f'
        do itran=1,Atom(ia)%nphot
          jtran = itran + Atom(ia)%pshift
          do iz=Rz0,Rz1
            write(800,'(i4,1x,i3,2(1x,es15.8))') &
                  itran,iz,LamP(1,jtran,1,iz)*ff,LamP(1,jtran,2,iz)*ff
          end do
        end do
      end do

      ! Close
      close(800)

      end subroutine dump_lambda

#endif
#ifdef DEBUGRHOKQ
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump rho solution into a file\n
      !!    Atom(Atom_class(:)): Structures with the atomic data\n
      !! folder(character(500)): Output folder\n
      !!          iter(integer): Iteration index
      subroutine dump_rho(Atom,folder,iter)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      character(len=500), intent(in):: folder
      integer, intent(in):: iter

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: ia,iz,it,iJ,iJ1,K,iQ
      double precision:: rJ,rJ1


      ! Sanity
      if (nA.lt.1) return

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_rho'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_rho_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.-3.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Write
      if (iter.eq.-3) then
        write(800,*) 'Initial rhoKQ'
      else if (iter.eq.-2) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'After correction rhoKQ'
      else if (iter.eq.-1) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'Final rhoKQ'
      else
        write(800,*) ''
        write(800,*) ''
        write(800,'("Iteration:",i3)') iter
      end if
      do ia=1,na
        write(800,'("  o Atom ",A)') Atom(ia)%element
        do it=1,Atom(ia)%nMulti
        do iJ=1,Atom(ia)%nJ(it)
          rJ = Atom(ia)%rJval(iJ,it)
          do iJ1=1,Atom(ia)%nJ(it)
            rJ1 = Atom(ia)%rJval(iJ1,it)
            do K=nint(abs(rJ-rJ1)), &
                 min(nint(rJ+rJ1),Atom(ia)%Kcut(it))
            do iQ=0,K
              if (iQ == 0) then
                do iz=Rz0,Rz1
                  write(800,'(i4,4(1x,i2),1x,i3,1x,es15.8)') &
                        it,iJ,iJ1,K,iQ,iz, &
                        dble(Atom(ia)%crho(Atom(ia)%irho(it)% &
                                           Jrho(iJ,iJ1)%kq(iQ,K),iz))
                end do
              else
                do iz=1,nZ
                  write(800,'(i4,4(1x,i2),1x,i3,2(1x,es15.8))') &
                      it,iJ,iJ1,K,iQ,iz, &
                      dble(Atom(ia)%crho(Atom(ia)%irho(it)% &
                                         Jrho(iJ,iJ1)%kq(iQ,K),iz)), &
                      dimag(Atom(ia)%crho(Atom(ia)%irho(it)% &
                                         Jrho(iJ,iJ1)%kq(iQ,K),iz))
                end do
              end if
            end do
            end do
          end do
        end do
        end do
      end do

      ! Close
      close(800)

      end subroutine dump_rho

#endif
#ifdef DEBUGJKQ
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump rho solution into a file\n
      !!     Atom(Atom_class(:)): Structures with the atomic data\n
      !!    Bfield(Bfield_blass): Structure with magnetic field
      !!                          data\n
      !!      Flgsg(Fctsg_class): Structure with factorials and
      !!                          signs\n
      !!  JKQ(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                             over absorption profile\n
      !! JKQS(dcomplex(:,:,:,:)): Radiation field tensors integrated
      !!                            over emission profile\n
      !!  folder(character(500)): Output folder\n
      !!           iter(integer): Iteration index
      subroutine dump_jkq(Atom,Bfield,Flgsg,JKQ,JKQS,folder,iter)

      ! I/O
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class):: Flgsg
      character(len=500), intent(in):: folder
      integer, intent(in):: iter
      complex(kind=8), &
             dimension(-2:2,0:2,nxtran,Rz0:Rz1), intent(in):: JKQ
      complex(kind=8), &
             dimension(-2:2,0:2,nxtran,Rz0:Rz1), intent(in):: JKQS

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: ia,iz,itran,jtran,K,iQ
      complex(kind=8), dimension(-2:2,0:2,nz):: JKQaux,JKQSaux


      ! Sanity
      if (nA.lt.1) return

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_jkq'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_jkq_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.-2.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Write
      if (iter.eq.-2) then
        write(800,*) 'Initial JKQ'
      else if (iter.eq.-1) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'Final JKQ'
      else
        write(800,*) ''
        write(800,*) ''
        write(800,'("Iteration:",i3)') iter
      end if
      do ia=1,na
        write(800,'("  o Atom ",A)') Atom(ia)%element
        do itran=1,Atom(ia)%ntran
          jtran = itran + Atom(ia)%tshift
          do iz=Rz0,Rz1
            JKQaux(:,:,iz) = JKQ(:,:,jtran,iz)
            JKQSaux(:,:,iz) = JKQS(:,:,jtran,iz)
            if (Bfield%Bstrength(iz).gt.TINYB) then
              call fieldB(JKQaux(:,:,iz),1,Flgsg,-Bfield%Btheta(iz), &
                          -Bfield%Bphi(iz),-1)
              call fieldB(JKQSaux(:,:,iz),1,Flgsg,-Bfield%Btheta(iz),&
                          -Bfield%Bphi(iz),-1)
            end if
          end do
          do K=0,2
            do iQ=0,K
              if (iQ.eq.0) then
                do iz=Rz0,Rz1
                  if (Bfield%Bstrength(iz).gt.TINYB) then
                    write(800,'(i4,2(1x,i2),1x,i3,2(1x,es15.8,'// &
                               '17x,es15.8))') &
                      itran,K,iQ,iz,dble(JKQaux(iQ,K,iz)), &
                                    dble(JKQSaux(iQ,K,iz)), &
                                    dble(JKQ(iQ,K,jtran,iz)), &
                                    dble(JKQS(iQ,K,jtran,iz))
                  else
                    write(800,'(i4,2(1x,i2),1x,i3,1x,es15.8,'// &
                               '17x,es15.8)') &
                      itran,K,iQ,iz,dble(JKQaux(iQ,K,iz)), &
                                    dble(JKQSaux(iQ,K,iz))
                  end if
                end do
              else
                do iz=Rz0,Rz1
                  if (Bfield%Bstrength(iz).gt.TINYB) then
                    write(800,'(i4,2(1x,i2),1x,i3,8(1x,es15.8))') &
                      itran,K,iQ,iz,dble(JKQaux(iQ,K,iz)), &
                                    dimag(JKQaux(iQ,K,iz)), &
                                    dble(JKQSaux(iQ,K,iz)), &
                                    dimag(JKQSaux(iQ,K,iz)), &
                                    dble(JKQ(iQ,K,jtran,iz)), &
                                    dimag(JKQ(iQ,K,jtran,iz)), &
                                    dble(JKQS(iQ,K,jtran,iz)), &
                                    dimag(JKQS(iQ,K,jtran,iz))
                  else
                    write(800,'(i4,2(1x,i2),1x,i3,4(1x,es15.8))') &
                      itran,K,iQ,iz,dble(JKQaux(iQ,K,iz)), &
                                    dimag(JKQaux(iQ,K,iz)), &
                                    dble(JKQSaux(iQ,K,iz)), &
                                    dimag(JKQSaux(iQ,K,iz))
                  end if
                end do
              end if
            end do
          end do
        end do
      end do

      ! Close
      close(800)

      end subroutine dump_jkq

#endif
#ifdef DEBUGATMO
!#####################################################################
!#####################################################################
!#####################################################################

      !> Dump rho solution into a file\n
      !!       Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Bfield(Bfield_blass): Structure with magnetic field
      !!                         data\n
      !! folder(character(500)): Output folder\n
      !!          iter(integer): Iteration index
      subroutine dump_atmo(Atmo,Bfield,folder,iter)

      ! I/O
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      character(len=500), intent(in):: folder
      integer, intent(in):: iter

      ! Local
      character(len=500):: filename
      logical:: exists
      integer:: iz,iRz0,iRz1
      double precision, dimension(1:nz):: z,tau

      ! Get file name for 1D
      if (run_mode.eq.0) then

        filename = trim(folder)//'/debug_atmo'

      ! Get file name for rest
      else

        write(filename,'(A,I0.7)') trim(folder)//'/debug_atmo_', &
                                   icoords(3)

      end if

      !
      ! Exists?
      inquire(file=trim(filename), exist=exists)
      if(.not.exists.or.(iter.eq.0.and.run_mode.ne.-1))then
        open(800,file=trim(filename))
      else
        open(800,file=trim(filename),position='append')
      endif

      ! Scale
      if (ztau) then
        if (iter.eq.1) then
          z = Atmo%zalt
        else
          tau = 0d0
        end if
        tau = Atmo%z
      else
        z = Atmo%z
        if (iter.eq.1) then
          tau = Atmo%zalt
        else
          tau = 0d0
        end if
      end if

      ! Write Initial
      if (iter.eq.0) then
        write(800,*) 'Initial Atmosphere'
        iRz0 = 1
        iRz1 = nz
      ! Write calculating
      else if (iter.eq.1) then
        write(800,*) ''
        write(800,*) ''
        write(800,*) 'Final Atmosphere'
        iRz0 = Rz0
        iRz1 = Rz1
      end if

      ! Standard or final
      if (Atmo%typo.eq.0.or.iter.eq.1) then
        write(800,*) 'Number densities'
        write(800,*) '               z', &
                     '             tau', &
                     '               T', &
                     '              vx', &
                     '              vy', &
                     '              vz', &
                     '              vm', &
                     '               B', &
                     '          Btheta', &
                     '            Bphi', &
                     '              ne', &
                     '             nH0', &
                     '             nH1', &
                     '             nH2', &
                     '             nH3', &
                     '             nH4', &
                     '              np'
        do iz=iRz0,iRz1
          write(800,'(17(1x,es15.8))') &
            z(iz),tau(iz),Atmo%T(iz),Atmo%vx(iz)*1d6*c, &
                                     Atmo%vy(iz)*1d6*c, &
                                     Atmo%vz(iz)*1d6*c, &
                                     Atmo%vmi(iz)*1d6*c, &
                                     Bfield%Bstrength(iz), &
                                     Bfield%Btheta(iz), &
                                     Bfield%Bphi(iz), &
                                     Atmo%ne(iz), &
                                     Atmo%nH(iz,:)
        end do
      else if (Atmo%typo.eq.1.or.Atmo%typo.eq.2.or. &
               Atmo%typo.eq.3) then
        write(800,*) 'Electron number density'
        write(800,*) '               z', &
                     '             tau', &
                     '               T', &
                     '              vx', &
                     '              vy', &
                     '              vz', &
                     '              vm', &
                     '               B', &
                     '          Btheta', &
                     '            Bphi', &
                     '              ne'
        do iz=iRz0,iRz1
          write(800,'(11(1x,es15.8))') &
            z(iz),tau(iz),Atmo%T(iz),Atmo%vx(iz)*1d6*c, &
                                     Atmo%vy(iz)*1d6*c, &
                                     Atmo%vz(iz)*1d6*c, &
                                     Atmo%vmi(iz)*1d6*c, &
                                     Bfield%Bstrength(iz), &
                                     Bfield%Btheta(iz), &
                                     Bfield%Bphi(iz), &
                                     Atmo%ne(iz)
        end do
      else if (Atmo%typo.eq.4) then
        write(800,*) 'Gas pressure'
        write(800,*) '               z', &
                     '             tau', &
                     '               T', &
                     '              vx', &
                     '              vy', &
                     '              vz', &
                     '              vm', &
                     '               B', &
                     '          Btheta', &
                     '            Bphi', &
                     '              Pg'
        do iz=iRz0,iRz1
          write(800,'(11(1x,es15.8))') &
            z(iz),tau(iz),Atmo%T(iz),Atmo%vx(iz)*1d6*c, &
                                     Atmo%vy(iz)*1d6*c, &
                                     Atmo%vz(iz)*1d6*c, &
                                     Atmo%vmi(iz)*1d6*c, &
                                     Bfield%Bstrength(iz), &
                                     Bfield%Btheta(iz), &
                                     Bfield%Bphi(iz), &
                                     Atmo%Pg(iz)
        end do
      else if (Atmo%typo.eq.4) then
        write(800,*) 'Gas density'
        write(800,*) '               z', &
                     '             tau', &
                     '               T', &
                     '              vx', &
                     '              vy', &
                     '              vz', &
                     '              vm', &
                     '               B', &
                     '          Btheta', &
                     '            Bphi', &
                     '             rho'
        do iz=iRz0,iRz1
          write(800,'(11(1x,es15.8))') &
            z(iz),tau(iz),Atmo%T(iz),Atmo%vx(iz)*1d6*c, &
                                     Atmo%vy(iz)*1d6*c, &
                                     Atmo%vz(iz)*1d6*c, &
                                     Atmo%vmi(iz)*1d6*c, &
                                     Bfield%Bstrength(iz), &
                                     Bfield%Btheta(iz), &
                                     Bfield%Bphi(iz), &
                                     Atmo%rho(iz)
        end do
      end if

      ! Close
      close(800)

      end subroutine dump_atmo

#endif
!#####################################################################
!#####################################################################
!#####################################################################

      end module iosolution_mod
