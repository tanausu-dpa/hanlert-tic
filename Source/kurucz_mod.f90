      !> Kurucz lines management
      module kurucz_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     22/03/2019
!  Last version:
!     20/08/2025 V4.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     20/08/2025:    V4.0.2 - Added the possibility of a CPU not
!                             having work to do (TdPA)
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
!  kurucz_get
!    Read and organize the data in Kurucz line files
!
!  kurucz_bb
!    Compute the absorptivity and emissivity for a given frequency
!  of the included Kurucz lines
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use chemicaux_mod
      use commons_mod
      use math_mod
      use parameters_mod , only : CfA , KDT , dopp , hplanck , c , &
                                  PI , fktoJ , kb , convF , kdif
      use profile_mod
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Read and organize the data in Kurucz line files\n
      !!         Atom(Atom_class(:)): Structures with atomic data\n
      !!        Atomb(Atom_class(:)): Structures with atomic data for
      !!                              background atoms\n
      !!            Atmo(Atmo_class): Structure with atmospheric
      !!                              data\n
      !!       LTE(LTEline_class(:)): Structures with LTE line data\n
      !!  filenames(strarr_class(:)): List of file names\n
      !!                 NK(integer): Filenames sizes\n
      !!            omega(double(:)): Frequency array\n
      !!             MPID(MPI_class): Structure with MPI data\n
      !!                  Dw(double): Data for Doppler width\n
      !!                lDw(logical): If the value must be
      !!                              normalized\n
      !!        kurucz(kurucz_class): Structure with Kurucz line data
      subroutine kurucz_get(Atom,Atomb,Atmo,LTE,filenames,NK,omega, &
                            MPID,Dw,lDw,kurucz)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(LTEline_class), dimension(:), intent(in):: LTE
      type(Atmo_class), intent(in):: Atmo
      type(strarr_class), dimension(:), intent(in):: filenames
      type(MPI_class), intent(in):: MPID
      logical, intent(in):: lDw
      integer, intent(in):: NK
      double precision, intent(in):: Dw
      double precision, dimension(:), intent(in):: omega
      type(kurucz_class), intent(out):: kurucz

      ! Local

      character(len=1):: c0
      character(len=2):: element
      character(len=11):: c3,c4
      character(len=160):: cdump
      character(len=500):: lfilename

      logical:: found,foundb
      logical, dimension(:), allocatable:: contribute

      integer:: ifile,ios,itran,iftran,A,Z,nlin,ilin
      integer:: ia,iterml,itermu,iJl,iJu
      integer, dimension(NK):: lnlin

      real:: wave

      double precision:: E1,E2,J1,J2,omg,O0,O1,loggf
      double precision:: Grad,Gstk,Gvdw,DwT,Domg
      double precision:: lfac,ufac,Dfreq,rJl,rJu


      ! Routine name
      urou = 'kurucz_get'

      !
      ! Get factors for energy limits
      !
      lfac = 1d0 - kdif
      ufac = 1d0 + kdif

      ! If no frequencies
      if (MPID%nf(pid).lt.1) then

        ! No transitions
        kurucz%ntran = 0

        ! Return
        return

      end if

      !
      ! Get omega limits
      !
      O0 = omega(MPID%if0(pid))
      O1 = omega(MPID%if1(pid))

      !
      ! First, count number of raw lines
      !
      nlin = 0

      !
      ! First, count number of lines
      !

      ! For each filename, count lines
      do ifile=1,NK

        ! Count lines in individual files
        lnlin(ifile) = 0

        ! Get filename
        lfilename = filenames(ifile)%str

        ! Open file with Kurucz data
        open (100, file=trim(lfilename), status='old', &
              iostat=ios, err=1000)

        ! Until finished
        do while (.True.)

          ! Try to read one byte
          read(100,'(A)',iostat=ios) c0

          ! Good
          if (ios.eq.0) then

            ! Count line
            nlin = nlin + 1
            lnlin(ifile) = lnlin(ifile) + 1

          ! EoF
          else if (ios.lt.0) then

            ! Leave
            exit

          ! Error
          else

            ! Issue error
            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,.True.,.True.)
            return

          end if ! Could read

        end do ! Reading loop

        ! Close
        close(100)

      end do ! File loop

      ! Initialize flag to contribution
      allocate(contribute(nlin))
      contribute = .False.

      ! Initialize counters
      kurucz%ntran = 0
      ilin = 0

      ! For each filename, count lines
      do ifile=1,NK

        ! Get filename
        lfilename = filenames(ifile)%str

        ! Open file with Kurucz data
        open (100, file=trim(lfilename), status='old', &
              iostat=ios, err=1000)

        ! Until finished
        do while (.True.)

          ! Read line in buffer
          read(100,'(A)',iostat=ios) cdump

          ! Good
          if (ios.eq.0) then

            ! Advance line index
            ilin = ilin + 1

            ! Read up to energies from buffer
            read(cdump,'(f11.4,f7.3,i3,A1,i2,f12.3,f5.1,A11,'// &
                       'f12.3,f5.1)', iostat=ios) &
                wave,loggf,A,c0,Z,E1,J1,c3,E2,J2

            ! If valid line
            if (ios.eq.0) then

              ! Ensure positive energies
              E1 = abs(E1)*1d-5
              E2 = abs(E2)*1d-5

              ! Make E1 lower frequency if inverse order
              if (E1.gt.E2) then
                omg = E2
                E2 = E1
                E1 = omg
                omg = J2
                J2 = J1
                J1 = omg
              end if

              ! Get resonance
              omg = E2 - E1

              ! Get Doppler width
              if (lDw) then

                ! Complete Doppler width
                DwT = recallmass_ind(A)
                DwT = dopp*Dw/sqrt(DwT)

              ! If comes from outside
              else

                ! Copy Doppler width
                DwT = Dw

              end if

              ! Distance to line center
              Domg = omg*KDT*DwT

              ! If resonance out of CPU limits
              if (omg+Domg.lt.O0.or.omg-Domg.gt.O1) cycle

              ! Element and charges
              element = atom_index2char(A)
              Z = Z + 1

              ! Check atom not in list
              found = .False.

              ! If no partition function data
              if (Z.gt.Atmo%ele(A)%nstg) cycle

              ! For each active atom
              do ia=1,nA

                ! If not same atom, exit
                if (element.ne.Atom(ia)%Element) cycle

                ! If stage not in model atom
                if (minval(Atom(ia)%stage).gt.Z.or. &
                    maxval(Atom(ia)%stage).lt.Z) cycle

                ! If ML atom
                if (Atom(ia)%ML) then

                  ! For each lower term
                  do iterml=1,Atom(ia)%nMulti-1

                    ! Check ion
                    if (Atom(ia)%stage(iterml).gt.Z) exit
                    if (Atom(ia)%stage(iterml).ne.Z) cycle

                    ! Check can be same level
                    if (abs(Atom(ia)%rJval(1,iterml)-J1).gt.0.25) &
                      cycle
                    if (Atom(ia)%TRfreq(iterml).lt.lfac*E1) cycle
                    if (Atom(ia)%TRfreq(iterml).gt.ufac*E1) exit

                    ! For each upper term
                    do itermu=iterml+1,Atom(ia)%nMulti

                      ! Check ion
                      if (Atom(ia)%stage(itermu).gt. &
                          Atom(ia)%stage(iterml)) exit

                      ! Check can be same level
                      if (abs(Atom(ia)%rJval(1,itermu)-J2).gt.0.25) &
                        cycle
                      if (Atom(ia)%TRfreq(itermu).lt.lfac*E2) cycle
                      if (Atom(ia)%TRfreq(itermu).gt.ufac*E2) exit

                      ! Line was found, leave
                      found = .True.
                      exit

                    end do ! Upper terms

                    ! If found, leave
                    if (found) exit

                  end do ! Lower terms

                ! If MT atom
                else

                  ! For each lower term
                  do iterml=1,Atom(ia)%nMulti-1

                    ! Check ion
                    if (Atom(ia)%stage(iterml).gt.Z) exit
                    if (Atom(ia)%stage(iterml).ne.Z) cycle

                    ! For each lower level
                    do iJl=1,Atom(ia)%nJ(iterml)

                      ! Check can be same level
                      if (abs(Atom(ia)%rJval(iterml,iterml)-J1).gt. &
                          0.25) cycle
                      if (Atom(ia)%FSfreq(iJl,iterml).lt.lfac*E1) &
                        cycle
                      if (Atom(ia)%FSfreq(Ijl,iterml).gt.ufac*E1) &
                        exit

                      ! For each upper term
                      do itermu=iterml+1,Atom(ia)%nMulti

                        ! Check ion
                        if (Atom(ia)%stage(itermu).gt. &
                            Atom(ia)%stage(iterml)) exit

                        ! For each upper level
                        do iJu=1,Atom(ia)%nJ(itermu)

                          ! Check can be same level
                          if (abs(Atom(ia)%rJval(itermu,itermu)-J2) &
                              .gt.0.25) cycle
                          if (Atom(ia)%FSfreq(iJu,itermu).lt. &
                              lfac*E2) cycle
                          if (Atom(ia)%FSfreq(iJu,itermu).gt. &
                              ufac*E2) exit

                          ! Found line
                          found = .True.

                        end do ! Upper level

                        ! If found, leave
                        if (found) exit

                      end do ! Upper term

                      ! If found, leave
                      if (found) exit

                    end do ! Lower level

                    ! If found, leave
                    if (found) exit

                  end do ! Lower term

                end if ! ML or MT

                ! If found, leave
                if (found) exit

                !
                ! Check transitions
                !

                foundb = .False.

                ! If ML atom
                if (Atom(ia)%ML) then

                  ! For each transition
                  do itran=1,Atom(ia)%ntran

                    ! Terms
                    iterml = Atom(ia)%fst(itran)%iterml

                    ! Check stage
                    if (Z.ne.Atom(ia)%stage(iterml)) cycle

                    ! Get frequency of FS transition
                    Dfreq = Atom(ia)%Dfreq(itran)

                    ! Check blend
                    if ((Dfreq+Domg.gt.omg-Domg.and. &
                         Dfreq+Domg.lt.omg+Domg).or. &
                        (Dfreq-Domg.gt.omg-Domg.and. &
                         Dfreq-Domg.gt.omg+Domg)) then

                      ! Found blend
                      foundb = .True.

                      ! Issue warning
                      write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                      ' # Warning: Kurucz line for atom '// &
                      element//' in stage',Z,' is blended '// &
                      'with active line',itran,' at ', &
                      1d2/Dfreq,'nm. Cannot promise it '// &
                      'is not the same'
                      call verbose
                      exit

                    end if ! Blend

                    ! If blend found, leave
                    if (foundb) exit

                  end do ! Term transitions

                ! If MT
                else

                  ! For each transition
                  do itran=1,Atom(ia)%ntran

                    ! Terms
                    iterml = Atom(ia)%fst(itran)%iterml
                    itermu = Atom(ia)%fst(itran)%itermu

                    ! Check stage
                    if (Z.ne.Atom(ia)%stage(iterml)) cycle

                    ! Check blended FS transition
                    do iftran=1,Atom(ia)%fst(itran)%nt

                      ! Idenfity involved levels
                      iJu = Atom(ia)%fst(itran)%ilevelu(iftran)
                      iJl = Atom(ia)%fst(itran)%ilevell(iftran)

                      ! Get frequency of FS transition
                      Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                              Atom(ia)%FSfreq(iJl,iterml)

                      ! Check blend
                      if ((Dfreq+Domg.gt.omg-Domg.and. &
                           Dfreq+Domg.lt.omg+Domg).or. &
                          (Dfreq-Domg.gt.omg-Domg.and. &
                           Dfreq-Domg.gt.omg+Domg)) then

                        ! Blend found
                        foundb = .True.

                        ! Issue warning
                        write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                        ' # Warning: Kurucz line for atom '// &
                        element//' in stage',Z,' is blended '// &
                        'with active line',itran,' at ', &
                        1d2/Dfreq,'nm. Cannot promise it '// &
                        'is not the same'
                        call verbose
                        exit

                      end if ! Blenc

                    end do ! FS transitions

                    ! If bend found, leave
                    if (foundb) exit

                  end do ! Term transitions

                end if ! ML/MT

              end do ! Active atoms

              ! For each pasive atom
              do ia=1,nAb

                ! If it was found already, leave
                if (found) exit

                ! If not same atom, exit
                if (element.ne.Atomb(ia)%Element) cycle

                ! If stage not in model atom
                if (minval(Atomb(ia)%stage).gt.Z.or. &
                    maxval(Atomb(ia)%stage).lt.Z) cycle

                ! If ML atom
                if (Atomb(ia)%ML) then

                  ! For each lower term
                  do iterml=1,Atomb(ia)%nMulti-1

                    ! Check ion
                    if (Atomb(ia)%stage(iterml).gt.Z) exit
                    if (Atomb(ia)%stage(iterml).ne.Z) cycle

                    ! Check can be same level
                    if (abs(Atomb(ia)%rJval(1,iterml)-J1).gt.0.25) &
                      cycle
                    if (Atomb(ia)%TRfreq(iterml).lt.lfac*E1) cycle
                    if (Atomb(ia)%TRfreq(iterml).gt.ufac*E1) exit

                    ! For each upper term
                    do itermu=iterml+1,Atomb(ia)%nMulti

                      ! Check ion
                      if (Atomb(ia)%stage(itermu).gt. &
                          Atomb(ia)%stage(iterml)) exit

                      ! Check can be same level
                      if (abs(Atomb(ia)%rJval(1,itermu)-J2).gt.0.25) &
                        cycle
                      if (Atomb(ia)%TRfreq(itermu).lt.lfac*E2) cycle
                      if (Atomb(ia)%TRfreq(itermu).gt.ufac*E2) exit

                      ! Flag as found and leave
                      found = .True.
                      exit

                    end do ! Upper terms

                    ! If found, leave
                    if (found) exit

                  end do ! Lower terms

                ! If MT atom
                else

                  ! For each lower term
                  do iterml=1,Atomb(ia)%nMulti-1

                    ! Check ion
                    if (Atomb(ia)%stage(iterml).gt.Z) exit
                    if (Atomb(ia)%stage(iterml).ne.Z) cycle

                    ! For each lower level
                    do iJl=1,Atomb(ia)%nJ(iterml)

                      ! Check can be same level
                      if (abs(Atomb(ia)%rJval(iterml,iterml)-J1).gt. &
                          0.25) cycle
                      if (Atomb(ia)%FSfreq(iJl,iterml).lt.lfac*E1) &
                        cycle
                      if (Atomb(ia)%FSfreq(Ijl,iterml).gt.ufac*E1) &
                        exit

                      ! For each upper term
                      do itermu=iterml+1,Atomb(ia)%nMulti

                        ! Check ion
                        if (Atomb(ia)%stage(itermu).gt. &
                            Atomb(ia)%stage(iterml)) exit

                        ! For each upper level
                        do iJu=1,Atomb(ia)%nJ(itermu)

                          ! Check can be same level
                          if (abs(Atomb(ia)%rJval(itermu,itermu)-J2) &
                              .gt.0.25) cycle
                          if (Atomb(ia)%FSfreq(iJu,itermu).lt. &
                              lfac*E2) cycle
                          if (Atomb(ia)%FSfreq(iJu,itermu).gt. &
                              ufac*E2) exit

                          ! Flag as found
                          found = .True.

                        end do ! Upper level

                        ! If found, leave
                        if (found) exit

                      end do ! Upper term

                      ! If found, leave
                      if (found) exit

                    end do ! Lower level

                    ! If found, leave
                    if (found) exit

                  end do ! Lower term

                end if ! ML or MT

                ! If found, leave
                if (found) exit

                !
                ! Check transitions
                !

                foundb = .False.

                ! If ML atom
                if (Atomb(ia)%ML) then

                  ! For each lower term
                  do iterml=1,Atomb(ia)%nMulti-1

                    ! For each upper term
                    do itermu=iterml+1,Atomb(ia)%nMulti

                      ! Check stage
                      if (Z.ne.Atomb(ia)%stage(iterml)) cycle

                      ! Find the transition
                      itran = Atomb(ia)%irad(iterml,itermu)

                      ! If no transition, skip
                      if (itran.lt.1) cycle

                      ! Get frequency of the transition
                      Dfreq = Atomb(ia)%FSfreq(1,itermu) - &
                              Atomb(ia)%FSfreq(1,iterml)

                      ! Check blend
                      if ((Dfreq+Domg.gt.omg-Domg.and. &
                           Dfreq+Domg.lt.omg+Domg).or. &
                          (Dfreq-Domg.gt.omg-Domg.and. &
                           Dfreq-Domg.gt.omg+Domg)) then

                        ! Blend found
                        foundb = .True.

                        ! Issue warning
                        write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                        ' # Warning: Kurucz line for atom '// &
                        element//' in stage',Z,' is blended '// &
                        'with pasive line',itran,' at ', &
                        1d2/Dfreq,'nm. Cannot promise it '// &
                        'is not the same'
                        call verbose
                        exit

                      end if ! Found blend

                    end do ! Termu

                    ! If blend found, leave
                    if (foundb) exit

                  end do ! Term l

                ! If MT
                else

                  ! For each lower term
                  do iterml=1,Atomb(ia)%nMulti-1

                    ! For each upper term
                    do itermu=iterml+1,Atomb(ia)%nMulti

                      ! Check stage
                      if (Z.ne.Atomb(ia)%stage(iterml)) cycle

                      ! Find the transition
                      itran = Atomb(ia)%irad(iterml,itermu)

                      ! If no valid transition, skip
                      if (itran.lt.1) cycle

                      ! For each upper level
                      do iJu=1,Atomb(ia)%nJ(itermu)

                        ! Get angular momentum
                        rJu = Atomb(ia)%rJval(iJu,itermu)

                        ! For each lower level
                        do iJl=1,Atomb(ia)%nJ(iterml)

                          ! Get angular momentum
                          rJl = Atomb(ia)%rJval(iJl,iterml)

                          ! Check if allowed FS transition
                          if (abs(rJl - rJu).gt.1.1d0.or. &
                              abs(rJl + rJu).lt..1d0) cycle

                          ! Get frequency of the transition
                          Dfreq = Atomb(ia)%FSfreq(iju,itermu) - &
                                  Atomb(ia)%FSfreq(ijl,iterml)

                          ! Check blend
                          if ((Dfreq+Domg.gt.omg-Domg.and. &
                               Dfreq+Domg.lt.omg+Domg).or. &
                              (Dfreq-Domg.gt.omg-Domg.and. &
                               Dfreq-Domg.gt.omg+Domg)) then

                            ! Flag found blend
                            foundb = .True.

                            ! Issue warning
                            write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                            ' # Warning: Kurucz line for atom '// &
                            element//' in stage',Z,' is blended '// &
                            'with pasive line',itran,' at ', &
                            1d2/Dfreq,'nm. Cannot promise it '// &
                            'is not the same'
                            call verbose
                            exit

                          end if ! Blend

                        end do ! Jl

                        ! If blend found, leave
                        if (foundb) exit

                      end do ! Ju

                      ! If blend found, leave
                      if (foundb) exit

                    end do ! Termu

                    ! If blend found, leave
                    if (foundb) exit

                  end do ! Term l

                end if ! ML/MT

              end do ! Pasive atoms

              ! For each LTE line
              do ia=1,nLTEl

                ! If already found, leave
                if (found) exit

                ! If not same atom, exit
                if (A.ne.LTE(ia)%ele) cycle

                ! If not same stage
                if (LTE(ia)%stage.ne.Z) cycle

                ! Check can be same level
                if (abs(LTE(ia)%Jl-J1).gt.0.25) cycle
                if (LTE(ia)%El.lt.lfac*E1) cycle
                if (LTE(ia)%El.gt.ufac*E1) exit

                ! Check can be same level
                if (abs(LTE(ia)%Ju-J2).gt.0.25) cycle
                if (LTE(ia)%Eu.lt.lfac*E2) cycle
                if (LTE(ia)%Eu.gt.ufac*E2) exit

                ! Flag found
                found = .True.
                exit

              end do ! LTE lines

              ! If found, continue
              if (found) cycle

              ! Add to the Kurucz transition
              kurucz%ntran = kurucz%ntran + 1
              contribute(ilin) = .True.

            end if ! valid line

          ! EoF
          else if (ios.lt.0) then

            ! Finished this file
            exit

          ! Error
          else

            ! Issue error
            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,.True.,.True.)
            return

          end if ! Could read

        end do ! Reading loop

        ! Close
        close(100)

      end do ! File loop

      ! If no transitions, leave
      if (kurucz%ntran.le.0) return

      !
      ! Second, read proper data
      !

      ! Allocate structure
      allocate(kurucz%tran(kurucz%ntran))

      ! Initialize index
      itran = 0
      ilin = 0

      ! Initialize maximum frequency distance
      kurucz%MDomg = 0d0

      ! For each filename, count lines
      do ifile=1,NK

        ! Get filename
        lfilename = filenames(ifile)%str

        ! Open file with Kurucz data
        open (100, file=trim(lfilename), status='old', &
              iostat=ios, err=1000)

        ! Run until EoF
        do while (.True.)

          ! Advance index
          ilin = ilin + 1

          ! If beyond
          if (ilin.gt.lnlin(ifile)) then

            ! Close and leave
            close(100)
            exit

          end if ! Beyond

          ! If no contribution
          if (.not.contribute(ilin)) then

            ! Skip line
            read(100,'(A1)') c0
            cycle

          end if

          ! Read data
          read(100,'(f11.4,f7.3,i3,A1,i2,f12.3,f5.1,A11,'// &
                   'f12.3,f5.1,A11,3f6.3)', iostat=ios) &
            wave,loggf,A,c0,Z,E1,J1,c3,E2,J2,c4,Grad,Gstk,Gvdw

          ! Good
          if (ios.eq.0) then

            ! Ensure positive energies
            E1 = abs(E1)*1d-5
            E2 = abs(E2)*1d-5

            ! Make E1 lower frequency
            if (E1.gt.E2) then
              omg = E2
              E2 = E1
              E1 = omg
              omg = J2
              J2 = J1
              J1 = omg
            end if

            ! Get resonance
            omg = E2 - E1

            ! Get Doppler width
            if (lDw) then

              ! Recalculate
              DwT = recallmass_ind(A)
              DwT = dopp*Dw/sqrt(DwT)

            else

              ! Get from input
              DwT = Dw

            end if

            ! Distance to line center
            Domg = omg*KDT*DwT

            ! If resonance not within CPU limits
            if (omg+Domg.lt.O0.or.omg-Domg.gt.O1) cycle

            ! Update maximum
            if (Domg.gt.kurucz%MDomg) kurucz%MDomg = Domg

            ! Advance indexes
            itran = itran + 1
            Z = Z + 1

            ! Add to RAM
            MRAMc = MRAMc + 1d-6*sizeof(kurucz%tran(itran))

            kurucz%tran(itran)%O0 = omg-Domg
            kurucz%tran(itran)%O1 = omg+Domg
            kurucz%tran(itran)%A = A
            kurucz%tran(itran)%Z = Z
            kurucz%tran(itran)%Ei = E1
            kurucz%tran(itran)%Dfreq = omg
            kurucz%tran(itran)%gu = (2d0*J2 + 1d0)
            kurucz%tran(itran)%gl = (2d0*J1 + 1d0)

            ! gf
            loggf = 10d0**loggf

            ! gf*Cfa/gu
            loggf = loggf*CfA/kurucz%tran(itran)%gu

            ! Get Aul
            kurucz%tran(itran)%Aul = loggf* &
                                     kurucz%tran(itran)%Dfreq* &
                                     kurucz%tran(itran)%Dfreq

            ! Get broadening
            kurucz%tran(itran)%Grad = 10d0**Grad
            kurucz%tran(itran)%Gstk = 10d0**Gstk
            kurucz%tran(itran)%Gvdw = 10d0**Gvdw

            ! Check abundance
            found = .False.

            ! Get string
            element = atom_index2char(A)

            ! For each active atom
            do ia=1,nA

              ! If same atom, exit
              if (element.eq.Atom(ia)%Element) then

                ! Get abundance from model
                kurucz%tran(itran)%abun = Atom(ia)%abun
                found = .True.
                exit

              end if

            end do ! Active atoms

            ! If no active
            if (.not.found) then

              ! For each pasive atom
              do ia=1,nAb

                ! If same atom, exit
                if (element.eq.Atomb(ia)%Element) then

                  ! Get abundance from model
                  kurucz%tran(itran)%abun = Atomb(ia)%abun
                  found = .True.
                  exit

                end if

              end do ! Pasive atoms

            end if ! Not active

            ! If did not find abundance from database
            if (.not.found) &
              kurucz%tran(itran)%abun = Atmo%abund(A)

          ! EoF
          else if (ios.lt.0) then

            ! Leave
            exit

          ! Error
          else

            ! Issue error
            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,.True.,.True.)
            return

          end if ! Read status

        end do ! Reading loop

        ! Close
        close(100)

      end do ! File loop

      return

1000  umsg = 'Error opening kurucz file '//trim(lfilename)
      call abortedS(umsg,urou,.True.,.True.)
      return

      end subroutine kurucz_get

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the absorptivity and emissivity for a given frequency
      !! of the included Kurucz lines\n
      !!          freq(double): Frequency\n
      !!  Kurucz(kurucz_class): Structure with Kurucz line data
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!          T(double(:)): Temperature\n
      !!         ne(double(:)): Electron number density\n
      !!       nH(double(:,:)): Hydrogen number density\n
      !!        vmi(double(:)): Microturbulent velocity\n
      !!       vfac(double(:)): Doppler shift factor\n
      !!          iz0(integer): First height index to consider\n
      !!          iz1(integer): Last height index to consider\n
      !!        fline(logical): Indicate if there was a line in this
      !!                        frequency\n
      !!        eta(double(:)): Absorptivity\n
      !!        eps(double(:)): Emissivity
      subroutine kurucz_bb(freq,Kurucz,Atmo,vfac,iz0,iz1, &
                           fline,eta,eps)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Kurucz_class), intent(in):: Kurucz
      logical, intent(out):: fline
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: vfac
      double precision, dimension(iz0:iz1), intent(out):: eta,eps

      ! Local

      type(catm_class):: Atom

      character(len=2):: element

      integer:: itran, iz, izp

      double precision:: rmass,Dw,adamp,prof,c0,c1,sqrtpi,nHI
      double precision:: arg,nl,nu,etal,epsl,Nion,exu,pEl,sig
      double precision, dimension(1):: frc


      ! Initialize the line found variable
      fline = .False.

      ! hc/4pi and square root of pi
      c0 = hplanck*c*1d12*.25d0/PI
      sqrtpi = sqrt(PI)

      ! Initialize RT coefficients
      eta = 0d0
      eps = 0d0

      ! Go through all the lines
      do itran=1,Kurucz%ntran

        ! Check within limits
        if (Kurucz%tran(itran)%O0.gt.freq*maxval(vfac).or. &
            Kurucz%tran(itran)%O1.lt.freq*minval(vfac)) cycle

        ! If not, we have found a line
        fline = .True.

        ! Get element
        element = atom_index2char(kurucz%tran(itran)%A)

        ! Free atomic data
        if (allocated(Atom%pf)) then
          MRAMc = MRAMc - 1d-6*sizeof(Atom%pf)
          MRAMc = MRAMc - 1d-6*sizeof(Atom%Eion)
          deallocate(Atom%pf)
          deallocate(Atom%Eion)
        end if

        ! Get partition function for this atom
        call getpf(element,Atom%nstg,Atom%pf, &
                   Atom%Eion,Atmo)

        ! Get mass
        rmass = recallmass_ind(Kurucz%tran(itran)%A)

        ! Get RT quantities
        pEl = convF*2d21*c*Kurucz%tran(itran)%Dfreq* &
                           Kurucz%tran(itran)%Dfreq* &
                           Kurucz%tran(itran)%Dfreq

        ! Constant part of the RT coefficients
        ! 1d14 = 1d9 (m^3 -> cm^3) * 1d5 (10^5 cm^-1 -> cm^-1)
        c1 = c0*Kurucz%tran(itran)%Dfreq*1d14* &
             Kurucz%tran(itran)%Aul/sqrtpi

        ! For each height
        do iz=iz0,iz1

          ! Corrected index
          izp = iz - iz0 + 1

          ! hydrogen number densities
          nHI = sum(Atmo%nH(iz,1:5))

          ! Damping
          Dw = dopp*sqrt(Atmo%T(iz))/sqrt(rmass)
          Dw = Kurucz%tran(itran)%Dfreq*sqrt(Dw*Dw + &
                                             Atmo%vmi(iz)**2d0)
          adamp = Kurucz%tran(itran)%Grad + &
                  Kurucz%tran(itran)%Gstk*Atmo%ne(iz) + &
                  Kurucz%tran(itran)%Gvdw*nHI
          adamp = adamp*1d-16/c/(4d0*PI)
          adamp = adamp/Dw

          ! Profile
          call voigtI((Kurucz%tran(itran)%Dfreq - freq*vfac(iz))/Dw, &
                       adamp,prof)

          ! Ionization fraction
          call getfrc(Atom%nstg,Atom%pf(:,izp),Atom%Eion, &
                      Atmo%T(iz),Atmo%ne(iz), &
                      Kurucz%tran(itran)%Z,frc)

          ! Compute density of ion
          Nion = frc(1)*Kurucz%tran(itran)%abun*Atmo%nHT(iz)

          !
          ! Compute level population
          !

          ! Exponential argument
          arg = -Kurucz%tran(itran)%Ei*fktoJ/kb/Atmo%T(iz) - &
                 Atom%pf(Kurucz%tran(itran)%Z,izp)

          ! Exponential
          if (arg.lt.0d0) then
            arg = -arg
            exu = diexp(arg)
          else
            exu = ddexp(arg)
          end if

          ! Population lower level times gu/gl
          nl = Nion*exu*Kurucz%tran(itran)%gu

          ! Population upper level
          exu = diexp(Kurucz%tran(itran)%Dfreq*fktoJ/kb/Atmo%T(iz))
          nu = nl*exu

          ! 'cross-section'
          sig = c1/Dw

          ! Absorptivity and emissivity
          etal = sig*(nl - nu)/pEl
          epsl = sig*nu

          ! Scattering
          ! TODO or not, because not sure I want to produce
          ! scattering polarization this way

          ! Add contribuions to absorptivity and emissivity
          eta(iz) = eta(iz) + etal*prof
          eps(iz) = eps(iz) + epsl*prof

        end do ! Heights
      end do ! Kurucz lines

      return

      end subroutine kurucz_bb

!#####################################################################
!#####################################################################
!#####################################################################

      end module kurucz_mod
