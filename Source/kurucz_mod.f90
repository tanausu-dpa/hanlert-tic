      !> Kurucz lines management
      module kurucz_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     03/22/2019
!  Last version:
!     08/07/2023 V3.0.4
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.4 - Exclude lines that have been loaded as
!                             LTE lines when checking LTE lines for
!                             the background (TdPA)
!
!     07/03/2023:    V3.0.3 - Bugfix: If there is an error reading,
!                             the routine needs to return (TdPA)
!
!     07/27/2022:    V3.0.2 - Renamed MPI to MPID (TdPA)
!
!     06/29/2022:    V3.0.1 - The calls to kurucz_get and kurucz_bb
!                             now require Atmo and not the resource
!                             argument (TdPA)
!                           - In kurucz_bb the specific thermodynamic
!                             variables have been removed from the
!                             arguments, as they are in Atmo (TdPA)
!                           - The call to checkpf is no longer needed
!                             because Atmo has the data (TdPA)
!                           - The call to getpf required Atmo instead
!                             of the resource argument (TdPA)
!                           - Abundances are extracted from Atmo and
!                             not from the hard-coded table in the
!                             chemicaux_mod module (TdPA)
!
!     06/29/2022:    V3.0.0 - Changed global version (TdPA)
!
!     12/16/2021:    V2.0.2 - Bugfix: With several files, the reading
!                             failed because the running index was
!                             compared against the total number of
!                             lines in all files, and not the lines
!                             in the file being currectly read (TdPA)
!
!     03/23/2021:    V2.0.1 - Changed call to abortedS (TdPA)
!
!     03/17/2021:    V2.0.0 - Changed global version (TdPA)
!
!     11/05/2020:    V1.0.3 - Passing full atmospheric arrays to
!                             kurucz_bb to avoid creation of temporal
!                             arrays with domain decomposition (TdPA)
!
!     03/05/2020:    V1.0.2 - Kurucz receives and uses total hydrogen
!                             density, not just atomic (TdPA)
!
!     11/19/2019:    V1.0.1 - Removed checks in allocate and
!                             deallocate calls (TdPA)
!
!     03/22/2019:    V1.0.0 - First version (TdPA)
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
!    This program handles Kurucz files
!
!  kurucz_get:
!    Reads and organized the data in Kurucz line files
!
!  kurucz_bb:
!    Computes the absorptivity and emissivity for a given frequency
!  of the Kurucz lines stored in the previously built structure
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

      !> Gets Kurucz line data\n
      !!          Atom(Atom_class): Structure with the atomic data\n
      !!         Atomb(Atom_class): Structure with the atomic data for
      !!                            background opacities\n
      !!          Atmo(Atmo_class): Structure with atmospheric data\n
      !!        LTE(LTEline_class): Structure with the LTE line data\n
      !!   filenames(strarr_class): List of file names\n
      !!               NK(integer): Filenames size\n
      !!          omega(dfloat(:)): Frequency array\n
      !!           MPID(MPI_class): Structure with MPI data\n
      !!                Dw(dfloat): Data for Doppler width\n
      !!              lDw(logical): If the value must be normalized\n
      !!      kurucz(kurucz_class): Structure with Kurucz data
      subroutine kurucz_get(Atom, Atomb, Atmo, LTE, filenames, NK, &
                            omega, MPID, Dw, lDw, kurucz)

      ! I/O

      type(Atom_class), dimension(:), intent(in):: Atom
      type(Atom_class), dimension(:), intent(in):: Atomb
      type(LTEline_class), dimension(:), intent(in):: LTE
      type(Atmo_class), intent(in):: Atmo
      type(strarr_class), dimension(:), intent(in):: filenames
      type(MPI_class):: MPID
      logical, intent(in):: lDw
      integer, intent(in):: NK
      double precision, intent(in):: Dw
      double precision, dimension(:), intent(in):: omega
      type(kurucz_class), intent(out):: kurucz

      ! Local

      character(len=500):: lfilename
      character(len=160):: cdump
      character(len=1):: c0
      character(len=11):: c3
      character(len=11):: c4
      character(len=2):: element
      logical:: found, foundb
      logical, dimension(:), allocatable:: contribute
      integer:: ifile, ios, itran, iftran, A, Z, nlin, ilin
      integer:: ia, iterml, itermu, iJl, iJu, i, i1
      integer, dimension(NK):: lnlin
      real:: wave
      double precision:: E1, E2, J1, J2, omg, O0, O1, loggf
      double precision:: Grad, Gstk, Gvdw, DwT, Domg
      double precision:: lfac, ufac, Dfreq, rJl, rJu

      urou = 'kurucz_get'

      !
      ! Get factors for energy limits
      !
      lfac = 1d0 - kdif
      ufac = 1d0 + kdif

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

        do while (.True.)

          read(100,'(A)',iostat=ios) c0

          ! Good
          if (ios.eq.0) then

            nlin = nlin + 1
            lnlin(ifile) = lnlin(ifile) + 1

          ! EoF
          else if (ios.lt.0) then

            exit

          ! Error
          else

            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,-1,.True.,.True.)
            return

          end if

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

        do while (.True.)

          read(100,'(A)',iostat=ios) cdump

          ! Good
          if (ios.eq.0) then

            ! Advance line index
            ilin = ilin + 1

            ! Read up to energies
            read(cdump,'(f11.4,f7.3,i3,A1,i2,f12.3,f5.1,A11,'// &
                       'f12.3,f5.1)', iostat=ios) &
                wave,loggf,A,c0,Z,E1,J1,c3,E2,J2

            ! If valid line
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
              omg = E2-E1

              ! Get Doppler width
              if (lDw) then

                DwT = recallmass_ind(A)
                DwT = dopp*Dw/sqrt(DwT)

              else

                DwT = Dw

              end if

              ! Distance to line center
              Domg = omg*KDT*DwT

              ! If resonance out of CPU limits
              if (omg+Domg.lt.O0.or.omg-Domg.gt.O1) cycle

              ! Check atom not in list
              element = atom_index2char(A)
              Z = Z + 1
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

                      found = .True.
                      exit

                    end do ! Upper terms
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

                          found = .True.

                        end do ! Upper level
                        if (found) exit
                      end do ! Upper term
                      if (found) exit
                    end do ! Lower level
                    if (found) exit
                  end do ! Lower term

                end if ! ML or MT

                if (found) exit

                !
                ! Check transitions
                !

                foundb = .False.

                ! If ML atom
                if (Atom(ia)%ML) then

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
                      foundb = .True.
                      write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                      ' # Warning: Kurucz line for atom '// &
                      element//' in stage',Z,' is blended '// &
                      'with active line',itran,' at ', &
                      1d2/Dfreq,'nm. Cannot promise it '// &
                      'is not the same'
                      call verbose
                      exit
                    end if

                    if (foundb) exit
                  end do ! Term transitions

                ! If MT
                else

                  do itran=1,Atom(ia)%ntran

                    ! Terms
                    iterml = Atom(ia)%fst(itran)%iterml
                    itermu = Atom(ia)%fst(itran)%itermu

                    ! Check stage
                    if (Z.ne.Atom(ia)%stage(iterml)) cycle

                    ! Check blended FS transition
                    do iftran=1,Atom(ia)%fst(itran)%nt

                      ! Idenfity involved levels
                      iJu = -1
                      do i=1,Atom(ia)%nJ(iterml)
                        do i1=1,Atom(ia)%nJ(itermu)
                          if (Atom(ia)%fst(itran)%irad(i1,i).eq. &
                              iftran) then
                            iJl = i
                            iJu = i1
                            exit
                          end if
                        end do
                        if (iJu.ge.0) exit
                      end do

                      ! Get frequency of FS transition
                      Dfreq = Atom(ia)%FSfreq(iJu,itermu) - &
                              Atom(ia)%FSfreq(iJl,iterml)

                      ! Check blend
                      if ((Dfreq+Domg.gt.omg-Domg.and. &
                           Dfreq+Domg.lt.omg+Domg).or. &
                          (Dfreq-Domg.gt.omg-Domg.and. &
                           Dfreq-Domg.gt.omg+Domg)) then
                        foundb = .True.
                        write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                        ' # Warning: Kurucz line for atom '// &
                        element//' in stage',Z,' is blended '// &
                        'with active line',itran,' at ', &
                        1d2/Dfreq,'nm. Cannot promise it '// &
                        'is not the same'
                        call verbose
                        exit
                      end if

                    end do ! FS transitions
                    if (foundb) exit
                  end do ! Term transitions

                end if ! ML/MT

              end do ! Active atoms

              ! For each pasive atom
              do ia=1,nAb

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

                      found = .True.
                      exit

                    end do ! Upper terms
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

                          found = .True.

                        end do ! Upper level
                        if (found) exit
                      end do ! Upper term
                      if (found) exit
                    end do ! Lower level
                    if (found) exit
                  end do ! Lower term

                end if ! ML or MT

                if (found) exit

                !
                ! Check transitions
                !

                foundb = .False.

                ! If ML atom
                if (Atomb(ia)%ML) then

                  do iterml=1,Atomb(ia)%nMulti-1
                    do itermu=iterml+1,Atomb(ia)%nMulti

                      ! Check stage
                      if (Z.ne.Atomb(ia)%stage(iterml)) cycle

                      ! Find the transition
                      itran = Atomb(ia)%irad(iterml,itermu)

                      if (itran.lt.1) cycle

                      ! Get frequency of the transition
                      Dfreq = Atomb(ia)%FSfreq(1,itermu) - &
                              Atomb(ia)%FSfreq(1,iterml)

                      ! Check blend
                      if ((Dfreq+Domg.gt.omg-Domg.and. &
                           Dfreq+Domg.lt.omg+Domg).or. &
                          (Dfreq-Domg.gt.omg-Domg.and. &
                           Dfreq-Domg.gt.omg+Domg)) then
                        foundb = .True.
                        write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                        ' # Warning: Kurucz line for atom '// &
                        element//' in stage',Z,' is blended '// &
                        'with pasive line',itran,' at ', &
                        1d2/Dfreq,'nm. Cannot promise it '// &
                        'is not the same'
                        call verbose
                        exit
                      end if

                    end do ! Termu
                    if (foundb) exit
                  end do ! Term l

                ! If MT
                else

                  do iterml=1,Atomb(ia)%nMulti-1
                    do itermu=iterml+1,Atomb(ia)%nMulti

                      ! Check stage
                      if (Z.ne.Atomb(ia)%stage(iterml)) cycle

                      ! Find the transition
                      itran = Atomb(ia)%irad(iterml,itermu)

                      if (itran.lt.1) cycle

                      ! For each pair of FS levels
                      do iJu=1,Atomb(ia)%nJ(itermu)

                        rJu = Atomb(ia)%rJval(iJu,itermu)

                        do iJl=1,Atomb(ia)%nJ(iterml)

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
                            foundb = .True.
                            write(umsg,'(A,1x,i2,A,1x,i5,A,f7.3,A)') &
                            ' # Warning: Kurucz line for atom '// &
                            element//' in stage',Z,' is blended '// &
                            'with pasive line',itran,' at ', &
                            1d2/Dfreq,'nm. Cannot promise it '// &
                            'is not the same'
                            call verbose
                            exit
                          end if

                        end do ! Jl
                        if (foundb) exit
                      end do ! Ju
                      if (foundb) exit
                    end do ! Termu
                    if (foundb) exit
                  end do ! Term l

                end if ! ML/MT

              end do ! Pasive atoms


              ! For each LTE line
              do ia=1,nLTEl

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

                found = .True.
                exit

              end do ! LTE lines


              if (found) cycle

              kurucz%ntran = kurucz%ntran + 1
              contribute(ilin) = .True.

            end if ! valid line

          ! EoF
          else if (ios.lt.0) then

            exit

          ! Error
          else

            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,-1,.True.,.True.)
            return

          end if

        end do ! Reading loop

        ! Close
        close(100)

      end do ! File loop

      !
      ! If no transitions, leave
      !
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

        do while (.True.)

          ilin = ilin + 1

          if (ilin.gt.lnlin(ifile)) then
            close(100)
            exit
          end if

          if (.not.contribute(ilin)) then

            read(100,'(A1)') c0
            cycle

          end if

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
            omg = E2-E1

            ! Get Doppler width
            if (lDw) then

              DwT = recallmass_ind(A)
              DwT = dopp*Dw/sqrt(DwT)

            else

              DwT = Dw

            end if

            ! Distance to line center
            Domg = omg*KDT*DwT

            ! If resonance not within CPU limits
            if (omg+Domg.lt.O0.or.omg-Domg.gt.O1) cycle

            ! Update maximum
            if (Domg.gt.kurucz%MDomg) kurucz%MDomg = Domg

            itran = itran + 1
            Z = Z + 1

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
            kurucz%tran(itran)%Aul = loggf* &
                                     kurucz%tran(itran)%Dfreq* &
                                     kurucz%tran(itran)%Dfreq
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

                  kurucz%tran(itran)%abun = Atomb(ia)%abun
                  found = .True.
                  exit

                end if

              end do ! Pasive atoms

            end if ! Not active

            ! If did not find abundance
            if (.not.found) &
              kurucz%tran(itran)%abun = Atmo%abund(A)

          ! EoF
          else if (ios.lt.0) then

            exit

          ! Error
          else

            umsg = 'Error reading kurucz file '//trim(lfilename)
            close(100)
            call abortedS(umsg,urou,-1,.True.,.True.)
            return

          end if

        end do ! Reading loop

        ! Close
        close(100)

      end do ! File loop

      return

1000  umsg = 'Error opening kurucz file '//trim(lfilename)
      call abortedS(umsg,urou,-1,.True.,.True.)
      return

      end subroutine kurucz_get

!#####################################################################
!#####################################################################
!#####################################################################

      !> Computes atomic bound-bound continuum radiation transfer
      !! coefficients from Kurucz line data\n
      !!             freq(dfloat): Frequency\n
      !!     Kurucz(Kurucz_class): Structure with Kurucz data\n
      !!         Atmo(Atmo_class): Structure with atmospheric data\n
      !!             T(dfloat(:)): Temperature\n
      !!            ne(dfloat(:)): Electron density\n
      !!          nH(dfloat(:,:)): Hydrogen density\n
      !!           vmi(dfloat(:)): Microturbulent velocity\n
      !!          vfac(dfloat(:)): Doppler shift factor\n
      !!             iz0(integer): First height index for this CPU\n
      !!             iz1(integer): Last height index for this CPU\n
      !!           fline(logical): Bool that tells if a line was found
      !!                           in this frequency\n
      !!           eta(dfloat(:)): Absorptivity\n
      !!           eps(dfloat(:)): Emissivity
      subroutine kurucz_bb(freq,Kurucz,Atmo,vfac,iz0,iz1, &
                           fline,eta,eps)
      ! I/O

      type(Atmo_class), intent(in):: Atmo
      type(Kurucz_class), intent(in):: Kurucz
      logical, intent(out):: fline
      integer, intent(in):: iz0,iz1
      double precision, intent(in):: freq
      double precision, dimension(iz0:iz1), intent(in):: vfac
      double precision, dimension(iz0:iz1), intent(out):: eta, eps

      ! Local

      type(catm_class):: Atom

      character(len=2):: element
      integer:: itran, iz, izp
      double precision:: rmass, Dw, adamp, prof, c0, c1, sqrtpi
      double precision:: arg, nl, nu, etal, epsl, Nion, exu, pEl, sig
      double precision:: nHI
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

        ! Get the partition function
        if (allocated(Atom%pf)) then
          deallocate(Atom%pf)
          deallocate(Atom%Eion)
        end if
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

          call getfrc(Atom%nstg,Atom%pf(:,izp),Atom%Eion, &
                      Atmo%T(iz),Atmo%ne(iz), &
                      Kurucz%tran(itran)%Z,frc)

          ! Compute density of ion
          Nion = frc(1)*Kurucz%tran(itran)%abun*Atmo%nHT(iz)

          ! Compute level population
          arg = -Kurucz%tran(itran)%Ei*fktoJ/kb/Atmo%T(iz) - &
                 Atom%pf(Kurucz%tran(itran)%Z,izp)
          if (arg.lt.0d0) then
            arg = -arg
            exu = diexp(arg)
          else
            exu = ddexp(arg)
          end if
          ! population lower level times gu/gl
          nl = Nion*exu*Kurucz%tran(itran)%gu
          exu = diexp(Kurucz%tran(itran)%Dfreq*fktoJ/kb/Atmo%T(iz))
          ! population upper level
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

        end do !heights
      end do ! Kurucz lines

      return

      end subroutine kurucz_bb

!#####################################################################
!#####################################################################
!#####################################################################

      end module kurucz_mod
