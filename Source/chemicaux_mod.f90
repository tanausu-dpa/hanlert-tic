      !> Abundances, partition functions, and ion fractions
      module chemicaux_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     19/04/2017
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
!    Implement better equation of state
!
!#####################################################################
!#####################################################################
!
!  Data:
!
!  recallnumber
!    Return number of atoms in database
!
!  recallabund
!    Return hard-coded abundance of an atom specified by name in
!  N/NH
!
!  recallabund_ind
!    Return hard-coded abundance of an atom specified by index in
!  N/NH
!
!  recallmass_ind
!    Return hard-coded mass of an atom specified by index in AMU
!
!  atom_char2index
!    Return the index of an atom specified by name
!
!  atom_index2char
!    Return name of an atom specified by index
!
!  getpf_T
!    Get partition function, ionization energy, and number of stages
!  of the atom specified by name at a given temperature
!
!  getpf
!    Get partition function, ionization energy, and number of stages
!  of the atom specified by name at all heights
!
!  getfrc
!    Get ionization fraction for the specified stages at a given
!  temperature and for a given electron density
!
!  getdfrc
!    Get ionization fraction derivative at a given temperature and
!  for a given electron density
!
!  eqstate_known
!    Get pressure and mass density for known number densities at every
!  height
!
!  eqstate_gas
!    Get pressure and mass density for known gas pressure at every
!  height
!
!  eqstate_ele
!    Get pressure and mass density for known electron pressure at
!  every height
!
!  partial_press_known
!    Compute partial pressures of ions from known ionization data
!
!  partial_press_Pg
!    Compute partial pressures of ions from known gas pressure
!
!  partial_press_Pe
!    Compute partial pressures of ions from known electron pressure
!
!  metal_fraction
!    Compute ion fraction from populations in an atomic model
!
!  chianti_fraction
!    Compute ionization fraction of stages in an atomic model from
!  CHIANTI data for a given temperature
!
!  H_fraction
!    Compute ion fraction from populations in Hydrogen model
!
!  LTEiz
!    Computes LTE populations at a given height
!
!  set_densities
!    Set the atmospheric quantities from the computed partial
!  pressures
!
!  set_Hdensities
!    Set the hydrogen density atmospheric quantities from the computed
!  partial pressures
!
!  fsaha
!    Compute Saha factor for equation of state
!
!  moldata_ind
!    Get molecular quantities for equation of state. Data taken from
!  the SIR code
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use inter_mod
      use math_mod
      use parameters_mod , only : hplanck , PI , fktoJ , me , kb , &
                                  TINYFRC , fktoev, qel, pi4eps0
      use types_mod

      ! Parameters

      ! Number of atoms contributing in partial_pressure
      integer, parameter:: natom = 28

      ! Number of molecules contributing in partial_pressure
      integer, parameter:: nmol = 2

      ! Maxumum iterations in eqstate_gas
      integer, parameter:: maxitereq = 50

      ! MRC tolerance in eqstate_gas
      double precision, parameter:: eps = 1d-2

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return number of atoms in database\n
      !!  x(integer):: Dummy argument
      integer function recallnumber(x)

      integer, intent(in):: x

      recallnumber = x
      recallnumber = 99

      end function recallnumber

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return hard-coded abundance of an atom specified by name in
      !! N/NH\n
      !!  atm(character(:)): ID name of the atom
      double precision function recallabund(atm)

      ! I/O

      character(len=2), intent(in):: atm


      ! Select return value given the name
      select case (atm)
        case (' H')
          recallabund = 12.00
        case ('HE')
          recallabund = 10.99
        case ('LI')
          recallabund =  1.16
        case ('BE')
          recallabund =  1.15
        case (' B')
          recallabund =  2.6
        case (' C')
          recallabund =  8.39
        case (' N')
          recallabund =  8.00
        case (' O')
          recallabund =  8.66
        case (' F')
          recallabund =  4.56
        case ('NE')
          recallabund =  8.09
        case ('NA')
          recallabund =  6.33
        case ('MG')
          recallabund =  7.58
        case ('AL')
          recallabund =  6.47
        case ('SI')
          recallabund =  7.55
        case (' P')
          recallabund =  5.45
        case (' S')
          recallabund =  7.21
        case ('CL')
          recallabund =  5.50
        case ('AR')
          recallabund =  6.56
        case (' K')
          recallabund =  5.12
        case ('CA')
          recallabund =  6.36
        case ('SC')
          recallabund =  3.10
        case ('TI')
          recallabund =  4.99
        case (' V')
          recallabund =  4.00
        case ('CR')
          recallabund =  5.67
        case ('MN')
          recallabund =  5.39
        case ('FE')
          recallabund =  7.44
        case ('CO')
          recallabund =  4.92
        case ('NI')
          recallabund =  6.25
        case ('CU')
          recallabund =  4.21
        case ('ZN')
          recallabund =  4.60
        case ('GA')
          recallabund =  2.88
        case ('GE')
          recallabund =  3.41
        case ('AS')
          recallabund =  2.37
        case ('SE')
          recallabund =  3.35
        case ('BR')
          recallabund =  2.63
        case ('KR')
          recallabund =  3.23
        case ('RB')
          recallabund =  2.60
        case ('SR')
          recallabund =  2.90
        case (' Y')
          recallabund =  2.24
        case ('ZR')
          recallabund =  2.60
        case ('NB')
          recallabund =  1.42
        case ('MO')
          recallabund =  1.92
        case ('TC')
          recallabund = -7.96
        case ('RU')
          recallabund =  1.84
        case ('RH')
          recallabund =  1.12
        case ('PD')
          recallabund =  1.69
        case ('AG')
          recallabund =  0.94
        case ('CD')
          recallabund =  1.86
        case ('IN')
          recallabund =  1.66
        case ('SN')
          recallabund =  2.00
        case ('SB')
          recallabund =  1.00
        case ('TE')
          recallabund =  2.24
        case (' I')
          recallabund =  1.51
        case ('XE')
          recallabund =  2.23
        case ('CS')
          recallabund =  1.12
        case ('BA')
          recallabund =  2.13
        case ('LA')
          recallabund =  1.22
        case ('CE')
          recallabund =  1.55
        case ('PR')
          recallabund =  0.71
        case ('ND')
          recallabund =  1.50
        case ('PM')
          recallabund = -7.96
        case ('SM')
          recallabund =  1.00
        case ('EU')
          recallabund =  0.51
        case ('GD')
          recallabund =  1.12
        case ('TB')
          recallabund = -0.10
        case ('DY')
          recallabund =  1.10
        case ('HO')
          recallabund =  0.26
        case ('ER')
          recallabund =  0.93
        case ('TM')
          recallabund =  0.00
        case ('YB')
          recallabund =  1.08
        case ('LU')
          recallabund =  0.76
        case ('HF')
          recallabund =  0.88
        case ('TA')
          recallabund =  0.13
        case (' W')
          recallabund =  1.11
        case ('RE')
          recallabund =  0.27
        case ('OS')
          recallabund =  1.45
        case ('IR')
          recallabund =  1.35
        case ('PT')
          recallabund =  1.80
        case ('AU')
          recallabund =  1.01
        case ('HG')
          recallabund =  1.09
        case ('TL')
          recallabund =  0.90
        case ('PB')
          recallabund =  1.85
        case ('BI')
          recallabund =  0.71
        case ('PO')
          recallabund = -7.96
        case ('AT')
          recallabund = -7.96
        case ('RN')
          recallabund = -7.96
        case ('FR')
          recallabund = -7.96
        case ('RA')
          recallabund = -7.96
        case ('AC')
          recallabund = -7.96
        case ('TH')
          recallabund =  0.12
        case ('PA')
          recallabund = -7.96
        case (' U')
          recallabund = -0.47
        case ('NP')
          recallabund = -7.96
        case ('PU')
          recallabund = -7.96
        case ('AM')
          recallabund = -7.96
        case ('CM')
          recallabund = -7.96
        case ('BK')
          recallabund = -7.96
        case ('CF')
          recallabund = -7.96
        case ('ES')
          recallabund = -7.96
        case default
          urou = 'recallabund'
          umsg = 'Atom '//atm//' not found in abundance list'
          call aborted
          return
      end select

      ! Convert to N/NH ratio
      recallabund = 1d1**(recallabund - 12d0)

      return

      end function recallabund

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return hard-coded abundance of an atom specified by index in
      !! N/NH\n
      !!  atm(integer): Index of the atom
      double precision function recallabund_ind(atm)

      ! I/O

      integer, intent(in):: atm

      ! Local

      ! Asplund (2009) abundances
      real, dimension(99), parameter:: aspl2009 = &
           (/ 12.00, 10.99, 1.16, 1.15, 2.6, 8.39, 8.00, 8.66, 4.56, &
              8.09, 6.33 , 7.58, 6.47, 7.55, 5.45, 7.21, 5.50, 6.56, &
              5.12, 6.36, 3.10, 4.99, 4.00, 5.67, 5.39, 7.44, 4.92, &
              6.25, 4.21, 4.60, 2.88, 3.41, 2.37, 3.35, 2.63, 3.23, &
              2.60, 2.90, 2.24, 2.60, 1.42, 1.92, -7.96, 1.84, 1.12, &
              1.69, 0.94, 1.86, 1.66, 2.00, 1.00, 2.24, 1.51, 2.23, &
              1.12, 2.13, 1.22, 1.55, 0.71, 1.50, -7.96, 1.00, 0.51, &
              1.12, -0.10, 1.10, 0.26, 0.93, 0.00, 1.08, 0.76, 0.88, &
              0.13, 1.11, 0.27, 1.45, 1.35, 1.80, 1.01, 1.09, 0.90, &
              1.85, 0.71, -7.96, -7.96, -7.96, -7.96, -7.96, -7.96, &
              0.12, -7.96, -0.47, -7.96, -7.96, -7.96, -7.96, -7.96, &
              -7.96, -7.96 /)

      ! Check if the index is out of bounds
      if (atm.lt.1.or.atm.gt.99) then

        ! Call error
        write(umsg,'(A,i2,A)') 'Index ',atm, &
                               ' not found in abundance list'
        urou = 'recallabund_ind'
        call aborted

        ! Dummy value
        recallabund_ind = 0d0

        ! Return
        return

      end if ! Invalid index

      ! Get abundance in N/NH
      recallabund_ind = 1d1**(dble(aspl2009(atm)) - 12d0)

      return

      end function recallabund_ind

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return hard-coded mass of an atom specified by index in AMU\n
      !!   atm(integer): Index of the atom
      double precision function recallmass_ind(atm)

      ! I/O

      integer, intent(in):: atm

      ! Local

      ! List of masses
      real, dimension(99), parameter:: mass = &
           (/ 1.008,4.003,6.939,9.013,10.810,12.010,14.010,16.000, &
              19.000,20.180,22.990,24.310,26.980,28.090,30.980, &
              32.070,35.450,39.950,39.100,40.080,44.960,47.900, &
              50.940,52.000,54.940,55.850,58.940,58.710,63.550, &
              65.370,69.720,72.600,74.920,78.960,79.910,83.800, &
              85.480,87.630,88.910,91.220,92.910,95.950,99.000, &
              101.100,102.900,106.400,107.900,112.400,114.800, &
              118.700,121.800,127.600,126.900,131.300,132.900, &
              137.400,138.900,140.100,140.900,144.300,147.000, &
              150.400,152.000,157.300,158.900,162.500,164.900, &
              167.300,168.900,173.000,175.000,178.500,181.000, &
              183.900,186.300,190.200,192.200,195.100,197.000, &
              200.600,204.400,207.200,209.000,210.000,211.000, &
              222.000,223.000,226.100,227.100,232.000,231.000, &
              238.000,237.000,244.000,243.000,247.000,247.000, &
              251.000,254.000 /)

      ! Check if the index is out of bounds
      if (atm.lt.1.or.atm.gt.99) then

        ! Call error
        write(umsg,'(A,i2,A)') 'Index ',atm, &
                               ' not found in mass list'
        urou = 'recallmass_ind'
        call aborted

        ! Dummy value
        recallmass_ind = -1d0

        ! Return
        return

      end if ! Invalid index

      ! Get mass in AMU
      recallmass_ind = dble(mass(atm))

      return

      end function recallmass_ind

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return the index of an atom specified by name\n
      !!  atm(character(:)): ID name of the atom
      integer function atom_char2index(atm)

      ! I/O

      character(len=2), intent(in):: atm


      ! Select return value given the name
      select case (atm)
        case (' H')
          atom_char2index = 1
        case ('HE')
          atom_char2index = 2
        case ('LI')
          atom_char2index = 3
        case ('BE')
          atom_char2index = 4
        case (' B')
          atom_char2index = 5
        case (' C')
          atom_char2index = 6
        case (' N')
          atom_char2index = 7
        case (' O')
          atom_char2index = 8
        case (' F')
          atom_char2index = 9
        case ('NE')
          atom_char2index = 10
        case ('NA')
          atom_char2index = 11
        case ('MG')
          atom_char2index = 12
        case ('AL')
          atom_char2index = 13
        case ('SI')
          atom_char2index = 14
        case (' P')
          atom_char2index = 15
        case (' S')
          atom_char2index = 16
        case ('CL')
          atom_char2index = 17
        case ('AR')
          atom_char2index = 18
        case (' K')
          atom_char2index = 19
        case ('CA')
          atom_char2index = 20
        case ('SC')
          atom_char2index = 21
        case ('TI')
          atom_char2index = 22
        case (' V')
          atom_char2index = 23
        case ('CR')
          atom_char2index = 24
        case ('MN')
          atom_char2index = 25
        case ('FE')
          atom_char2index = 26
        case ('CO')
          atom_char2index = 27
        case ('NI')
          atom_char2index = 28
        case ('CU')
          atom_char2index = 29
        case ('ZN')
          atom_char2index = 30
        case ('GA')
          atom_char2index = 31
        case ('GE')
          atom_char2index = 32
        case ('AS')
          atom_char2index = 33
        case ('SE')
          atom_char2index = 34
        case ('BR')
          atom_char2index = 35
        case ('KR')
          atom_char2index = 36
        case ('RB')
          atom_char2index = 37
        case ('SR')
          atom_char2index = 38
        case (' Y')
          atom_char2index = 39
        case ('ZR')
          atom_char2index = 40
        case ('NB')
          atom_char2index = 41
        case ('MO')
          atom_char2index = 42
        case ('TC')
          atom_char2index = 43
        case ('RU')
          atom_char2index = 44
        case ('RH')
          atom_char2index = 45
        case ('PD')
          atom_char2index = 46
        case ('AG')
          atom_char2index = 47
        case ('CD')
          atom_char2index = 48
        case ('IN')
          atom_char2index = 49
        case ('SN')
          atom_char2index = 50
        case ('SB')
          atom_char2index = 51
        case ('TE')
          atom_char2index = 52
        case (' I')
          atom_char2index = 53
        case ('XE')
          atom_char2index = 54
        case ('CS')
          atom_char2index = 55
        case ('BA')
          atom_char2index = 56
        case ('LA')
          atom_char2index = 57
        case ('CE')
          atom_char2index = 58
        case ('PR')
          atom_char2index = 59
        case ('ND')
          atom_char2index = 60
        case ('PM')
          atom_char2index = 61
        case ('SM')
          atom_char2index = 62
        case ('EU')
          atom_char2index = 63
        case ('GD')
          atom_char2index = 64
        case ('TB')
          atom_char2index = 65
        case ('DY')
          atom_char2index = 66
        case ('HO')
          atom_char2index = 67
        case ('ER')
          atom_char2index = 68
        case ('TM')
          atom_char2index = 69
        case ('YB')
          atom_char2index = 70
        case ('LU')
          atom_char2index = 71
        case ('HF')
          atom_char2index = 72
        case ('TA')
          atom_char2index = 73
        case (' W')
          atom_char2index = 74
        case ('RE')
          atom_char2index = 75
        case ('OS')
          atom_char2index = 76
        case ('IR')
          atom_char2index = 77
        case ('PT')
          atom_char2index = 78
        case ('AU')
          atom_char2index = 79
        case ('HG')
          atom_char2index = 80
        case ('TL')
          atom_char2index = 81
        case ('PB')
          atom_char2index = 82
        case ('BI')
          atom_char2index = 83
        case ('PO')
          atom_char2index = 84
        case ('AT')
          atom_char2index = 85
        case ('RN')
          atom_char2index = 86
        case ('FR')
          atom_char2index = 87
        case ('RA')
          atom_char2index = 88
        case ('AC')
          atom_char2index = 89
        case ('TH')
          atom_char2index = 90
        case ('PA')
          atom_char2index = 91
        case (' U')
          atom_char2index = 92
        case ('NP')
          atom_char2index = 93
        case ('PU')
          atom_char2index = 94
        case ('AM')
          atom_char2index = 95
        case ('CM')
          atom_char2index = 96
        case ('BK')
          atom_char2index = 97
        case ('CF')
          atom_char2index = 98
        case ('ES')
          atom_char2index = 99
        case default
          urou = 'atom_char2index'
          umsg = 'Atom '//atm//' not found in index list'
          call aborted
          return
      end select

      ! Return
      return

      end function atom_char2index

!#####################################################################
!#####################################################################
!#####################################################################

      !> Return name of an atom specified by index\n
      !!  atm(integer): Index of the atom
      character(len=2) function atom_index2char(atm)

      ! I/O
      integer, intent(in):: atm

      ! Local

      ! List of names
      character(len=198), parameter:: symbols = &
          ' HHELIBE B C N O FNENAMGALSI P SCLAR KCASCTI VCRMNFE'// &
          'CONICUZNGAGEASSEBRKRRBSR YZRNBMOTCRURHPDAGCDINSNSBTE'// &
          ' IXECSBALACEPRNDPMSMEUGDTBDYHOERTMYBLUHFTA WREOSIRPT'// &
          'AUHGTLPBBIPOATRNFRRAACTHPA UNPPUAMCMBKCFES'

      ! Check if the index is out of bounds
      if (atm.lt.1.or.atm.gt.99) then

        ! Call error
        urou = 'atom_index2char'
        write(umsg,'(A,i2,A)') 'Index ',atm, &
                               ' not found in atom symbol list'
        call aborted

        ! Return
        return

      end if ! Invalid index

      ! Get atom ID name
      atom_index2char = symbols(atm*2-1:atm*2)

      end function atom_index2char

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get partition function, ionization energy, and number of
      !! stages of the atom specified by name at a given temperature\n
      !!  Ele(character(:)): ID name of the atom\n
      !!      nstg(integer): Number of stages\n
      !!    pf(double(:,:)): Partition function\n
      !!      Ei(double(:)): Ionization energy\n
      !!   Atmo(Atmo_class): Structure with PF data\n
      !!          T(double): Temperature
      subroutine getpf_T(Ele,nstg,pf,Ei,Atmo,T)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      character(len=2), intent(in):: Ele
      integer, intent(out):: nstg
      double precision, intent(in):: T
      double precision, dimension(:), allocatable, intent(inout):: Ei
      double precision, dimension(:), allocatable, intent(inout):: pf

      ! Local

      logical:: interp

      integer:: istg,ii,iele,ind

      double precision:: x0,x1,y0,y1,dx


      ! Allocate the partition function and ionization energies
      ! They are not added to memory count because they live
      ! only in the chemeq() subroutine
      iele = atom_char2index(Ele)
      nstg = Atmo%ele(iele)%nstg
      allocate(pf(nstg))
      allocate(Ei(nstg))

      ! If temperature below lower limit
      if (T.le.Atmo%pT(1)) then

        ! Take constant in first data position
        interp = .False.
        ind = 1

      ! If temperature above upper limit
      else if (T.ge.(Atmo%pT(Atmo%NT) - 1d-6)) then

        ! Take constant in last data position
        interp = .False.
        ind = Atmo%NT

      ! If within boundaries
      else

        ! Run over all tabulated temperatures
        do ii=1,Atmo%NT-1

          ! If found an exact value
          if (abs(T - Atmo%pT(ii)).le.1d-6) then

            ! Take constant in found position
            interp = .False.
            ind = ii
            exit

          ! If temperature between two tabulated temperatures
          else if (T.gt.Atmo%pT(ii).and. &
                   T.lt.Atmo%pT(ii+1)) then

            ! Linear interpolation between the two points
            interp = .True.
            ind = ii
            dx = 1d0/(Atmo%pT(ii+1) - Atmo%pT(ii))
            exit

          end if ! Exact or found between

        end do ! Tabulated temperatures

      end if ! Out/in boundaries

      ! For each stage
      do istg=1,nstg

        ! If interpolating
        if (interp) then

          ! Get x and y axis data for the two involved points
          y1 = Atmo%ele(iele)%pf(ind + 1,istg)
          y0 = Atmo%ele(iele)%pf(ind,istg)
          x1 = Atmo%pT(ind + 1)
          x0 = Atmo%pT(ind)

          ! Linear interpolation (dx was pre-calculated above)
          pf(istg) = ((y1 - y0)*T + y0*x1 - y1*x0)*dx

        ! If not interpolating
        else

          ! Just take the value at the suitable location
          pf(istg) = Atmo%ele(iele)%pf(ind,istg)

        end if ! Interpolating

        ! Ionization energy
        Ei(istg) = Atmo%ele(iele)%Ei(istg)

      end do ! Ionization stages

      return

      end subroutine getpf_T

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get partition function, ionization energy, and number of
      !! stages of the atom specified by name at all heights\n
      !!  Ele(character(:)): ID name of the atom\n
      !!      nstg(integer): Number of stages\n
      !!    pf(double(:,:)): Partition function\n
      !!      Ei(double(:)): Ionization energy\n
      !!   Atmo(Atmo_class): Structure with atmospheric data
      subroutine getpf(Ele,nstg,pf,Ei,Atmo)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      character(len=2), intent(in):: Ele
      integer, intent(out):: nstg
      double precision, dimension(:), allocatable, intent(inout):: Ei
      double precision, dimension(:,:), &
                        allocatable, intent(inout):: pf

      ! Local

      logical:: interp

      integer:: istg,ii,iele,iz,ind

      double precision:: x0,x1,y0,y1,dx


      ! Allocate the partition function and ionization energies
      ! They are not added to memory count because they live
      ! only in the chemeq() subroutine
      iele = atom_char2index(Ele)
      nstg = Atmo%ele(iele)%nstg
      allocate(pf(nstg,NZ))
      allocate(Ei(nstg))

      ! For each height
      do iz=1,NZ

        ! If temperature below lower limit
        if (Atmo%T(iz).le.Atmo%pT(1)) then

          ! Take constant in first data position
          interp = .False.
          ind = 1

        ! If temperature above upper limit
        else if (Atmo%T(iz).ge.(Atmo%pT(Atmo%NT) - 1d-6)) then

          ! Take constant in last data position
          interp = .False.
          ind = Atmo%NT

        ! If within boundaries
        else

          ! Run over all tabulated temperatures
          do ii=1,Atmo%NT-1

            ! If found an exact value
            if (abs(Atmo%T(iz) - Atmo%pT(ii)).le.1d-6) then

              ! Take constant in found position
              interp = .False.
              ind = ii
              exit

            ! If temperature between two tabulated temperatures
            else if (Atmo%T(iz).gt.Atmo%pT(ii).and. &
                     Atmo%T(iz).lt.Atmo%pT(ii+1)) then

              ! Linear interpolation between the two points
              interp = .True.
              ind = ii
              dx = 1d0/(Atmo%pT(ii+1) - Atmo%pT(ii))
              exit

            end if ! Exact or found between

          end do ! Tabulated temperatures

        end if ! Out/in boundaries

        ! For each stage
        do istg=1,nstg

          ! If interpolating
          if (interp) then

            ! Get x and y axis data for the two involved points
            y1 = Atmo%ele(iele)%pf(ind + 1,istg)
            y0 = Atmo%ele(iele)%pf(ind,istg)
            x1 = Atmo%pT(ind + 1)
            x0 = Atmo%pT(ind)

            ! Linear interpolation (dx was pre-calculated above)
            pf(istg,iz) = ((y1 - y0)*Atmo%T(iz) + &
                           y0*x1 - y1*x0)*dx

          ! If not interpolating
          else

            ! Just take the value at the suitable location
            pf(istg,iz) = Atmo%ele(iele)%pf(ind,istg)

          end if ! Interpolating

          ! Ionization energy
          Ei(istg) = Atmo%ele(iele)%Ei(istg)

        end do ! Ionization stages
      end do ! Heights

      return

      end subroutine getpf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get ionization fraction for the specified stages at a given
      !! temperature and for a given electron density\n
      !!   nstg(integer): Number of stages to consider\n
      !!   pf(double(:)): Partition function\n
      !!   Ei(double(:)): Ionization energy\n
      !!       T(double): Temperature\n
      !!      ne(double): Electron density\n
      !!    ion(integer): Requested ion\n
      !!  frc(double(:)): Ionization fraction
      subroutine getfrc(nstg,pf,Ei,T,ne,ion,frc)

      ! I/O

      integer, intent(in):: nstg,ion
      double precision, intent(in):: T,ne
      double precision, dimension(:), intent(in):: Ei,pf
      double precision, dimension(:), intent(inout):: frc

      ! Local

      integer:: istg

      double precision:: U0, U1, C0, C1, S, ikT, arg, exu
      double precision, dimension(nstg):: frcl


      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      C1 = 1d0/(.5d6*ne*((C0/T)**(1.5d0)))
      ikT = fktoJ/kb/T

      ! Initializations
      S = 1d0
      frcl(1) = 1d0

      ! Neutral partition function
      U0 = pf(1)

      ! For each stage
      do istg=2,nstg

        ! Stage partition function
        U1 = pf(istg)

        ! Argument of the exponential
        arg = U1 - U0 - Ei(istg-1)*ikT

        ! Calculate exponential for the numerator of the population
        ! solution
        if (arg.lt.0d0) then
          arg = -arg
          exu = diexp(arg)
        else
          exu = ddexp(arg)
        end if

        ! Numerator of the ionization fraction for this stage
        frcl(istg) = frcl(istg-1)*C1*exu

        ! Avoid very small fractions
        if (frcl(istg).lt.TINYFRC) frcl(istg) = 0d0

        ! Accumulate fractions
        S = S + frcl(istg)

        ! Shift the partition function
        U0 = U1

      end do ! Stages

      ! Compute the fractions
      frcl = frcl/S

      ! If called from chemic subroutine
      if (ion.lt.0) then

        ! We only want the first two
        frc = frcl(1:2)

      ! If calling from redo_ne subroutine
      else if (ion.eq.0) then

        ! We want them all
        frc = frcl

      ! If calling from Kurucz subroutine
      else

        ! We only want one of them
        frc = frcl(ion)

      end if

      return

      end subroutine getfrc

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get ionization fraction derivative at a given temperature and
      !! for a given electron density\n
      !!    nstg(integer): Number of stages to consider\n
      !!    pf(double(:)): Partition function\n
      !!    Ei(double(:)): Ionization energy\n
      !!        T(double): Temperature\n
      !!       ne(double): Electron density\n
      !!   frc(double(:)): Ionization fraction\n
      !!  dfrc(double(:)): Ionization fraction derivative
      subroutine getdfrc(nstg,pf,Ei,T,ne,frc,dfrc)

      ! I/O

      integer, intent(in):: nstg
      double precision, intent(in):: T,ne
      double precision, dimension(:), intent(in):: Ei,pf
      double precision, dimension(:), intent(out):: frc,dfrc

      ! Local

      integer:: istg

      double precision:: U0,U1,C0,C1,S,S2,ikT,arg,exu


      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      C1 = 1d0/(.5d6*ne*((C0/T)**(1.5d0)))
      ikT = fktoJ/kb/T

      ! Initializations
      S = 1d0
      S2 = 0d0
      frc(1) = 1d0
      dfrc(1) = 0d0

      ! Neutral partition function
      U0 = pf(1)

      ! For each stage
      do istg=2,nstg

        ! Stage partition function
        U1 = pf(istg)

        ! Argument of the exponential
        arg = U1 - U0 - Ei(istg-1)*ikT

        ! Exponential for the numerator of the population solution
        if (arg.lt.0d0) then
          arg = -arg
          exu = diexp(arg)
        else
          exu = ddexp(arg)
        end if

        ! Numerator of the fraction for this stage
        frc(istg) = frc(istg-1)*C1*exu

        ! Avoid very small fractions
        if (frc(istg).lt.TINYFRC) frc(istg) = 0d0

        ! For derivative
        dfrc(istg) = -dble(istg-1)*frc(istg)/ne

        ! Accumulate fractions
        S = S + frc(istg)
        S2 = S2 + dfrc(istg)

        ! Shift the partition function
        U0 = U1

      end do ! Stages

      ! For each stage
      do istg=1,nstg

        ! Compute the fractions
        frc(istg) = frc(istg)/S
        dfrc(istg) = (dfrc(istg) - frc(istg)*S2)/S

      end do ! Stages

      return

      end subroutine getdfrc

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get pressure and mass density for known number densities
      !! at every height\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      subroutine eqstate_known(Atmo,Atom,Atomb,nlte,depar,atoms)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(Atmo_class), intent(inout):: Atmo
      type(catm_class), dimension(:), intent(in):: atoms
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      integer:: iz


      ! For each height
      do iz=1,nz

        ! Compute gas and electron pressure for known number
        ! densities
        call partial_press_known(Atmo,iz,Atom,Atomb,atoms,nlte,depar)

      end do ! Heights

      end subroutine eqstate_known

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get pressure and mass density for known gas pressure at every
      !! height\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      subroutine eqstate_gas(Atmo,Atom,Atomb,nlte,depar,atoms)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(Atmo_class), intent(inout):: Atmo
      type(catm_class), dimension(:), intent(in):: atoms
      integer, dimension(:), intent(in):: nlte, depar

      ! Local

      integer:: iz

      ! First guess for electron pressure
      Atmo%Pe = 0.3d0*Atmo%Pg
      Atmo%ne = Atmo%Pe*1d-7/(Atmo%T*kb)

      ! For each height
      do iz=1,nz

        ! Call partial pressure calculator
        call partial_press_Pg(Atmo,iz,Atom,Atomb,atoms,nlte,depar)

        ! Set atmo if LTE H
        if (nlte(1).eq.0.and.depar(1).eq.0) &
          call set_Hdensities(Atmo,iz)

      end do ! Heights

      end subroutine eqstate_gas

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get pressure and mass density for known electron pressure at
      !! every height\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      subroutine eqstate_ele(Atmo,Atom,Atomb,nlte,depar,atoms)

      ! I/O

      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(Atmo_class), intent(inout):: Atmo
      type(catm_class), dimension(:), intent(in):: atoms
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      integer:: iz

      ! For each height
      do iz=1,nz

        ! Compute partial pressures
        call partial_press_Pe(Atmo,iz,Atom,Atomb,atoms,nlte,depar)

        ! Set atmo if LTE H
        if (nlte(1).eq.0.and.depar(1).eq.0) &
          call set_Hdensities(Atmo,iz)

      end do

      end subroutine eqstate_ele

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute partial pressures of ions from known ionization
      !! data\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!           iz(integer): Height index\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      subroutine partial_press_known(Atmo,iz,Atom,Atomb,atoms, &
                                     nlte,depar)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(catm_class), dimension(:), intent(in):: atoms
      integer, intent(in):: iz
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      integer:: im,iatom,ia,ii
      double precision:: t,Pe,Pa,PHtot
      double precision:: a,b,c,d,e
      double precision:: f1,f2,f3,f4,f5,fe
      double precision:: g1,g2,g3,g4,g5
      double precision:: c1,c2,c3,c6,c7

      ! Variables dependent on parameters
      double precision, dimension(nmol):: cmol


      ! Temperature in eV and gas known electron pressure
      t = 5040d0/Atmo%T(iz)
      Pe = Atmo%Pe(iz)

      ! Get molecular data
      do im=1,nmol
        cmol(im) = moldata_ind(t,im)
      end do

      ! H2+ and H2 dissociation coefficient inverses, multiplied by
      ! electron pressure
      g4 = Pe*10d0**cmol(1)
      g5 = Pe*10d0**cmol(2)

      !
      ! Hydrogen
      !

      ! P(H+)/P(H)
      g2 = Atmo%nH(iz,6)/sum(Atmo%nH(iz,1:5))

      ! P(H-)/P(H)
      g3 = 1d0/fsaha(t,0.754d0,1d0,atoms(1)%pf(1,iz),Pe)

      ! Initialize P(atoms)/P(H)
      g1 = 0d0

      ! For every atom taken into account
      do iatom=2,natom

        ! If full LTE
        if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

          ! Ionization fraction first ion to neutral
          a = fsaha(t,atoms(iatom)%Eion(1),atoms(iatom)%pf(1,iz), &
                    atoms(iatom)%pf(2,iz),Pe)

          ! Ionization fraction second ion to first ion
          b = fsaha(t,atoms(iatom)%Eion(2),atoms(iatom)%pf(2,iz), &
                    atoms(iatom)%pf(3,iz),Pe)

        ! If NLTE
        else

          ! Departure coefficients
          if (depar(iatom).ne.0) then

            ! Atom index
            ia = abs(depar(iatom))

            ! Active
            if (depar(iatom).gt.0) then

              ! Get fraction
              call metal_fraction(Atom(ia),Atmo,iz,a,b,.True.)

            ! Passive
            else if (depar(iatom).lt.0) then

              ! Get fraction
              call metal_fraction(Atomb(ia),Atmo,iz,a,b,.True.)

            end if ! Active/Passive

          ! Populations
          else

            ! Atom index
            ia = abs(nlte(iatom))

            ! Active
            if (nlte(iatom).gt.0) then

              ! Get fraction
              call metal_fraction(Atom(ia),Atmo,iz,a,b,.False.)

            ! Passive
            else if (nlte(iatom).lt.0) then

              ! Get fraction
              call metal_fraction(Atomb(ia),Atmo,iz,a,b,.False.)

            end if ! Active/Passive
          end if ! Departure/populations

          ! If invalid fractions
          if (a.lt.0d0) then

            !
            ! Do LTE
            !

            ! Ionization fraction first ion to neutral
            a = fsaha(t,atoms(iatom)%Eion(1),atoms(iatom)%pf(1,iz), &
                      atoms(iatom)%pf(2,iz),Pe)

            ! Ionization fraction second ion to first ion
            b = fsaha(t,atoms(iatom)%Eion(2),atoms(iatom)%pf(2,iz), &
                      atoms(iatom)%pf(3,iz),Pe)

          end if ! Invalid fraction
        end if ! LTE/NLTE

        ! Normalization
        c = 1d0 + a*(1d0 + b)

        ! Partial pressure for the atom
        Pa = atoms(iatom)%abun/c

        ! Add to total atomic pressure
        g1 = g1 + Pa*a*(1d0 + 2d0*b)

      end do ! Atoms

      ! H + H+ + H- pressure/p(H)
      a = 1d0 + g2 + g3

      ! Compute magical coefficients
      c = g5
      d = g2 - g3
      e = g4*g2/g5
      b = 2d0*(1d0 + e)
      c1 = c*b*b + a*d*b - e*a*a
      c2 = 2d0*a*e - d*b + a*b*g1
      c3 = -(e + b*g1)
      f1 = 0.5d0*c2/c1
      f1 = sign(1d0,c1)*sqrt(f1*f1 - c3/c1) - f1
      f5 = (1d0 - a*f1)/b
      f4 = e*f5
      f3 = g3*f1
      f2 = g2*f1
      fe = f2 - f3 + f4 + g1

      ! H pressure
      PHtot = Pe/fe

      ! If whatever is f5 is too small
      if (f5.le.1d-4) then

        ! Whatever this is
        c6 = g5*f1*f1/Pe
        c7 = g1 + f2 - f3

        ! No idea what this is for
        do ii=1,5
          f5 = PHtot*c6
          f4 = e*f5
          fe = c7 + f4
          PHtot = Pe/fe
        end do

      end if ! Dimension f5

      ! Solve for gas pressure
      Atmo%Pg(iz) = Pe*(1d0 + (f1 + f2 + f3 + f4 + f5 + 0.1014d0)/fe)

      end subroutine partial_press_known

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute partial pressures of ions from known gas pressure\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!           iz(integer): Height index\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      subroutine partial_press_Pg(Atmo,iz,Atom,Atomb,atoms,nlte,depar)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(catm_class), dimension(:), intent(in):: atoms
      integer, intent(in):: iz
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      integer:: im,iatom,ia,ii,iter

      double precision:: t,Pe,Pe0,Pa,PHtot,Pg
      double precision:: a,b,c,d,e
      double precision:: f1,f2,f3,f4,f5,fe
      double precision:: g1,g2,g3,g4,g5
      double precision:: c1,c2,c3,c6,c7
      double precision:: ikbcgsT,diff
      double precision, dimension(6):: NH

      ! Variables dependent on parameters
      double precision, dimension(nmol):: cmol


      ! Constants
      ikbcgsT = 1d-7/kb/Atmo%T(iz)

      ! Input quantities
      t = 5040d0/Atmo%T(iz)
      NH = Atmo%nH(iz,:)
      Pg = Atmo%Pg(iz)
      Pe = Atmo%Pe(iz)
      Pe0 = Pe

      ! Get molecular data
      do im=1,nmol
        cmol(im) = moldata_ind(t,im)
      end do

      ! Initialize electron pressure change
      diff = 1d0

      ! Iterate
      do iter=1,maxitereq

        ! New pressure average of previous and current
        Pe0 = (Pe + Pe0)*0.5d0
        Pe = Pe0

        ! Other molecular quantities
        g4 = Pe*10d0**cmol(1)
        g5 = Pe*10d0**cmol(2)

        ! Hydrogen

        ! LTE
        if (nlte(1).eq.0.and.depar(1).eq.0) then

          ! P(H+)/P(H)
          g2 = fsaha(t,atoms(1)%Eion(1),atoms(1)%pf(1,iz), &
                     atoms(1)%pf(2,iz),Pe)

        ! Departure coefficients
        else if (depar(1).ne.0) then

          ! Atom index
          ia = abs(depar(1))

          ! Active
          if (depar(1).gt.0) then

            ! Get fraction
            call H_fraction(Atom(ia),Atmo,iz,g2,.True.)

          ! Passive
          else if (depar(1).lt.0) then

            ! Get fraction
            call H_fraction(Atomb(ia),Atmo,iz,g2,.True.)

          end if ! Active/Passive

        ! NLTE populations
        else

          ! P(H+)/P(H)
          g2 = NH(6)/sum(NH(1:5))

        end if

        ! P(H-)/P(H)
        g3 = 1d0/fsaha(t,0.754d0,1d0,atoms(1)%pf(1,iz),Pe)

        ! Initialize P(atoms)/P(H)
        g1 = 0d0

        ! For every atom taken into account
        do iatom=2,natom

          ! If full LTE
          if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

            ! Ionization fraction first ion to neutral
            a = fsaha(t,atoms(iatom)%Eion(1),atoms(iatom)%pf(1,iz), &
                      atoms(iatom)%pf(2,iz),Pe)

            ! Ionization fraction second ion to first ion
            b = fsaha(t,atoms(iatom)%Eion(2),atoms(iatom)%pf(2,iz), &
                      atoms(iatom)%pf(3,iz),Pe)

          ! If NLTE
          else

            ! If departure coefficients
            if (depar(iatom).ne.0) then

              ! Atom index
              ia = abs(depar(iatom))

              ! Active
              if (depar(iatom).gt.0) then

                ! Get fraction
                call metal_fraction(Atom(ia),Atmo,iz,a,b,.True.)

              ! Passive
              else if (depar(iatom).lt.0) then

                ! Get fraction
                call metal_fraction(Atomb(ia),Atmo,iz,a,b,.True.)

              end if ! Active/Passive

            ! If populations
            else

              ! Atom index
              ia = abs(nlte(iatom))

              ! Active
              if (nlte(iatom).gt.0) then

                ! Get fraction
                call metal_fraction(Atom(ia),Atmo,iz,a,b,.False.)

              ! Passive
              else if (nlte(iatom).lt.0) then

                ! Get fraction
                call metal_fraction(Atomb(ia),Atmo,iz,a,b,.False.)

              end if ! Active/Passive
            end if ! Departure or populations

            ! If invalid fractions
            if (a.lt.0d0) then

              !
              ! Do LTE
              !

              ! Fraction first ion to neutral
              a = fsaha(t,atoms(iatom)%Eion(1), &
                        atoms(iatom)%pf(1,iz), &
                        atoms(iatom)%pf(2,iz),Pe)

              ! Fraction second ion to first ion
              b = fsaha(t,atoms(iatom)%Eion(2), &
                        atoms(iatom)%pf(2,iz), &
                        atoms(iatom)%pf(3,iz),Pe)

            end if ! Invalid fraction
          end if ! LTE/NLTE

          ! Normalization
          c = 1d0 + a*(1d0 + b)

          ! Partial pressure for the atom
          Pa = atoms(iatom)%abun/c

          ! Add to total atomic pressure
          g1 = g1 + Pa*a*(1d0 + 2d0*b)

        end do ! Atoms

        ! H + H+ + H- pressure
        a = 1d0 + g2 + g3

        ! Compute magical coefficients
        c = g5
        d = g2 - g3
        e = g4*g2/g5
        b = 2d0*(1d0 + e)
        c1 = c*b*b + a*d*b - e*a*a
        c2 = 2d0*a*e - d*b + a*b*g1
        c3 = -(e + b*g1)
        f1 = 0.5d0*c2/c1
        f1 = sign(1d0,c1)*sqrt(f1*f1 - c3/c1) - f1
        f5 = (1d0 - a*f1)/b
        f4 = e*f5
        f3 = g3*f1
        f2 = g2*f1
        fe = f2 - f3 + f4 + g1

        ! H pressure
        PHtot = Pe/fe

        ! If whatever is f5 is too small
        if (f5.le.1d-4) then

          ! Whatever this is
          c6 = g5*f1*f1/Pe
          c7 = g1 + f2 - f3

          ! No idea what this is for
          do ii=1,5
            f5 = PHtot*c6
            f4 = e*f5
            fe = c7 + f4
            PHtot = Pe/fe
          end do

        end if ! Dimension f5

        ! Solve for electron pressure
        Pe = Pg/(1d0 + (f1 + f2 + f3 + f4 + f5 + 0.1014d0)/fe)

        ! Compute difference
        diff = abs(Pe0-Pe)/Pe0

        ! Scale densities
        PHtot = PHtot*ikbcgsT
        NH(6) = f2*PHtot

        !
        ! Neutral hydrogen
        !

        ! If LTE
        if (nlte(1).eq.0.and.depar(1).eq.0) then

          ! Everything in ground level
          NH(1) = f1*PHtot
          NH(2:5) = 0d0

        ! If NLTE
        else

          ! Adjust ionization without changing relative population
          NH(1:5) = NH(1:5)*f1*PHtot/sum(NH(1:5))

        end if

        ! Check convergence
        if (diff.lt.eps) exit

      end do ! Iterations

      ! Electron pressure and density
      Atmo%Pe(iz) = Pe
      Atmo%ne(iz) = Atmo%Pe(iz)*ikbcgsT

      ! Atmospheric hydrogen
      Atmo%nH(iz,:) = NH
      Atmo%nHa(iz) = sum(NH)
      Atmo%nHt(iz) = Phtot

      end subroutine partial_press_Pg

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute partial pressures of ions from known electron
      !! pressure\n
      !!      Atmo(Atmo_class): Structure with atmospheric data\n
      !!           iz(integer): Height index\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atins\n
      !!     atoms(catm_class): Array of structures with
      !!                        ionization and partition function
      !!                        data
      !!      nlte(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      !!     depar(integer(:)): Array with indication about the type
      !!                        of ionization balance\n
      subroutine partial_press_Pe(Atmo,iz,Atom,Atomb,atoms,nlte,depar)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), &
                        allocatable, intent(in):: Atomb
      type(catm_class), dimension(:), intent(in):: atoms
      integer, intent(in):: iz
      integer, dimension(:), intent(in):: nlte,depar

      ! Local

      integer:: im,iatom,ia,ii
      double precision:: t,Pe,Pg,Pa,PHtot
      double precision:: a,b,c,d,e
      double precision:: f1,f2,f3,f4,f5,fe
      double precision:: g1,g2,g3,g4,g5
      double precision:: c1,c2,c3,c6,c7
      double precision:: ikbcgsT
      double precision, dimension(6):: NH

      ! Variables dependent on parameters
      double precision, dimension(nmol):: cmol


      ! Constants
      ikbcgsT = 1d-7/kb/Atmo%T(iz)

      ! Input quantities
      t = 5040d0/Atmo%T(iz)
      NH = Atmo%nH(iz,:)
      Pe = Atmo%Pe(iz)

      ! Get molecular data
      do im=1,nmol
        cmol(im) = moldata_ind(t,im)
      end do

      ! Other molecular quantities
      g4 = Pe*10d0**cmol(1)
      g5 = Pe*10d0**cmol(2)

      !
      ! Hydrogen
      !

      ! LTE
      if (nlte(1).eq.0.and.depar(1).eq.0) then

        ! P(H+)/P(H)
        g2 = fsaha(t,atoms(1)%Eion(1),atoms(1)%pf(1,iz), &
                   atoms(1)%pf(2,iz),Pe)

      ! Departure coefficients
      else if (depar(1).ne.0) then

        ! Atom index
        ia = abs(depar(1))

        ! Active
        if (depar(1).gt.0) then

          ! Get fraction
          call H_fraction(Atom(ia),Atmo,iz,g2,.True.)

        ! Passive
        else if (depar(1).lt.0) then

          ! Get fraction
          call H_fraction(Atomb(ia),Atmo,iz,g2,.True.)

        end if ! Active/Passive

      ! NLTE
      else

        ! P(H+)/P(H)
        g2 = Atmo%nH(iz,6)/sum(Atmo%nH(iz,1:5))

      end if

      ! P(H-)/P(H)
      g3 = 1d0/fsaha(t,0.754d0,1d0,atoms(1)%pf(1,iz),Pe)

      ! Initialize P(atoms)/P(H)
      g1 = 0d0

      ! For every atom taken into account
      do iatom=2,natom

        ! If LTE
        if (nlte(iatom).eq.0.and.depar(iatom).eq.0) then

          ! Fraction first ion to neutral
          a = fsaha(t,atoms(iatom)%Eion(1),atoms(iatom)%pf(1,iz), &
                    atoms(iatom)%pf(2,iz),Pe)

          ! Fraction second ion to neutral
          b = fsaha(t,atoms(iatom)%Eion(2),atoms(iatom)%pf(2,iz), &
                    atoms(iatom)%pf(3,iz),Pe)

        ! If NLTE
        else

          ! If departure coefficients
          if (depar(iatom).ne.0) then

            ! Atom index
            ia = abs(depar(iatom))

            ! Active
            if (depar(iatom).gt.0) then

              ! Get fraction
              call metal_fraction(Atom(ia),Atmo,iz,a,b,.True.)

            ! Passive
            else if (depar(iatom).lt.0) then

              ! Get fraction
              call metal_fraction(Atomb(ia),Atmo,iz,a,b,.True.)

            end if ! Active/Passive

          ! If populations
          else

            ! Atom index
            ia = abs(nlte(iatom))

            ! Active
            if (nlte(iatom).gt.0) then

              ! Get fraction
              call metal_fraction(Atom(ia),Atmo,iz,a,b,.False.)

            ! Passive
            else if (nlte(iatom).lt.0) then

              ! Get fraction
              call metal_fraction(Atomb(ia),Atmo,iz,a,b,.False.)

            end if ! Active/Passive
          end if ! population/departure

          ! If invalid fractions
          if (a.lt.0d0) then

            !
            ! Do LTE
            !

            ! Ionization fraction first ion to neutral
            a = fsaha(t,atoms(iatom)%Eion(1),atoms(iatom)%pf(1,iz), &
                      atoms(iatom)%pf(2,iz),Pe)

            ! Ionization fraction second ion to first ion
            b = fsaha(t,atoms(iatom)%Eion(2),atoms(iatom)%pf(2,iz), &
                      atoms(iatom)%pf(3,iz),Pe)

          end if ! Invalid fraction
        end if ! LTE/NLTE

        ! Normalization
        c = 1d0 + a*(1d0 + b)

        ! Partial pressure for the atom
        Pa = atoms(iatom)%abun/c

        ! Add to total atomic pressure
        g1 = g1 + Pa*a*(1d0 + 2d0*b)

      end do ! Atoms

      ! H + H+ + H- pressure
      a = 1d0 + g2 + g3

      ! Compute magical coefficients
      c = g5
      d = g2 - g3
      e = g4*g2/g5
      b = 2d0*(1d0 + e)
      c1 = c*b*b + a*d*b - e*a*a
      c2 = 2d0*a*e - d*b + a*b*g1
      c3 = -(e + b*g1)
      f1 = 0.5d0*c2/c1
      f1 = sign(1d0,c1)*sqrt(f1*f1 - c3/c1) - f1
      f5 = (1d0 - a*f1)/b
      f4 = e*f5
      f3 = g3*f1
      f2 = g2*f1
      fe = f2 - f3 + f4 + g1

      ! H pressure
      PHtot = Pe/fe

      ! If whatever is f5 is too small
      if (f5.le.1d-4) then

        c6 = g5*f1*f1/Pe
        c7 = g1 + f2 - f3

        ! No idea what this is for
        do ii=1,5
          f5 = Phtot*c6
          f4 = e*f5
          fe = c7 + f4
          PHtot = Pe/fe
        end do

      end if ! Dimension f5

      ! Solve for gas pressure
      Pg = Pe*(1d0 + (f1 + f2 + f3 + f4 + f5 + 0.1014d0)/fe)

      ! Scale densities
      PHtot = PHtot*ikbcgsT

      ! Protons
      NH(6) = f2*PHtot

      !
      ! Neutral hydrogen
      !

      ! If LTE
      if (nlte(1).eq.0) then

        ! Everything in ground level
        NH(1) = f1*PHtot
        NH(2:5) = 0d0

      ! If NLTE
      else

        ! Adjust ionization without changing relative population
        NH(1:5) = NH(1:5)*f1*PHtot/sum(NH(1:5))

      end if

      ! Gas pressure
      Atmo%Pg(iz) = Pg

      ! Atmospheric hydrogen
      Atmo%nH(iz,:) = NH
      Atmo%nHa(iz) = sum(NH)
      Atmo%nHt(iz) = PHtot

      end subroutine partial_press_Pe

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute ion fraction from populations in an atomic model\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iz(integer): Height index\n
      !!        f1(double): Ion II to I\n
      !!        f2(double): Ion III to II\n
      !!    depar(logical): If using departure coefficients
      subroutine metal_fraction(Atom,Atmo,iz,f1,f2,depar)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: depar
      integer, intent(in):: iz
      double precision, intent(out):: f1,f2

      ! Local
      integer:: iterm,iJ,ilevel,istg
      double precision:: f0
      double precision, dimension(:), allocatable:: ltepopu


      ! If not neutral or first stage in the model
      if (minval(Atom%stage).ge.2.or. &
          maxval(Atom%stage).le.1) then

        ! Return invalid values
        f1 = -1d0
        f2 = -1d0
        return

      end if ! No neutral or first stage

      ! Initialize quantities
      f0 = 0d0
      f1 = 0d0
      f2 = 0d0

      ! If departure coefficients
      if (depar) then

        ! Allocate ltepop
        allocate(ltepopu(Atom%nlevel))

        ! Compute LTE
        call LTEiz(Atom,Atmo,iz,ltepopu)

        ! Apply departure coefficients
        ltepopu = ltepopu*Atom%depar(:,iz)

        ! Initialize level
        ilevel = 0

        ! For each level
        do iterm=1,Atom%nmulti
          do iJ=1,Atom%nJ(iterm)

            ! Advance index
            ilevel = ilevel + 1

            ! Stage
            istg = Atom%stage(iterm)

            ! Neutral
            if (istg.eq.1) then
              f0 = f0 + ltepopu(ilevel)

            ! Ion 1
            else if (istg.eq.2) then
              f1 = f1 + ltepopu(ilevel)

            ! Ion 2
            else if (istg.eq.3) then
              f2 = f2 + ltepopu(ilevel)

            end if

          end do ! FS level
        end do ! Term

        ! Deallocate ltepopu
        deallocate(ltepopu)

      ! If populations
      else

        ! Initialize level
        ilevel = 0

        ! For each level
        do iterm=1,Atom%nmulti
          do iJ=1,Atom%nJ(iterm)

            ! Advance index
            ilevel = ilevel + 1

            ! Stage
            istg = Atom%stage(iterm)

            ! Neutral
            if (istg.eq.1) then
              f0 = f0 + Atom%popu(ilevel,iz)

            ! Ion 1
            else if (istg.eq.2) then
              f1 = f1 + Atom%popu(ilevel,iz)

            ! Ion 2
            else if (istg.eq.3) then
              f2 = f2 + Atom%popu(ilevel,iz)

            end if ! Stage

          end do ! FS level
        end do ! Term

      end if ! Populations or departure

      ! If valid first ion
      if (f1.gt.0d0) then

        ! Second to first ion
        f2 = f2/f1

      ! Nothing in first ion
      else

        ! Cannot be anything in second ion either
        f2 = 0d0

      end if

      ! If valid neutral
      if (f0.gt.0d0) then

        ! First ion to neutral
        f1 = f1/f0

      ! Nothing in neutral
      else

        ! Cannot be anythin in first ion either
        f1 = 0d0

      end if

      end subroutine metal_fraction

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute ionization fraction of stages in an atomic model from
      !! CHIANTI data for a given temperature\n
      !!        Atom(Atom_class): Structure with atomic data\n
      !!               T(double): Temperature\n
      !!  chianti(chianti_class): Structure with the CHIANTI data\n
      !!             frc(double): Ionization fraction
      subroutine chianti_fraction(Atom,T,chianti,frc)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(chianti_class), intent(in):: chianti
      double precision, intent(in):: T
      double precision, intent(out):: frc

      ! Local

      integer:: iatom,istg,minstg,maxstg
      double precision:: cT,lfrc

      ! Get atom index
      iatom = atom_char2index(Atom%Element)

      ! If index outside of CHIANTI data, abort
      if (iatom.gt.chianti%nE) then

        ! Report error
        urou = 'chianti_fraction'
        umsg = 'The element '//Atom%Element//' is '// &
               'not included in the CHIANTI data. '// &
               'aborting the run'
        call verbose
        frc = 0d0
        laborted = .True.

        ! return
        return

      end if ! Atom not in CHIANTI database

      ! Get minimum and maximum stages in the atom
      minstg = minval(Atom%stage)
      maxstg = maxval(Atom%stage)

      ! If zeroing last stage, reduce maximum
      if (Atom%zero_ion) maxstg = maxstg - 1

      ! Initialize fraction
      frc = 0d0

      ! Current temperature in CHIANTI units
      cT = log10(T)

      ! For each stage
      do istg=minstg,maxstg

        ! Evaluate spline interpolation in CHIANTI tabulation
        lfrc = ispline(cT, &
                       dble(chianti%ioneq_T), &
                       dble(chianti%ioneq(iatom)%stage(istg)%p), &
                       chianti%ioneq(iatom)%b(:,istg), &
                       chianti%ioneq(iatom)%c(:,istg), &
                       chianti%ioneq(iatom)%d(:,istg), &
                       chianti%nT)

        ! If negative fraction
        if (lfrc.lt.0d0) then

          ! Do linear interpolation instead
          call linear(dble(chianti%ioneq_T), &
                      dble(chianti%ioneq(iatom)%stage(istg)%p), &
                      cT,lfrc)

        end if ! Wrong spline interpolation

        ! Add to total fraction
        frc = frc + lfrc

      end do ! Stages

      return

      end subroutine chianti_fraction

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute ion fraction from populations in Hydrogen model\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iz(integer): Height index\n
      !!        f1(double): Ion II to I\n
      !!    depar(logical): If using departure coefficients
      subroutine H_fraction(Atom,Atmo,iz,f1,depar)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: depar
      integer, intent(in):: iz
      double precision, intent(out):: f1

      ! Local

      integer:: iterm, iJ, ilevel, istg

      double precision:: f0
      double precision, dimension(:), allocatable:: ltepopu


      ! If not neutral or first stage in the model
      if (minval(Atom%stage).ge.1.or. &
          maxval(Atom%stage).le.2) then

        ! Return invalid values
        f1 = -1d0
        return

      end if ! No neutral or first stage

      ! Initialize quantities
      f0 = 0d0
      f1 = 0d0

      ! If departure coefficients
      if (depar) then

        ! Allocate ltepop
        allocate(ltepopu(Atom%nlevel))

        ! Compute LTE
        call LTEiz(Atom,Atmo,iz,ltepopu)

        ! Apply departure
        ltepopu = ltepopu*Atom%depar(:,iz)

        ! Initialize level
        ilevel = 0

        ! For each level
        do iterm=1,Atom%nmulti
          do iJ=1,Atom%nJ(iterm)

            ! Advance index
            ilevel = ilevel + 1

            ! Stage
            istg = Atom%stage(iterm)

            ! Neutral
            if (istg.eq.1) then
              f0 = f0 + ltepopu(ilevel)

            ! Ion 1
            else if (istg.eq.2) then
              f1 = f1 + ltepopu(ilevel)

            end if

          end do ! FS level
        end do ! Term

        ! Deallocate ltepopu
        deallocate(ltepopu)

      ! If populations
      else

        ! Initialize level
        ilevel = 0

        ! For each level
        do iterm=1,Atom%nmulti
          do iJ=1,Atom%nJ(iterm)

            ! Advance index
            ilevel = ilevel + 1

            ! Stage
            istg = Atom%stage(iterm)

            ! Neutral
            if (istg.eq.1) then
              f0 = f0 + Atom%popu(ilevel,iz)

            ! Ion 1
            else if (istg.eq.2) then
              f1 = f1 + Atom%popu(ilevel,iz)

            end if

          end do ! FS level
        end do ! Term

      end if ! Populations or departure

      ! If valid neutral
      if (f0.gt.0d0) then

        ! First ion to neutral
        f1 = f1/f0

      ! Nothing in neutral
      else

        ! Cannot be anythin in first ion either
        f1 = 0d0

      end if

      end subroutine H_fraction

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute LTE populations at a given height\n
      !!  Atom(Atom_class): Structure with atomic data\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iz(integer): Height index\n
      !!   popu(double(:)): LTE populations
      subroutine LTEiz(Atom,Atmo,iz,popu)

      ! I/O

      type(Atom_class), intent(in):: Atom
      type(Atmo_class), intent(in):: Atmo
      integer, intent(in):: iz
      double precision, dimension(:), intent(out):: popu

      ! Local
      integer:: iterm,iJ,ilevel,zm,iim,ii,dZ
      double precision:: C0,C1,C2,Tin,dby,S,dE,arg,ne,T,gi0
      double precision, dimension(:), allocatable:: debey

      ! Get variables from atmosphere
      T = Atmo%T(iz)
      ne = Atmo%ne(iz)

      ! Allocate debey
      allocate(debey(Atom%nlevel))

      ! Constants
      C0 = hplanck*hplanck/2d0/PI/me/kb
      ! 1d3 sqrt(cm^-3 -> m^-3) electron density below
      C1 = sqrt(8d0*PI/kb)*((qel*qel/pi4eps0)**(1.5d0))*1d3

      ! Debey correction to the ionization potential
      debey = 0d0

      ! Initialize the index
      ilevel = 1

      ! Run over all the levels using the term-J indexes
      do iterm=1,Atom%nMulti
        do iJ=1,Atom%nJ(iterm)

          ! Ignore the ground level
          if (iterm.eq.1.and.iJ.eq.1) cycle

          ! Determine the level index, the charge of the
          ! shell nucleus + rest of electrons and the
          ! change of stages between this level and the
          ! ground level of the model
          ilevel = ilevel + 1
          zm = Atom%stage(iterm) - 1
          iim = Atom%stage(iterm) - Atom%stage(1)

          ! Add the contribution to the Debey correction
          do ii=1,iim

            ! Contribution to debey value
            debey(ilevel) = debey(ilevel) + zm
            zm = zm + 1

          end do ! Stage difference
        end do ! J levels
      end do ! Term

      !
      ! Compute LTE populations
      !

      ! Inverse of temperature
      Tin = 1d0/T
      ! Multiplicative factor for Debey correction
      dby = C1*sqrt(ne*Tin)
      ! Non atomic part of the Saha function
      C2 = .5d6*ne*((C0*Tin)**(1.5d0))
      ! Cumulative factor reset
      S = 1d0

      ! Initialize index
      ilevel = 1

      ! For each level running through term and J indexes
      do iterm=1,Atom%nMulti
        do iJ=1,Atom%nJ(iterm)

          ! Skip ground level
          if (iterm.eq.1.and.iJ.eq.1) cycle

          ! Determine the index and the ionization potential
          ilevel = ilevel + 1
          dE = (Atom%FSfreq(iJ,iterm) - Atom%FSfreq(1,1))*fktoJ

          ! If hard-coded hydrogen model
          if (Atom%cust) then

            ! Use the degeneration variable
            gi0 = Atom%deg(ilevel)/Atom%deg(1)

          ! If it is a normal atom
          else

            ! Calculate degenerations
            gi0 = (2d0*Atom%rJval(iJ,iterm) + 1d0)/ &
                  (2d0*Atom%rJval(1,1) + 1d0)

          end if

          ! Determine the difference in charge with the ground state
          dZ = Atom%stage(iterm) - Atom%stage(1)

          ! Calculate the argument of the exponential, ionization
          ! potential with Debey correction
          arg = (debey(ilevel)*dby - dE)*Tin/kb

          ! Numerator of the population solution
          if (arg.lt.0d0) then
            arg = -arg
            popu(ilevel) = gi0*diexp(arg)
          else
            popu(ilevel) = gi0*ddexp(arg)
          end if

          ! Denominator
          do ii=1,dZ
            popu(ilevel) = popu(ilevel)/C2
          end do

          ! Accumulate for the normalization equation
          S = S + popu(ilevel)

        end do ! J levels
      end do ! Terms

      ! Get the ground level population from the accumulated
      ! factor
      popu(1) = Atom%abun_mod*Atom%abun*Atmo%nht(iz)/S

      ! Determine the population of each level from the
      ! population of the ground state
      do ilevel=2,Atom%nlevel
        popu(ilevel) = popu(ilevel)*popu(1)
      end do

      ! Deallocate debey
      deallocate(debey)

      end subroutine LTEiz

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set the atmospheric quantities from the pressures\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!       iz(integer): Height index\n
      !!   Pall(double(:)): Array of pressures
      subroutine set_densities(Atmo,iz,Pall)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      integer, intent(in):: iz
      double precision, dimension(:), intent(in):: Pall

      ! Local

      integer:: il

      double precision:: ikbcgsT,ikbT,PHtot,S,nH,gi0,dE
      double precision, dimension(5):: deg,EH


      ! Constants
      ikbT = 1d0/kb/Atmo%T(iz)
      ikbcgsT = 1d-7*ikbT

      ! Pressures
      Atmo%Pg(iz) = Pall(1)
      Atmo%Pe(iz) = Pall(2)

      ! Electron density
      Atmo%ne(iz) = Pall(2)*ikbcgsT

      ! Get Hydrogen densities
      PHtot = Pall(3)*ikbcgsT
      nH = Pall(4)*PHtot
      Atmo%nh(iz,6) = Pall(5)*PHtot
      Atmo%nHa(iz) = nH + Atmo%nh(iz,6)
      Atmo%nHt(iz) = PHtot

      !
      ! LTE excitation balance of neutral hydrogen
      !

      ! Cumulative factor reset
      S = 1d0

      ! Hydrogen quantities
      deg = (/ 2d0, 8d0, 18d0, 32d0, 50d0 /)
      EH = (/ 0d0, .82258211d0, .97491219d0, &
              1.02822766d0, 1.05290508d0 /)

      ! Run over bound-bound levels
      do il=2,5

        ! Difference of energy divided by kT and ratio of weights
        dE = (EH(il) - EH(1))*fktoJ*ikbT
        gi0 = deg(il)/deg(1)

        ! Numerator of the population solution
        Atmo%nH(iz,il) = gi0*diexp(dE)

        ! Accumulate for the normalization equation
        S = S + Atmo%nH(iz,il)

      end do ! Levels

      ! Get the ground level population from the accumulated
      ! factor
      Atmo%nH(iz,1) = nH/S

      ! Determine the population of each level from the
      ! population of the ground state
      do il=2,5
        Atmo%nH(iz,il) = Atmo%nH(iz,il)*Atmo%nH(iz,1)
      end do

      ! Check if everything is fine
      call control

      end subroutine set_densities

!#####################################################################
!#####################################################################
!#####################################################################

      !> Set the hydrogen density atmospheric quantities from the
      !! computed partial pressures
      !!   Atmo(Atmo_class): Structure with atmospheric data\n
      !!        iz(integer): Height index
      subroutine set_Hdensities(Atmo,iz)

      ! I/O

      type(Atmo_class), intent(inout):: Atmo
      integer, intent(in):: iz

      ! Local

      integer:: il

      double precision:: ikbT,S,nH,gi0,dE
      double precision, dimension(5):: deg,EH


      ! Constant
      ikbT = 1d0/kb/Atmo%T(iz)

      !
      ! LTE excitation balance of neutral hydrogen
      !

      ! Total neutral H
      nH = Atmo%nH(iz,1)

      ! Cumulative factor reset
      S = 1d0

      ! Hydrogen quantities
      deg = (/ 2d0, 8d0, 18d0, 32d0, 50d0 /)
      EH = (/ 0d0, .82258211d0, .97491219d0, &
              1.02822766d0, 1.05290508d0 /)

      ! Run over bound-bound levels
      do il=2,5

        ! Difference of energy divided by kT and ratio of weights
        dE = (EH(il) - EH(1))*fktoJ*ikbT
        gi0 = deg(il)/deg(1)

        ! Numerator of the population solution
        Atmo%nH(iz,il) = gi0*diexp(dE)

        ! Accumulate for the normalization equation
        S = S + Atmo%nH(iz,il)

      end do ! Levels

      ! Get the ground level population from the accumulated
      ! factor
      Atmo%nH(iz,1) = nH/S

      ! Determine the population of each level from the
      ! population of the ground state
      do il=2,5
        Atmo%nH(iz,il) = Atmo%nH(iz,il)*Atmo%nH(iz,1)
      end do

      ! Check if everything is fine
      call control

      end subroutine set_Hdensities

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute Saha factor for equation of state\n
      !! Saha-Eggert equation\n
      !!   t(double): Temperature in eV\n
      !!  Ei(double): Ionization energy\n
      !!  U0(double): lower stage partition function\n
      !!  U1(double): upper stage partition function\n
      !!  Pe(double): electron pressure
      double precision function fsaha(t,Ei,U0,U1,Pe)

      double precision, intent(in):: t,Ei,U0,U1,Pe

      double precision:: den

      ! Denominator
      den = U0*Pe*(t**2.5d0)

      ! Factor
      fsaha = U1*10d0**(9.0805126d0 - t*Ei)/den

      end function fsaha

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get molecular quantities for equation of state. Data taken
      !! from the SIR code\n
      !!   t(double): Temperature in eV\n
      !!  i(integer): Molecule index
      double precision function moldata_ind(t,i)

      ! I/O

      integer:: i
      double precision:: t


      ! Choose case depending on index
      select case (i)
        case (1) ! H2+
          moldata_ind = -11.206998d0 + t*(2.7942767d0 + &
                                          t*(7.9196803d-2 - &
                                             t*2.4790744d-2))
        case (2) ! H2
          moldata_ind = -12.533505d0 + t*(4.9251644d0 + &
                                          t*(-5.6191273d-2 + &
                                             t*3.2687661d-3))
        case (3) ! CH
          moldata_ind = -11.89d0 + t*(3.8084d0 - t*2.4402d-2)
        case (4) ! NH
          moldata_ind = -11.85d0 + t*(4.1411d0 + t*(-6.7847d-2 + &
                                                    4.9178d-3*t))
        case (5) ! OH
          moldata_ind = -12.199d0 + t*(4.884d0 + t*(-7.2794d-2 + &
                                                    5.1747d-3*t))
        case (6) ! SiH
          moldata_ind = -11.205d0 + 2.0112d0*t
        case (7) ! MgH
          moldata_ind = -10.514d0 + t*(2.2206d0 - 1.3654d-2*t)
        case (8) ! CaH
          moldata_ind = -10.581d0 + t*(2.1713d0 + t*(-7.7446d-2 + &
                                                     6.5014d-3*t))
        case (9) ! AlH
          moldata_ind = -11.711d0 + t*(3.1571d0 - 1.5205d-2*t)
        case (10) ! HCl
          moldata_ind = -12.24d0 + t*(4.8009d0 - 2.6693d-2*t)
        case (11) ! HS
          moldata_ind = -11.849d0 + t*(4.1621d0 - 4.1213d-2*t)
        case (12) ! C2
          moldata_ind = -12.6d0 + t*(6.3336d0 - 1.2019d-2*t)
        case (13) ! CN
          moldata_ind = -12.667d0 + t*(7.779d0 - 1.1674d-2*t)
        case (14) ! CO
          moldata_ind = -13.641d0 + t*(11.591d0 + &
                                       t*(-8.8848d-2 + &
                                          t*7.3465d-3))
        case (15) ! SiC
          moldata_ind = -12.07d0 + t*(4.7321d0 - 1.4647d-2*t)
        case (16) ! CS
          moldata_ind = -13.122d0 + t*(8.1528d0 - 2.5302d-2*t)
        case (17) ! N2
          moldata_ind = -13.435d0 + t*(10.541d0 + &
                                       t*(-2.8061d-1 + &
                                          t*(5.883d-2 + &
                                             t*4.6609d-3)))
        case (18) ! NO
          moldata_ind = -12.606d0 + t*(6.9634d0 + &
                                       t*(-8.2021d-2 + &
                                          t*6.8654d-3))
        case (19) ! SiN
          moldata_ind = -12.43d0 + t*(5.409d0 - 2.341d-2*t)
        case (20) ! SN
          moldata_ind = -12.172d0 + t*(5.2755d0 -2.3355d-2*t)
        case (21) ! O2
          moldata_ind = -13.087d0 + t*(5.3673d0 - 1.4064d-2*t)
        case (22) ! SiO
          moldata_ind = -13.034d0 + 8.3616d0*t
        case (23) ! MgO
          moldata_ind = -11.39d0 + t*(4.267d0 -1.7738d-2*t)
        case (24) ! AlO
          moldata_ind = -12.493d0 + t*(5.3382d0 + &
                                       t*(-5.8504d-2 + &
                                          4.8035d-3*t))
        case (25) ! TiO
          moldata_ind = -13.367d0 + t*(8.869d0 + &
                                       t*(-7.1389d-1 + &
                                          t*(1.5259d-1 - &
                                             1.1909d-2*t)))
        case (26) ! VO
          moldata_ind = -12.876d0 + t*(6.8208d0 + &
                                       t*(-7.5817d-2 + &
                                          6.3995d-3*t))
        case (27) ! ZrO
          moldata_ind = -13.368d0 + t*(9.0764d0 + &
                                       t*(-2.8354d-1 + &
                                          2.5784d-2*t))
        case (28) ! SO
          moldata_ind = -12.645d0 + t*(5.6644d0 - 2.2882d-2*t)
        case (29) ! NaCl
          moldata_ind = -11.336d0 + t*(4.4639d0 - 1.6552d-2*t)
        case (30) ! SiS
          moldata_ind = -12.515d0 + 6.7906d0*t
        case (31) ! CaCl
          moldata_ind = -10.638d0 + 4.2139d0*t
        case (32) ! AlCl
          moldata_ind = -12.06d0 + t*(5.3309d0 - 1.6459d-2*t)
        case (33) ! Cl2
          moldata_ind = -12.508d0 + t*(2.727d0 - 1.7697d-2*t)
        case (34) ! S2
          moldata_ind = -12.651d0 + t*(4.697d0 - 2.5267d-2*t)
        case (35) ! CH2
          moldata_ind = -24.883d0 + t*(8.2225d0 - 2.6757d-1*t)
        case (36) ! NH2
          moldata_ind = -24.82d0 + t*(8.4594d0 + &
                                      t*(-1.2208d-1 + &
                                         8.6774d-3*t))
        case (37) ! H2O
          moldata_ind = -25.206d0 + t*(10.311d0 + &
                                       t*(-9.0573d-2 + &
                                          5.3739d-3*t))
        case (38) ! H2S
          moldata_ind = -24.314d0 + t*(8.1063d0 - 3.4079d-2*t)
        case (39) ! HCN
          moldata_ind = -25.168d0 + 13.401d0*t
        case (40) ! HCO
          moldata_ind = -25.103d0 + t*(12.87d0 - 3.8336d-2*t)
        case (41) ! HNO
          moldata_ind = -25.078d0 + t*(9.3435d0 + &
                                       t*(-1.06d-1 + &
                                          8.2469d-3*t))
        case (42) ! HO2
          moldata_ind = -25.161d0 + t*(7.8472d0 + &
                                       t*(-1.0399d-1 + &
                                          7.8209d-3*t))
        case (43) ! C3
          moldata_ind = -27.038d0 + t*(14.376d0 + &
                                       t*(6.8899d-2 + &
                                          6.0371d-2*t))
        case (44) ! SiC2
          moldata_ind = -25.889d0 + 13.317d0*t
        case (45) ! CO2
          moldata_ind = -27.261d0 + t*(16.866d0 - 1.0144d-2*t)
        case (46) ! N2O
          moldata_ind = -26.25d0 + t*(11.83d0 + &
                                      t*(-4.2021d-2 + &
                                         3.4397d-3*t))
        case (47) ! NO2
          moldata_ind = -26.098d0 + t*(10.26d0 + &
                                       t*(-1.0101d-1 + &
                                          8.4813d-3*t))
        case (48) ! O3
          moldata_ind = -26.115d0 + t*(6.5385d0 - 1.9332d-2*t)
        case (49) ! TiO2
          moldata_ind = -27.496d0 + t*(13.549d0 + &
                                       t*(-2.205d-2 + &
                                          2.1407d-2*t))
        case (50) ! ZrO2
          moldata_ind = -27.494d0 + t*(15.99d0 + &
                                       t*(-2.3605d-1 + &
                                          2.1644d-2*t))
        case (51) ! Al2O
          moldata_ind = -25.244d0 + t*(11.065d0 - 3.7924d-2*t)
        case (52) ! AlCl2
          moldata_ind = -23.748d0 + t*(9.5556d0 - 2.4867d-2*t)
        case (53) ! CH3
          moldata_ind = -37.194d0 + t*(13.371d0 - 3.4229d-2*t)
        case (54) ! NH3
          moldata_ind = -37.544d0 + t*(12.895d0 - 4.9012d-2*t)
        case (55) ! C2H2
          moldata_ind = -37.931d0 + 17.216d0*t
        case (56) ! HCOH
          moldata_ind = -38.274d0 + t*(16.264d0 - 3.2379d-2*t)
        case (57) ! HCNO
          moldata_ind = -38.841d0 + t*(18.907d0 - 3.5705d-2*t)
        case (58) ! CH4
          moldata_ind = -50.807d0 + t*(17.956d0 - 3.6861d-2*t)
        case (59) ! NaH
          moldata_ind = -11.4575d0 + t*(3.1080922d0 + &
                                        t*(-3.3159806d-1 + &
                                           4.314945d-2*t))
        case (60) ! KH
          moldata_ind = -10.964723d0 + t*(2.270225d0 + &
                                          t*(-7.66888d-2 + &
                                             6.519213d-3*t))
        case (61) ! BeH
          moldata_ind = -10.807839d0 + t*(2.744854d0 + &
                                          t*(5.758024d-2 + &
                                             3.315373d-3*t))
        case (62) ! SrH
          moldata_ind = -10.491008d0 + t*(2.051217d0 + &
                                          t*(-7.643729d-2 + &
                                             6.425358d-3*t))
        case (63) ! SrO
          moldata_ind = -11.01929d0 + t*(3.13829d0 + &
                                         t*(1.214975d0 - &
                                            1.77077d-1*t))
        case (64) ! BaH
          moldata_ind = -10.446909d0 + t*(2.024548d0 + &
                                          t*(-7.680739d-2 + &
                                             6.471443d-3*t))
        case (65) ! BaO
          moldata_ind = -10.921254d0 + t*(3.847116d0 + &
                                          t*(1.189653d0 - &
                                             1.662815d-1*t))
        case (66) ! ScO
          moldata_ind = -13.561415d0 + t*(7.528355d0 + &
                                          t*(-5.031809d-1 + &
                                             6.787392d-2*t))
        case (67) ! YO
          moldata_ind = -14.107593d0 + t*(12.226315d0 + &
                                          t*(-1.019148d0 + &
                                             1.02046d-1*t))
        case (68) ! LaO
          moldata_ind = -14.231303d0 + t*(11.28907d0 + &
                                          t*(-1.108545d0 + &
                                             1.274056d-1*t))
        case (69) ! Si2
          moldata_ind = -12.03193d0 + t*(3.012432d0 + &
                                         t*(1.798894d-1 - &
                                            1.79236d-2*t))
        case (70) ! LiH
          moldata_ind = -11.344906d0 + t*(2.836732d0 + &
                                          t*(-1.134115d-1 + &
                                             1.99661d-2*t))
        case (71) ! VO2
          moldata_ind = -25.913244d0 + t*(11.856324d0 + &
                                          t*(1.05407d0 - &
                                             1.541093d-1*t))
        case (72) ! SiO2
          moldata_ind = -26.934577d0 + t*(13.421189d0 + &
                                          t*(2.671335d-1 - &
                                             3.475224d-2*t))
        case (73) ! SiH2
          moldata_ind = -23.225499d0 + t*(4.820973d0 + &
                                          t*(6.722119d-1 - &
                                             5.903448d-2*t))
        case (74) ! Si3
          moldata_ind = -25.079417d0 + t*(7.196588d0 + &
                                          t*(1.196713d-1 + &
                                             1.0484d-2*t))
        case (75) ! C2H
          moldata_ind = -26.331315d0 + t*(12.500392d0 + &
                                          t*(6.531014d-1 - &
                                             1.162598d-1*t))
        case (76) ! BH
          moldata_ind = -11.673776d0 + t*(3.245147d0 + &
                                          t*(1.334288d-1 - &
                                             1.524113d-2*t))
        case (77) ! BO
          moldata_ind = -12.972588d0 + t*(7.80983d0 - &
                                          t*(6.263376d-2 - &
                                             4.763338d-3*t))
        case (78) ! HF
          moldata_ind = -12.654d0 + t*(6.2558d0 - 3.0868d-2*t)
        case (79) ! LiO
          moldata_ind = -11.639d0 + t*(3.9742d0 - &
                                       t*(1.7229d-1 - 1.687d-2*t))
        case (80) ! LiF
          moldata_ind = -11.92d0 + t*(6.1426d0 - 1.8981d-2*t)
        case (81) ! LiCl
          moldata_ind = -11.576d0 + t*(5.1447d0 + &
                                       t*(1.4625d-2 - &
                                          8.9675d-3*t))
        case (82) ! FeO
          moldata_ind = -12.579d0 + t*(4.8824d0 + &
                                       t*(-1.3848d-1 + &
                                          1.25d-2*t))
        case (83) ! NaF
          moldata_ind = -11.666d0 + t*(5.1748d0 - 1.9602d-2*t)
        case (84) ! MgF
          moldata_ind = -11.292d0 + t*(4.851d0 - 1.8104d-2*t)
        case (85) ! AlF
          moldata_ind = -12.453d0 + t*(7.1023d0 - 1.6086d-2*t)
        case (86) ! KF
          moldata_ind = -11.8d0 + t*(5.7274d0 - &
                                     t*(1.9201d-1 - &
                                        2.1306d-2*t))
        case (87) ! MgCl
          moldata_ind = -10.798d0 + t*(2.6664d0 + &
                                       t*(1.2037d-1 - &
                                          1.778d-2*t))
        case (88) ! KCl
          moldata_ind = -11.453d0 + t*(4.9299d0 + &
                                       t*(-1.4643d-1 + &
                                          1.4143d-2*t))
        case (89) ! FeCl
          moldata_ind = -11.484d0 + t*(3.2399d0 - 2.2356d-2*t)
        case (90) ! LiOH
          moldata_ind = -24.304d0 + t*(9.5257d0 - 4.2841d-2*t)
        case default
          urou = 'moldata_ind'
          moldata_ind = 0d0
          write(umsg,'(A,i2,A)') 'Index ',i, &
                                 ' not found in molecular data list'
          call aborted
          return
      end select

      end function moldata_ind

!#####################################################################
!#####################################################################
!#####################################################################

      end module chemicaux_mod
