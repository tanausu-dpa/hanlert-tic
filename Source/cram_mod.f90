      !> Count memory
      module cram_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     24/10/2024
!  Last version:
!     15/05/2025 V4.0.3
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     15/05/2025:    V4.0.3 - Generalized declarations of Atom, Atomb,
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
!#####################################################################
!#####################################################################
!
!  Data:
!
!  cram
!    Count RAM allocated in all structures entering hanle()
!
!  cram_report
!    Report RAM allocated in all structures entering hanle() and
!    contribute to MRAMc
!
!  cram_solf
!    Count memory allocated in the SolF variable
!
!  cram_atom
!    Count memory allocated in the Atom variable
!
!  cram_ltelines
!    Count memory allocated in the LTElines variable
!
!  cram_mol
!    Count memory allocated in the Mol variable
!
!  cram_atmo
!    Count memory allocated in the Atmo variable
!
!  cram_mpi
!    Count memory allocated in the MPID variable
!
!  cram_input
!    Count memory allocated in the Input variable
!
!  cram_geom
!    Count memory allocated in the Geom variable
!
!  cram_bfield
!    Count memory allocated in the Bfield variable
!
!  cram_frec
!    Count memory allocated in the Frec variable
!
!  cram_flgsg
!    Count memory allocated in the Flgsg variable
!
!  cram_estimate_norm
!    Estimate the amount of memory that would be used for
!  normalization data
!
!  cram_add
!    Compute the amount of RAM that is already reserved
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use commons_mod
      use parameters_mod , only : TINYVEL, TINYB
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count RAM allocated in all structures entering hanle()\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atomb(Atom_class(:)): Structures with atomic data for
      !!                             background atoms\n
      !! LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!          Mol(Mol_class(:)): Structures with molecular data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with geometric data for
      !!                             the intensity problem\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!           JKQin(double(:)): Data with ad-hoc JKQ tensors\n
      !!     SolF(Solution_F_class): Structure with the solution of
      !!                             the self-consistent problem and
      !!                             the corresponding emergent
      !!                             profiles, contribution function,
      !!                             and height for optical depth
      !!                             equal to one
      subroutine cram(Atom,Atomb,LTElines,Mol,Atmo,MPID,Input, &
                      GeomI,Geom,Bfield,Frec,Flgsg,fudge,kurucz, &
                      JKQin,SolF)

      ! I/O

      type(Atom_class), dimension(:), allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), allocatable, intent(in):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Mol_class), dimension(:), allocatable, intent(in):: Mol
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(in):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: GeomI,Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      type(Solution_F_class), intent(in):: SolF
      double precision, dimension(:), allocatable, intent(in):: JKQin

      ! Local

      double precision:: num


      ! Reset memory nums
      BRAMc = 0d0
      PRAMc = 0d0
      VRAMc = 0d0
      WRAMc = 0d0
      RRAMc = 0d0
      SRAMc = 0d0
      MRAMc = 0d0
      TRAMc = 0d0
      ORAMc = 0d0
      FRAMc = 0d0
      DRAMc = 0d0

      ! Add size of SolF to radiation
      call cram_solf(SolF,num)
      SRAMc = SRAMc + num

      !
      ! Count memory misc.
      !

      ! Atoms
      if (nA.gt.0) then
        call cram_atom(Atom,num)
        MRAMc = MRAMc + num
      end if
      if (allocated(Atomb)) then
        call cram_atom(Atomb,num)
        MRAMc = MRAMc + num
      end if

      ! LTE lines
      if (allocated(LTElines)) then
        call cram_ltelines(LTElines,num)
        MRAMc = MRAMc + num
      end if

      ! Molecules
      if (allocated(Mol)) then
        call cram_mol(Mol,num)
        MRAMc = MRAMc + num
      end if

      ! Atmosphere
      call cram_atmo(Atmo,num)
      MRAMc = MRAMc + num

      ! MPID
      call cram_mpi(MPID,num)
      MRAMc = MRAMc + num

      ! Input
      call cram_input(Input,num)
      MRAMc = MRAMc + num

      ! Geometry
      call cram_geom(GeomI,num)
      MRAMc = MRAMc + num
      call cram_geom(Geom,num)
      MRAMc = MRAMc + num

      ! Bfield
      call cram_bfield(Bfield,num)
      MRAMc = MRAMc + num

      ! Frec
      call cram_frec(Frec,num)
      MRAMc = MRAMc + num

      ! Flgsg
      call cram_flgsg(Flgsg,num)
      MRAMc = MRAMc + num

      ! Memoization
     !ERAMc = ERAMc + 8d-6*dble(nJs)

      ! Fudge
      call cram_fudge(fudge,num)
      MRAMc = MRAMc + num

      ! Kurucz
      call cram_kurucz(kurucz,num)
      MRAMc = MRAMc + num

      ! JKQin
      RRAMc = RRAMc + 1d-6*sizeof(JKQin)

      end subroutine cram

!#####################################################################
!#####################################################################
!#####################################################################

      !> Report RAM allocated in all structures entering hanle()
      !! and contribute to MRAMc\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !!       Atomb(Atom_class(:)): Structures with atomic data for
      !!                             background atoms\n
      !! LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!          Mol(Mol_class(:)): Structures with molecular data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!            MPID(MPI_class): Structure with MPI data\n
      !!         Input(Input_class): Structure with configuration
      !!                             data\n
      !!      GeomI(Geometry_class): Structure with geometric data for
      !!                             the intensity problem\n
      !!       Geom(Geometry_class): Structure with geometric data\n
      !!       Bfield(Bfield_class): Structure with magnetic field
      !!                             data\n
      !!      Frec(Frequency_class): Structure with frequency data\n
      !!         Flgsg(Fctsg_class): Structure with factorials, signs,
      !!                             and J-symbols\n
      !!         fudge(fudge_class): Structure with fudge data\n
      !!       kurucz(kurucz_class): Structure with Kurucz line data\n
      !!              init(logical): If initializing file
      subroutine cram_report(Atom,Atomb,LTElines,Mol,Atmo,MPID, &
                             Input,GeomI,Geom,Bfield,Frec,Flgsg, &
                             fudge,kurucz,init)

      ! I/O

      type(Atom_class), dimension(:), allocatable, intent(in):: Atom
      type(Atom_class), dimension(:), allocatable, intent(in):: Atomb
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Mol_class), dimension(:), allocatable, intent(in):: Mol
      type(Atmo_class), intent(in):: Atmo
      type(Bfield_class), intent(in):: Bfield
      type(Fctsg_class), intent(in):: Flgsg
      type(fudge_class), intent(in):: fudge
      type(kurucz_class), intent(in):: kurucz
      type(Frequency_class), intent(in):: Frec
      type(Geometry_class), intent(in):: GeomI,Geom
      type(Input_class), intent(in):: Input
      type(MPI_class), intent(inout):: MPID
      logical, intent(in):: init

      ! Local

      character(len=9):: FIL

      double precision:: num,tot


      ! Filename
      write(FIL,'("LOG_",i0.5)') gpid

      ! File
      if (init) then
        open(800,file=FIL,action='write')
      else
        open(800,file=FIL,position='append')
      end if

      ! Initialize total
      tot = 0d0

      !
      ! Count memory misc.
      !

      ! Atoms
      write(800,*) ''
      if (nA.gt.0) then
        call cram_atom(Atom,num)
        write(800,'("Atom      ",es15.8)') num
        tot = tot + num
      end if
      if (allocated(Atomb)) then
        write(800,'("Atomb     ",es15.8)') num
        call cram_atom(Atomb,num)
        tot = tot + num
      end if

      ! LTE lines
      if (allocated(LTElines)) then
        call cram_ltelines(LTElines,num)
        write(800,'("LTElines  ",es15.8)') num
        tot = tot + num
      end if

      ! Molecules
      if (allocated(Mol)) then
        call cram_mol(Mol,num)
        write(800,'("LTElines  ",es15.8)') num
        tot = tot + num
      end if

      ! Atmosphere
      call cram_atmo(Atmo,num)
      write(800,'("Atmo      ",es15.8)') num
      tot = tot + num

      ! MPID
      call cram_mpi(MPID,num)
      write(800,'("MPID      ",es15.8)') num
      tot = tot + num

      ! Input
      call cram_input(Input,num)
      write(800,'("Input     ",es15.8)') num
      tot = tot + num

      ! Geometry
      call cram_geom(GeomI,num)
      write(800,'("GeomI     ",es15.8)') num
      tot = tot + num
      call cram_geom(Geom,num)
      write(800,'("Geom      ",es15.8)') num
      tot = tot + num

      ! Bfield
      call cram_bfield(Bfield,num)
      write(800,'("Bfield    ",es15.8)') num
      tot = tot + num

      ! Frec
      call cram_frec(Frec,num)
      write(800,'("Frec      ",es15.8)') num
      tot = tot + num

      ! Flgsg
      call cram_flgsg(Flgsg,num)
      write(800,'("Flgsg     ",es15.8)') num
      tot = tot + num

      ! Fudge
      call cram_fudge(fudge,num)
      write(800,'("fudge     ",es15.8)') num
      tot = tot + num

      ! Kurucz
      call cram_kurucz(kurucz,num)
      write(800,'("kurucz    ",es15.8)') num
      tot = tot + num

      write(800,'("Total     ",es15.8)') tot

      end subroutine cram_report

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the SolF variable\n
      !!  SolF(Solution_F_class): Structure with the solution of the
      !!                          self-consistent problem and the
      !!                          corresponding emergent profiles,
      !!                          contribution function, and height
      !!                          for optical depth equal to one\n
      !!             num(double): Memory count
      subroutine cram_solf(SolF,num)

      ! I/O

      type(Solution_F_class), intent(in):: SolF
      double precision, intent(out):: num

      ! Initialize
      num = 0d0

      ! Check arrays
      if (allocated(SolF%e_tau1)) num = num + &
                                        1d-6*sizeof(SolF%e_tau1)
      if (allocated(SolF%e_Ctr)) num = num + &
                                       1d-6*sizeof(SolF%e_Ctr)
      if (allocated(SolF%i_J00)) num = num + &
                                       1d-6*sizeof(SolF%i_J00)
      if (allocated(SolF%i_J00C)) num = num + &
                                        1d-6*sizeof(SolF%i_J00C)
      if (allocated(SolF%i_J00P)) num = num + &
                                        1d-6*sizeof(SolF%i_J00P)
      if (allocated(SolF%e_Stk)) num = num + &
                                       1d-6*sizeof(SolF%e_Stk)
      if (allocated(SolF%i_StkI)) num = num + &
                                        1d-6*sizeof(SolF%i_StkI)
      if (allocated(SolF%i_Stk)) num = num + &
                                       1d-6*sizeof(SolF%i_Stk)
      if (allocated(SolF%i_JKQ)) num = num + &
                                       1d-6*sizeof(SolF%i_JKQ)
      if (allocated(SolF%i_JKQS)) num = num + &
                                        1d-6*sizeof(SolF%i_JKQS)
      if (allocated(SolF%i_JKQC)) num = num + &
                                        1d-6*sizeof(SolF%i_JKQC)
      if (allocated(SolF%e_tau1_t)) num = num + &
                                          1d-6*sizeof(SolF%e_tau1_t)
      if (allocated(SolF%e_Ctr_t)) num = num + &
                                         1d-6*sizeof(SolF%e_Ctr_t)
      if (allocated(SolF%i_J00_t)) num = num + &
                                         1d-6*sizeof(SolF%i_J00_t)
      if (allocated(SolF%i_J00C_t)) num = num + &
                                          1d-6*sizeof(SolF%i_J00C_t)
      if (allocated(SolF%i_J00P_t)) num = num + &
                                          1d-6*sizeof(SolF%i_J00P_t)
      if (allocated(SolF%e_Stk_t)) num = num + &
                                         1d-6*sizeof(SolF%e_Stk_t)
      if (allocated(SolF%i_StkI_t)) num = num + &
                                          1d-6*sizeof(SolF%i_StkI_t)
      if (allocated(SolF%i_Stk_t)) num = num + &
                                         1d-6*sizeof(SolF%i_Stk_t)
      if (allocated(SolF%i_JKQ_t)) num = num + &
                                         1d-6*sizeof(SolF%i_JKQ_t)
      if (allocated(SolF%i_JKQS_t)) num = num + &
                                          1d-6*sizeof(SolF%i_JKQS_t)
      if (allocated(SolF%i_JKQC_t)) num = num + &
                                          1d-6*sizeof(SolF%i_JKQC_t)
      if (allocated(SolF%e_tau1_b)) num = num + &
                                          1d-6*sizeof(SolF%e_tau1_b)
      if (allocated(SolF%e_Ctr_b)) num = num + &
                                         1d-6*sizeof(SolF%e_Ctr_b)
      if (allocated(SolF%i_J00_b)) num = num + &
                                         1d-6*sizeof(SolF%i_J00_b)
      if (allocated(SolF%i_J00C_b)) num = num + &
                                          1d-6*sizeof(SolF%i_J00C_b)
      if (allocated(SolF%i_J00P_b)) num = num + &
                                          1d-6*sizeof(SolF%i_J00P_b)
      if (allocated(SolF%e_Stk_b)) num = num + &
                                         1d-6*sizeof(SolF%e_Stk_b)
      if (allocated(SolF%i_StkI_b)) num = num + &
                                          1d-6*sizeof(SolF%i_StkI_b)
      if (allocated(SolF%i_Stk_b)) num = num + &
                                         1d-6*sizeof(SolF%i_Stk_b)
      if (allocated(SolF%i_JKQ_b)) num = num + &
                                         1d-6*sizeof(SolF%i_JKQ_b)
      if (allocated(SolF%i_JKQS_b)) num = num + &
                                          1d-6*sizeof(SolF%i_JKQS_b)
      if (allocated(SolF%i_JKQC_b)) num = num + &
                                          1d-6*sizeof(SolF%i_JKQC_b)
      if (allocated(SolF%i_rhoes)) num = num + &
                                         1d-6*sizeof(SolF%i_rhoes)
      if (allocated(SolF%i_rhoes_t)) num = num + &
                                           1d-6*sizeof(SolF%i_rhoes_t)
      if (allocated(SolF%i_rhoes_b)) num = num + &
                                          1d-6*sizeof(SolF%i_rhoes_b)

      end subroutine cram_solf

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Atom variable\n
      !!  Atom(Atom_class(:)): Structures with atomic data\n
      !!          num(double): Memory count
      subroutine cram_atom(Atom,num)

      ! I/O

      type(Atom_class), dimension(:), allocatable, intent(in):: Atom
      double precision, intent(out):: num

      ! Local

      integer:: ia,ios,i1,i2

      type(tmp_col_box_class), pointer:: next1
      type(Tbox_class), pointer:: next2


      ! Initialize
      num = 0d0

      ! Nullify
      nullify(next1,next2)

      ! Atoms
      do ia=lbound(Atom,1),ubound(Atom,1)

        ! Non arrays
        num = num + 1d-6*sizeof(Atom(ia))

        !
        ! Arrays
        !

        ! fflag
        if (allocated(Atom(ia)%fflag)) then
          do ios=lbound(Atom(ia)%fflag,1),ubound(Atom(ia)%fflag,1)
            num = num + 1d-6*sizeof(Atom(ia)%fflag(ios))
            if (allocated(Atom(ia)%fflag(ios)%Mabsent)) &
              num = num + 1d-6*sizeof(Atom(ia)%fflag(ios)%Mabsent)
            if (allocated(Atom(ia)%fflag(ios)%Vabsent)) &
              num = num + 1d-6*sizeof(Atom(ia)%fflag(ios)%Vabsent)
          end do
        end if

        ! fst
        if (allocated(Atom(ia)%fst)) then
          do ios=lbound(Atom(ia)%fst,1),ubound(Atom(ia)%fst,1)
            num = num + 1d-6*sizeof(Atom(ia)%fst(ios))
            if (allocated(Atom(ia)%fst(ios)%ilevell)) &
                num = num + 1d-6*sizeof(Atom(ia)%fst(ios)%ilevell)
            if (allocated(Atom(ia)%fst(ios)%ilevelu)) &
                num = num + 1d-6*sizeof(Atom(ia)%fst(ios)%ilevelu)
            if (allocated(Atom(ia)%fst(ios)%Aul)) &
                num = num + 1d-6*sizeof(Atom(ia)%fst(ios)%Aul)
            if (allocated(Atom(ia)%fst(ios)%Blu)) &
                num = num + 1d-6*sizeof(Atom(ia)%fst(ios)%Blu)
          end do
        end if

        ! phot
        if (allocated(Atom(ia)%phot)) then
          do ios=lbound(Atom(ia)%phot,1),ubound(Atom(ia)%phot,1)
            num = num + 1d-6*sizeof(Atom(ia)%phot(ios))
            if (allocated(Atom(ia)%phot(ios)%Mif0)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%Mif0)
            if (allocated(Atom(ia)%phot(ios)%Mif1)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%Mif1)
            if (allocated(Atom(ia)%phot(ios)%MW0)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%MW0)
            if (allocated(Atom(ia)%phot(ios)%MW1)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%MW1)
            if (allocated(Atom(ia)%phot(ios)%alpha)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%alpha)
            if (allocated(Atom(ia)%phot(ios)%TEI)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%TEI)
            if (allocated(Atom(ia)%phot(ios)%infreq)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%infreq)
            if (allocated(Atom(ia)%phot(ios)%inalpha)) &
                num = num + 1d-6*sizeof(Atom(ia)%phot(ios)%inalpha)
          end do
        end if

        ! trano
        if (allocated(Atom(ia)%trano)) then
          do ios=lbound(Atom(ia)%trano,1),ubound(Atom(ia)%trano,1)
            num = num + 1d-6*sizeof(Atom(ia)%trano(ios))
            if (allocated(Atom(ia)%trano(ios)%indT)) &
                num = num + 1d-6*sizeof(Atom(ia)%trano(ios)%indT)
            if (allocated(Atom(ia)%trano(ios)%indB)) &
                num = num + 1d-6*sizeof(Atom(ia)%trano(ios)%indB)
            if (allocated(Atom(ia)%trano(ios)%indNB)) &
                num = num + 1d-6*sizeof(Atom(ia)%trano(ios)%indNB)
            if (allocated(Atom(ia)%trano(ios)%trani)) then
              do i1=lbound(Atom(ia)%trano(ios)%trani,1), &
                    ubound(Atom(ia)%trano(ios)%trani,1)
                num = num + 1d-6*sizeof(Atom(ia)%trano(ios)%trani(i1))
                if (allocated(Atom(ia)%trano(ios)%trani(i1)%indB)) &
                  num = num + &
                        1d-6*sizeof(Atom(ia)%trano(ios)% &
                                    trani(i1)%indB)
                if (allocated(Atom(ia)%trano(ios)%trani(i1)%indNB)) &
                  num = num + &
                        1d-6*sizeof(Atom(ia)%trano(ios)% &
                                    trani(i1)%indNB)
              end do
            end if
          end do
        end if

        ! tranoI
        if (allocated(Atom(ia)%tranoI)) then
          do ios=lbound(Atom(ia)%tranoI,1),ubound(Atom(ia)%tranoI,1)
            num = num + 1d-6*sizeof(Atom(ia)%tranoI(ios))
            if (allocated(Atom(ia)%trano(ios)%indT)) &
                num = num + 1d-6*sizeof(Atom(ia)%trano(ios)%indT)
          end do
        end if

        ! Ccoeff_special
        if (associated(Atom(ia)%Ccoeff_special)) then
          num = num + 1d-6*sizeof(Atom(ia)%Ccoeff_special)
          if (allocated(Atom(ia)%Ccoeff_special%C)) &
            num = num + 1d-6*sizeof(Atom(ia)%Ccoeff_special%C)
          next1 => Atom(ia)%Ccoeff_special
          do while (associated(next1%next))
            next1 => next1%next
            num = num + 1d-6*sizeof(next1)
            if (allocated(next1%C)) num = num + 1d-6*sizeof(next1%C)
          end do
          nullify(next1)
        end if

        ! inelas
        if (allocated(Atom(ia)%inelas)) then
          do ios=lbound(Atom(ia)%inelas,1),ubound(Atom(ia)%inelas,1)
            num = num + 1d-6*sizeof(Atom(ia)%inelas(ios))
            if (allocated(Atom(ia)%inelas(ios)%Cul)) &
                num = num + 1d-6*sizeof(Atom(ia)%inelas(ios)%Cul)
          end do
        end if

        ! elas
        if (allocated(Atom(ia)%elas)) then
          do ios=lbound(Atom(ia)%elas,1),ubound(Atom(ia)%elas,1)
            num = num + 1d-6*sizeof(Atom(ia)%elas(ios))
            if (allocated(Atom(ia)%elas(ios)%datum)) &
                num = num + 1d-6*sizeof(Atom(ia)%elas(ios)%datum)
          end do
        end if

        ! Tbox
        if (associated(Atom(ia)%Tbox)) then
          num = num + 1d-6*sizeof(Atom(ia)%Tbox)
          if (allocated(Atom(ia)%Tbox%temp)) &
            num = num + 1d-6*sizeof(Atom(ia)%Tbox%temp)
          next2 => Atom(ia)%Tbox
          do while (associated(next2%next))
            next2 => next2%next
            num = num + 1d-6*sizeof(next2)
            if (allocated(next2%temp)) &
              num = num + 1d-6*sizeof(next2%temp)
          end do
          nullify(next2)
        end if

        ! rdip
        if (allocated(Atom(ia)%rdip)) then
          do ios=lbound(Atom(ia)%rdip,1),ubound(Atom(ia)%rdip,1)
            num = num + 1d-6*sizeof(Atom(ia)%rdip(ios))
            if (allocated(Atom(ia)%rdip(ios)%rdip)) &
                num = num + 1d-6*sizeof(Atom(ia)%rdip(ios)%rdip)
          end do
        end if

        ! rdipev
        if (allocated(Atom(ia)%rdipev)) then
          do ios=lbound(Atom(ia)%rdipev,1),ubound(Atom(ia)%rdipev,1)
            num = num + 1d-6*sizeof(Atom(ia)%rdipev(ios))
            if (allocated(Atom(ia)%rdipev(ios)%rdipev)) &
                num = num + 1d-6*sizeof(Atom(ia)%rdipev(ios)%rdipev)
          end do
        end if

        ! irho
        if (allocated(Atom(ia)%irho)) then
          do ios=lbound(Atom(ia)%irho,1),ubound(Atom(ia)%irho,1)
            num = num + 1d-6*sizeof(Atom(ia)%irho(ios))
            if (allocated(Atom(ia)%irho(ios)%Jrho)) then
              num = num + 1d-6*sizeof(Atom(ia)%irho(ios)%Jrho)
              do i1=lbound(Atom(ia)%irho(ios)%Jrho,2), &
                    ubound(Atom(ia)%irho(ios)%Jrho,2)
                do i2=lbound(Atom(ia)%irho(ios)%Jrho,1), &
                      ubound(Atom(ia)%irho(ios)%Jrho,1)
                  if (allocated(Atom(ia)%irho(ios)% &
                                 Jrho(i1,i2)%kq)) &
                    num = num + 1d-6*sizeof(Atom(ia)% &
                                            irho(ios)%Jrho(i2,i1)%kq)
                end do
              end do
            end if

            if (allocated(Atom(ia)%irho(ios)%irho_ij)) &
              num = num + 1d-6*sizeof(Atom(ia)%irho(ios)%irho_ij)
            if (allocated(Atom(ia)%irho(ios)%jM)) &
              num = num + 1d-6*sizeof(Atom(ia)%irho(ios)%jM)
          end do
        end if

        ! bbspecin
        if (allocated(Atom(ia)%bbspecin)) &
          num = num + 1d-6*sizeof(Atom(ia)%bbspecin)

        ! bfspecin
        if (allocated(Atom(ia)%bfspecin)) &
          num = num + 1d-6*sizeof(Atom(ia)%bfspecin)

        ! lemiss2
        if (allocated(Atom(ia)%lemiss2)) &
          num = num + 1d-6*sizeof(Atom(ia)%lemiss2)

        ! splitf
        if (allocated(Atom(ia)%splitf)) &
          num = num + 1d-6*sizeof(Atom(ia)%splitf)

        ! rhonull
        if (allocated(Atom(ia)%rhonull)) &
          num = num + 1d-6*sizeof(Atom(ia)%rhonull)

        ! NCHLT
        if (allocated(Atom(ia)%NCHLT)) &
          num = num + 1d-6*sizeof(Atom(ia)%NCHLT)

        ! nfreqt
        if (allocated(Atom(ia)%nfreqt)) &
          num = num + 1d-6*sizeof(Atom(ia)%nfreqt)

        ! nfreqtc
        if (allocated(Atom(ia)%nfreqtc)) &
          num = num + 1d-6*sizeof(Atom(ia)%nfreqtc)

        ! nJ
        if (allocated(Atom(ia)%nJ)) &
          num = num + 1d-6*sizeof(Atom(ia)%nJ)

        ! stage
        if (allocated(Atom(ia)%stage)) &
          num = num + 1d-6*sizeof(Atom(ia)%stage)

        ! broad_type
        if (allocated(Atom(ia)%broad_type)) &
          num = num + 1d-6*sizeof(Atom(ia)%broad_type)

        ! broad_type
        if (allocated(Atom(ia)%term)) &
          num = num + 1d-6*sizeof(Atom(ia)%term)

        ! sublevel
        if (allocated(Atom(ia)%sublevel)) &
          num = num + 1d-6*sizeof(Atom(ia)%sublevel)

        ! nfreqph
        if (allocated(Atom(ia)%nfreqph)) &
          num = num + 1d-6*sizeof(Atom(ia)%nfreqph)

        ! col_type
        if (allocated(Atom(ia)%col_type)) &
          num = num + 1d-6*sizeof(Atom(ia)%col_type)

        ! if0
        if (allocated(Atom(ia)%if0)) &
          num = num + 1d-6*sizeof(Atom(ia)%if0)

        ! if1
        if (allocated(Atom(ia)%if1)) &
          num = num + 1d-6*sizeof(Atom(ia)%if1)

        ! ifst
        if (allocated(Atom(ia)%ifst)) &
          num = num + 1d-6*sizeof(Atom(ia)%ifst)

        ! itrano
        if (allocated(Atom(ia)%itrano)) &
          num = num + 1d-6*sizeof(Atom(ia)%itrano)

        ! rif0
        if (allocated(Atom(ia)%rif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%rif0)

        ! rif1
        if (allocated(Atom(ia)%rif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%rif1)

        ! rif20
        if (allocated(Atom(ia)%rif20)) &
          num = num + 1d-6*sizeof(Atom(ia)%rif20)

        ! rif21
        if (allocated(Atom(ia)%rif21)) &
          num = num + 1d-6*sizeof(Atom(ia)%rif21)

        ! sbif0
        if (allocated(Atom(ia)%sbif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%sbif0)

        ! sbif1
        if (allocated(Atom(ia)%sbif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%sbif1)

        ! sfif0
        if (allocated(Atom(ia)%sfif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%sfif0)

        ! sfif1
        if (allocated(Atom(ia)%sfif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%sfif1)

        ! ilf0
        if (allocated(Atom(ia)%ilf0)) &
          num = num + 1d-6*sizeof(Atom(ia)%ilf0)

        ! ilf1
        if (allocated(Atom(ia)%ilf1)) &
          num = num + 1d-6*sizeof(Atom(ia)%ilf1)

        ! ipf0
        if (allocated(Atom(ia)%ipf0)) &
          num = num + 1d-6*sizeof(Atom(ia)%ipf0)

        ! ipf1
        if (allocated(Atom(ia)%ipf1)) &
          num = num + 1d-6*sizeof(Atom(ia)%ipf1)

        ! Kcut
        if (allocated(Atom(ia)%Kcut)) &
          num = num + 1d-6*sizeof(Atom(ia)%Kcut)

        ! Krad
        if (allocated(Atom(ia)%Krad)) &
          num = num + 1d-6*sizeof(Atom(ia)%Krad)

        ! tif0
        if (allocated(Atom(ia)%tif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%tif0)

        ! tif1
        if (allocated(Atom(ia)%tif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%tif1)

        ! irad
        if (allocated(Atom(ia)%irad)) &
          num = num + 1d-6*sizeof(Atom(ia)%irad)

        ! icol
        if (allocated(Atom(ia)%icol)) &
          num = num + 1d-6*sizeof(Atom(ia)%icol)

        ! iphot
        if (allocated(Atom(ia)%iphot)) &
          num = num + 1d-6*sizeof(Atom(ia)%iphot)

        ! Mif0
        if (allocated(Atom(ia)%Mif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%Mif0)

        ! Mif1
        if (allocated(Atom(ia)%Mif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%Mif1)

        ! CMif0
        if (allocated(Atom(ia)%CMif0)) &
          num = num + 1d-6*sizeof(Atom(ia)%CMif0)

        ! CMif1
        if (allocated(Atom(ia)%CMif1)) &
          num = num + 1d-6*sizeof(Atom(ia)%CMif1)

        ! ifst_ij
        if (allocated(Atom(ia)%ifst_ij)) &
          num = num + 1d-6*sizeof(Atom(ia)%ifst_ij)

        ! nblk
        if (allocated(Atom(ia)%nblk)) &
          num = num + 1d-6*sizeof(Atom(ia)%nblk)

        ! fcflag
        if (allocated(Atom(ia)%fcflag)) &
          num = num + 1d-6*sizeof(Atom(ia)%fcflag)

        ! iJval
        if (allocated(Atom(ia)%iJval)) &
          num = num + 1d-6*sizeof(Atom(ia)%iJval)

        ! rLval
        if (allocated(Atom(ia)%rLval)) &
          num = num + 1d-6*sizeof(Atom(ia)%rLval)

        ! Sval
        if (allocated(Atom(ia)%Sval)) &
          num = num + 1d-6*sizeof(Atom(ia)%Sval)

        ! TRfreq
        if (allocated(Atom(ia)%TRfreq)) &
          num = num + 1d-6*sizeof(Atom(ia)%TRfreq)

        ! Dwvl
        if (allocated(Atom(ia)%Dwvl)) &
          num = num + 1d-6*sizeof(Atom(ia)%Dwvl)

        ! Dwvlc
        if (allocated(Atom(ia)%Dwvlc)) &
          num = num + 1d-6*sizeof(Atom(ia)%Dwvlc)

        ! Dfreq
        if (allocated(Atom(ia)%Dfreq)) &
          num = num + 1d-6*sizeof(Atom(ia)%Dfreq)

        ! n
        if (allocated(Atom(ia)%n)) &
          num = num + 1d-6*sizeof(Atom(ia)%n)

        ! broad_stark
        if (allocated(Atom(ia)%broad_stark)) &
          num = num + 1d-6*sizeof(Atom(ia)%broad_stark)

        ! deg
        if (allocated(Atom(ia)%deg)) &
          num = num + 1d-6*sizeof(Atom(ia)%deg)

        ! W0
        if (allocated(Atom(ia)%W0)) &
          num = num + 1d-6*sizeof(Atom(ia)%W0)

        ! W1
        if (allocated(Atom(ia)%W1)) &
          num = num + 1d-6*sizeof(Atom(ia)%W1)

        ! gL
        if (allocated(Atom(ia)%gL)) &
          num = num + 1d-6*sizeof(Atom(ia)%gL)

        ! damp
        if (allocated(Atom(ia)%damp)) &
          num = num + 1d-6*sizeof(Atom(ia)%damp)

        ! FSfreq
        if (allocated(Atom(ia)%FSfreq)) &
          num = num + 1d-6*sizeof(Atom(ia)%FSfreq)

        ! rJval
        if (allocated(Atom(ia)%rJval)) &
          num = num + 1d-6*sizeof(Atom(ia)%rJval)

        ! Ecoeff
        if (allocated(Atom(ia)%Ecoeff)) &
          num = num + 1d-6*sizeof(Atom(ia)%Ecoeff)

        ! broad_args
        if (allocated(Atom(ia)%broad_args)) &
          num = num + 1d-6*sizeof(Atom(ia)%broad_args)

        ! ldamp
        if (allocated(Atom(ia)%ldamp)) &
          num = num + 1d-6*sizeof(Atom(ia)%ldamp)

        ! popu
        if (allocated(Atom(ia)%popu)) &
          num = num + 1d-6*sizeof(Atom(ia)%popu)

        ! populte
        if (allocated(Atom(ia)%populte)) &
          num = num + 1d-6*sizeof(Atom(ia)%populte)

        ! depar
        if (allocated(Atom(ia)%depar)) &
          num = num + 1d-6*sizeof(Atom(ia)%depar)

        ! MW0
        if (allocated(Atom(ia)%MW0)) &
          num = num + 1d-6*sizeof(Atom(ia)%MW0)

        ! MW1
        if (allocated(Atom(ia)%MW1)) &
          num = num + 1d-6*sizeof(Atom(ia)%MW1)

        ! qel
        if (allocated(Atom(ia)%qel)) &
          num = num + 1d-6*sizeof(Atom(ia)%qel)

        ! Ccoeff
        if (allocated(Atom(ia)%Ccoeff)) &
          num = num + 1d-6*sizeof(Atom(ia)%Ccoeff)

        ! CcoeffJ
        if (allocated(Atom(ia)%CcoeffJ)) &
          num = num + 1d-6*sizeof(Atom(ia)%CcoeffJ)

        ! eval
        if (allocated(Atom(ia)%eval)) &
          num = num + 1d-6*sizeof(Atom(ia)%eval)

        ! gk
        if (allocated(Atom(ia)%gk)) &
          num = num + 1d-6*sizeof(Atom(ia)%gk)

        ! evec
        if (allocated(Atom(ia)%evec)) &
          num = num + 1d-6*sizeof(Atom(ia)%evec)

        ! crho
        if (allocated(Atom(ia)%crho)) &
          num = num + 1d-6*sizeof(Atom(ia)%crho)

      end do

      return

      end subroutine cram_atom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the LTElines variable\n
      !!  LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!                 num(double): Memory count
      subroutine cram_ltelines(LTElines,num)

      ! I/O

      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      double precision, intent(out):: num

      ! Local

      integer:: ia


      ! Initialize
      num = 0d0

      ! Lines
      do ia=lbound(LTElines,1),ubound(LTElines,1)

        ! Non arrays
        num = num + 1d-6*sizeof(LTElines(ia))

        !
        ! Arrays
        !

        ! damp
        if (allocated(LTElines(ia)%damp)) &
          num = num + 1d-6*sizeof(LTElines(ia)%damp)

        ! broad_args
        if (allocated(LTElines(ia)%broad_args)) &
          num = num + 1d-6*sizeof(LTElines(ia)%broad_args)

        ! nl
        if (allocated(LTElines(ia)%nl)) &
          num = num + 1d-6*sizeof(LTElines(ia)%nl)

        ! nu
        if (allocated(LTElines(ia)%nu)) &
          num = num + 1d-6*sizeof(LTElines(ia)%nu)

        ! n
        if (allocated(LTElines(ia)%n)) &
          num = num + 1d-6*sizeof(LTElines(ia)%n)

      end do

      end subroutine cram_ltelines

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Mol variable\n
      !!  Mol(Mol_class(:)): Structures with molecular data\n
      !!        num(double): Memory count
      subroutine cram_mol(Mol,num)

      ! I/O

      type(Mol_class), dimension(:), allocatable, intent(in):: Mol
      double precision, intent(out):: num

      ! Local

      integer:: im,ios


      ! Initialize
      num = 0d0

      ! Molecules
      do im=lbound(Mol,1),ubound(Mol,1)

        ! Non arrays
        num = num + 1d-6*sizeof(Mol(im))

        !
        ! Arrays
        !

        ! catom
        if (allocated(Mol(im)%catom)) then
          do ios=lbound(Mol(im)%catom,1),ubound(Mol(im)%catom,1)
            num = num + 1d-6*sizeof(Mol(im)%catom(ios)%s)
          end do
        end if

        ! Molecule
        if (allocated(Mol(im)%Molecule)) &
          num = num + 1d-6*sizeof(Mol(im)%Molecule)

        ! natom
        if (allocated(Mol(im)%natom)) &
          num = num + 1d-6*sizeof(Mol(im)%natom)

        ! iatom
        if (allocated(Mol(im)%iatom)) &
          num = num + 1d-6*sizeof(Mol(im)%iatom)

        ! pfcoeff
        if (allocated(Mol(im)%pfcoeff)) &
          num = num + 1d-6*sizeof(Mol(im)%pfcoeff)

        ! eqcoeff
        if (allocated(Mol(im)%eqcoeff)) &
          num = num + 1d-6*sizeof(Mol(im)%eqcoeff)

        ! pf
        if (allocated(Mol(im)%pf)) &
          num = num + 1d-6*sizeof(Mol(im)%pf)

        ! eq
        if (allocated(Mol(im)%eq)) &
          num = num + 1d-6*sizeof(Mol(im)%eq)

        ! n
        if (allocated(Mol(im)%n)) &
          num = num + 1d-6*sizeof(Mol(im)%n)

      end do

      end subroutine cram_mol

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Atmo variable\n
      !!  Atmo(Atmo_class): Structure with atmospheric data\n
      !!       num(double): Memory count
      subroutine cram_atmo(Atmo,num)

      ! I/O

      type(Atmo_class), intent(in):: Atmo
      double precision, intent(out):: num


      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(Atmo)

      !
      ! Arrays
      !

      ! ele
      if (allocated(Atmo%ele)) &
        num = num + 1d-6*sizeof(Atmo%ele)

      ! z
      if (associated(Atmo%z)) &
        num = num + 1d-6*sizeof(Atmo%z)

      ! T
      if (associated(Atmo%T)) &
        num = num + 1d-6*sizeof(Atmo%T)

      ! vmi
      if (associated(Atmo%vmi)) &
        num = num + 1d-6*sizeof(Atmo%vmi)

      ! vx
      if (associated(Atmo%vx)) &
        num = num + 1d-6*sizeof(Atmo%vx)

      ! vy
      if (associated(Atmo%vy)) &
        num = num + 1d-6*sizeof(Atmo%vy)

      ! vz
      if (associated(Atmo%vz)) &
        num = num + 1d-6*sizeof(Atmo%vz)

      ! The vxa, vya, and vza pointer point either to
      ! vx, vy, and vz or to a common array of zeros

      ! Bx
      if (associated(Atmo%Bx)) &
        num = num + 1d-6*sizeof(Atmo%Bx)

      ! By
      if (associated(Atmo%By)) &
        num = num + 1d-6*sizeof(Atmo%By)

      ! Bz
      if (associated(Atmo%Bz)) &
        num = num + 1d-6*sizeof(Atmo%Bz)

      ! nHT
      if (allocated(Atmo%nHT)) &
        num = num + 1d-6*sizeof(Atmo%nHT)

      ! nHm
      if (allocated(Atmo%nHm)) &
        num = num + 1d-6*sizeof(Atmo%nHm)

      ! Pg
      if (allocated(Atmo%Pg)) &
        num = num + 1d-6*sizeof(Atmo%Pg)

      ! rho
      if (allocated(Atmo%rho)) &
        num = num + 1d-6*sizeof(Atmo%rho)

      ! Pe
      if (allocated(Atmo%Pe)) &
        num = num + 1d-6*sizeof(Atmo%Pe)

      ! zalt
      if (allocated(Atmo%zalt)) &
        num = num + 1d-6*sizeof(Atmo%zalt)

      ! nHa
      if (allocated(Atmo%nHa)) &
        num = num + 1d-6*sizeof(Atmo%nHa)

      ! pT
      if (allocated(Atmo%pT)) &
        num = num + 1d-6*sizeof(Atmo%pT)

      ! abund
      if (allocated(Atmo%abund)) &
        num = num + 1d-6*sizeof(Atmo%abund)

      ! JKQin
      if (allocated(Atmo%JKQin)) &
        num = num + 1d-6*sizeof(Atmo%JKQin)

      ! ne
      if (allocated(Atmo%ne)) &
        num = num + 1d-6*sizeof(Atmo%ne)

      ! vlos
      if (allocated(Atmo%vlos)) &
        num = num + 1d-6*sizeof(Atmo%vlos)

      ! vpos
      if (allocated(Atmo%vpos)) &
        num = num + 1d-6*sizeof(Atmo%vpos)

      ! vphi
      if (allocated(Atmo%vphi)) &
        num = num + 1d-6*sizeof(Atmo%vphi)

      ! chi500
      if (allocated(Atmo%chi500)) &
        num = num + 1d-6*sizeof(Atmo%chi500)

      ! zeros
      if (associated(Atmo%zeros)) &
        num = num + 1d-6*sizeof(Atmo%zeros)

      ! nh
      if (allocated(Atmo%nh)) &
        num = num + 1d-6*sizeof(Atmo%nh)

      ! nhe
      if (allocated(Atmo%nhe)) &
        num = num + 1d-6*sizeof(Atmo%nhe)

      end subroutine cram_atmo

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the MPID variable\n
      !!  MPID(MPI_class): Structure with MPI data\n
      !!      num(double): Memory count
      subroutine cram_mpi(MPID,num)

      ! I/O

      type(MPI_class), intent(in):: MPID
      double precision, intent(out):: num


      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(MPID)

      !
      ! Arrays
      !

      ! nf
      if (allocated(MPID%nf)) &
        num = num + 1d-6*sizeof(MPID%nf)

      ! if0
      if (allocated(MPID%if0)) &
        num = num + 1d-6*sizeof(MPID%if0)

      ! if1
      if (allocated(MPID%if1)) &
        num = num + 1d-6*sizeof(MPID%if1)

      ! size1
      if (allocated(MPID%size1)) &
        num = num + 1d-6*sizeof(MPID%size1)

      ! size2
      if (allocated(MPID%size2)) &
        num = num + 1d-6*sizeof(MPID%size2)

      ! size3
      if (allocated(MPID%size3)) &
        num = num + 1d-6*sizeof(MPID%size3)

      ! size4
      if (allocated(MPID%size4)) &
        num = num + 1d-6*sizeof(MPID%size4)

      ! size5
      if (allocated(MPID%size5)) &
        num = num + 1d-6*sizeof(MPID%size5)

      ! size6
      if (allocated(MPID%size6)) &
        num = num + 1d-6*sizeof(MPID%size6)

      ! size7
      if (allocated(MPID%size7)) &
        num = num + 1d-6*sizeof(MPID%size7)

      ! size8
      if (allocated(MPID%size8)) &
        num = num + 1d-6*sizeof(MPID%size8)

      ! size10
      if (allocated(MPID%size10)) &
        num = num + 1d-6*sizeof(MPID%size10)

      ! sizei0
      if (allocated(MPID%sizei0)) &
        num = num + 1d-6*sizeof(MPID%sizei0)

      ! sizei1
      if (allocated(MPID%sizei1)) &
        num = num + 1d-6*sizeof(MPID%sizei1)

      ! sizei2
      if (allocated(MPID%sizei2)) &
        num = num + 1d-6*sizeof(MPID%sizei2)

      ! sizei3
      if (allocated(MPID%sizei3)) &
        num = num + 1d-6*sizeof(MPID%sizei3)

      ! sizei4
      if (allocated(MPID%sizei4)) &
        num = num + 1d-6*sizeof(MPID%sizei4)

      ! sizei5
      if (allocated(MPID%sizei5)) &
        num = num + 1d-6*sizeof(MPID%sizei5)

      ! sizei6
      if (allocated(MPID%sizei6)) &
        num = num + 1d-6*sizeof(MPID%sizei6)

      ! sizei7
      if (allocated(MPID%sizei7)) &
        num = num + 1d-6*sizeof(MPID%sizei7)

      ! sizei8
      if (allocated(MPID%sizei8)) &
        num = num + 1d-6*sizeof(MPID%sizei8)

      ! sizei9
      if (allocated(MPID%sizei9)) &
        num = num + 1d-6*sizeof(MPID%sizei9)

      ! sizei10
      if (allocated(MPID%sizei10)) &
        num = num + 1d-6*sizeof(MPID%sizei10)

      ! sizei11
      if (allocated(MPID%sizei11)) &
        num = num + 1d-6*sizeof(MPID%sizei11)

      ! sizei12
      if (allocated(MPID%sizei12)) &
        num = num + 1d-6*sizeof(MPID%sizei12)

      ! sizei13
      if (allocated(MPID%sizei13)) &
        num = num + 1d-6*sizeof(MPID%sizei13)

      ! sizei14
      if (allocated(MPID%sizei14)) &
        num = num + 1d-6*sizeof(MPID%sizei14)

      ! sizei4b
      if (allocated(MPID%sizei4b)) &
        num = num + 1d-6*sizeof(MPID%sizei4b)

      ! ltslave
      if (allocated(MPID%ltslave)) &
        num = num + 1d-6*sizeof(MPID%ltslave)

      ! inf
      if (allocated(MPID%inf)) &
        num = num + 1d-6*sizeof(MPID%inf)

      ! iif0
      if (allocated(MPID%iif0)) &
        num = num + 1d-6*sizeof(MPID%iif0)

      ! iif1
      if (allocated(MPID%iif1)) &
        num = num + 1d-6*sizeof(MPID%iif1)

      ! CMf0
      if (allocated(MPID%CMf0)) &
        num = num + 1d-6*sizeof(MPID%CMf0)

      ! CMf1
      if (allocated(MPID%CMf1)) &
        num = num + 1d-6*sizeof(MPID%CMf1)

      end subroutine cram_mpi

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Input variable\n
      !!  Input(Input_class): Structure with configuration data\n
      !!         num(double): Memory count
      subroutine cram_input(Input,num)

      ! I/O

      type(Input_class), intent(in):: Input
      double precision, intent(out):: num

      ! Local

      integer:: ios

      double precision:: num2


      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(Input)

      !
      ! Arrays
      !

      ! lim_stk
      if (allocated(Input%lim_stk%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_stk%sindx)
      if (allocated(Input%lim_stk%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_stk%nbuff)
      if (allocated(Input%lim_stk%indx)) &
        num = num + 1d-6*sizeof(Input%lim_stk%indx)
      if (allocated(Input%lim_stk%doub)) &
        num = num + 1d-6*sizeof(Input%lim_stk%doub)

      ! lim_ctr
      if (allocated(Input%lim_ctr%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_ctr%sindx)
      if (allocated(Input%lim_ctr%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_ctr%nbuff)
      if (allocated(Input%lim_ctr%indx)) &
        num = num + 1d-6*sizeof(Input%lim_ctr%indx)
      if (allocated(Input%lim_ctr%doub)) &
        num = num + 1d-6*sizeof(Input%lim_ctr%doub)

      ! lim_tau
      if (allocated(Input%lim_tau%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_tau%sindx)
      if (allocated(Input%lim_tau%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_tau%nbuff)
      if (allocated(Input%lim_tau%indx)) &
        num = num + 1d-6*sizeof(Input%lim_tau%indx)
      if (allocated(Input%lim_tau%doub)) &
        num = num + 1d-6*sizeof(Input%lim_tau%doub)

      ! lim_cols_tt
      if (allocated(Input%lim_cols_tt%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_cols_tt%sindx)
      if (allocated(Input%lim_cols_tt%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_cols_tt%nbuff)
      if (allocated(Input%lim_cols_tt%indx)) &
        num = num + 1d-6*sizeof(Input%lim_cols_tt%indx)
      if (allocated(Input%lim_cols_tt%doub)) &
        num = num + 1d-6*sizeof(Input%lim_cols_tt%doub)

      ! lim_cols_ll
      if (allocated(Input%lim_cols_ll%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_cols_ll%sindx)
      if (allocated(Input%lim_cols_ll%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_cols_ll%nbuff)
      if (allocated(Input%lim_cols_ll%indx)) &
        num = num + 1d-6*sizeof(Input%lim_cols_ll%indx)
      if (allocated(Input%lim_cols_ll%doub)) &
        num = num + 1d-6*sizeof(Input%lim_cols_ll%doub)

      ! lim_damp
      if (allocated(Input%lim_damp%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_damp%sindx)
      if (allocated(Input%lim_damp%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_damp%nbuff)
      if (allocated(Input%lim_damp%indx)) &
        num = num + 1d-6*sizeof(Input%lim_damp%indx)
      if (allocated(Input%lim_damp%doub)) &
        num = num + 1d-6*sizeof(Input%lim_damp%doub)

      ! lim_back
      if (allocated(Input%lim_back%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_back%sindx)
      if (allocated(Input%lim_back%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_back%nbuff)
      if (allocated(Input%lim_back%indx)) &
        num = num + 1d-6*sizeof(Input%lim_back%indx)
      if (allocated(Input%lim_back%doub)) &
        num = num + 1d-6*sizeof(Input%lim_back%doub)

      ! lim_pop
      if (allocated(Input%lim_pop%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_pop%sindx)
      if (allocated(Input%lim_pop%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_pop%nbuff)
      if (allocated(Input%lim_pop%indx)) &
        num = num + 1d-6*sizeof(Input%lim_pop%indx)
      if (allocated(Input%lim_pop%doub)) &
        num = num + 1d-6*sizeof(Input%lim_pop%doub)

      ! lim_atmo
      if (allocated(Input%lim_atmo%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_atmo%sindx)
      if (allocated(Input%lim_atmo%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_atmo%nbuff)
      if (allocated(Input%lim_atmo%indx)) &
        num = num + 1d-6*sizeof(Input%lim_atmo%indx)
      if (allocated(Input%lim_atmo%doub)) &
        num = num + 1d-6*sizeof(Input%lim_atmo%doub)

      ! lim_qel
      if (allocated(Input%lim_qel%sindx)) &
        num = num + 1d-6*sizeof(Input%lim_qel%sindx)
      if (allocated(Input%lim_qel%nbuff)) &
        num = num + 1d-6*sizeof(Input%lim_qel%nbuff)
      if (allocated(Input%lim_qel%indx)) &
        num = num + 1d-6*sizeof(Input%lim_qel%indx)
      if (allocated(Input%lim_qel%doub)) &
        num = num + 1d-6*sizeof(Input%lim_qel%doub)

      ! lim_fwhm
      if (allocated(Input%lim_fwhm)) then
        do ios=lbound(Input%lim_fwhm,1),ubound(Input%lim_fwhm,1)
          num = num + 1d-6*sizeof(Input%lim_fwhm(ios))
          if (allocated(Input%lim_fwhm(ios)%sindx)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%sindx)
          if (allocated(Input%lim_fwhm(ios)%indx1)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%indx1)
          if (allocated(Input%lim_fwhm(ios)%indx2)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%indx2)
          if (allocated(Input%lim_fwhm(ios)%doub)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%doub)
          if (allocated(Input%lim_fwhm(ios)%wave)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%wave)
          if (allocated(Input%lim_fwhm(ios)%kernel)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%kernel)
          if (allocated(Input%lim_fwhm(ios)%idx)) &
            num = num + 1d-6*sizeof(Input%lim_fwhm(ios)%idx)
        end do
      end if

      ! LTEline
      if (allocated(Input%LTEline)) then
        call cram_ltelines(Input%LTEline,num2)
        num = num + num2
      end if

      ! atom
      if (allocated(Input%atom)) then
        do ios=lbound(Input%atom,1),ubound(Input%atom,1)
          num = num + 1d-6*sizeof(Input%atom(ios))
        end do
      end if

      ! popu
      if (allocated(Input%popu)) then
        do ios=lbound(Input%popu,1),ubound(Input%popu,1)
          num = num + 1d-6*sizeof(Input%popu(ios))
        end do
      end if

      ! atomback
      if (allocated(Input%atomback)) then
        do ios=lbound(Input%atomback,1),ubound(Input%atomback,1)
          num = num + 1d-6*sizeof(Input%atomback(ios))
        end do
      end if

      ! popuback
      if (allocated(Input%popuback)) then
        do ios=lbound(Input%popuback,1),ubound(Input%popuback,1)
          num = num + 1d-6*sizeof(Input%popuback(ios))
        end do
      end if

      ! mol
      if (allocated(Input%mol)) then
        do ios=lbound(Input%mol,1),ubound(Input%mol,1)
          num = num + 1d-6*sizeof(Input%mol(ios))
        end do
      end if

      ! kurucz
      if (allocated(Input%kurucz)) then
        do ios=lbound(Input%kurucz,1),ubound(Input%kurucz,1)
          num = num + 1d-6*sizeof(Input%kurucz(ios))
        end do
      end if

      ! waves
      if (allocated(Input%waves)) then
        do ios=lbound(Input%waves,1),ubound(Input%waves,1)
          num = num + 1d-6*sizeof(Input%waves(ios))
        end do
      end if

      ! asym_fil
      if (allocated(Input%asym_fil)) then
        do ios=lbound(Input%asym_fil,1),ubound(Input%asym_fil,1)
          num = num + 1d-6*sizeof(Input%asym_fil(ios))
        end do
      end if

      ! fwhm_fil
      if (allocated(Input%fwhm_fil)) then
        do ios=lbound(Input%fwhm_fil,1),ubound(Input%fwhm_fil,1)
          num = num + 1d-6*sizeof(Input%fwhm_fil(ios))
        end do
      end if

      ! ionf
      if (allocated(Input%ionf)) &
        num = num + 1d-6*sizeof(Input%ionf)

      ! fixp
      if (allocated(Input%fixp)) &
        num = num + 1d-6*sizeof(Input%fixp)

      ! zero_ion
      if (allocated(Input%zero_ion)) &
        num = num + 1d-6*sizeof(Input%zero_ion)

      ! skip_wave
      if (allocated(Input%skip_wave)) &
        num = num + 1d-6*sizeof(Input%skip_wave)

      ! fixplt
      if (allocated(Input%fixplt)) &
        num = num + 1d-6*sizeof(Input%fixplt)

      ! sol_box
      if (allocated(Input%sol_box)) &
        num = num + 1d-6*sizeof(Input%sol_box)

      ! Kcut_input
      if (allocated(Input%Kcut_input)) &
        num = num + 1d-6*sizeof(Input%Kcut_input)

      ! excl
      if (allocated(Input%excl)) &
        num = num + 1d-6*sizeof(Input%excl)

      ! L_mu
      if (allocated(Input%L_mu)) &
        num = num + 1d-6*sizeof(Input%L_mu)

      ! L_phi
      if (allocated(Input%L_phi)) &
        num = num + 1d-6*sizeof(Input%L_phi)

      ! asym_num
      if (allocated(Input%asym_num)) &
        num = num + 1d-6*sizeof(Input%asym_num)

      ! Node
      if (allocated(Input%Node)) then
        do ios=lbound(Input%Node,1),ubound(Input%Node,1)
          num = num + 1d-6*sizeof(Input%Node(ios))
          if (allocated(Input%Node(ios)%Tau_Indx)) &
            num = num + 1d-6*sizeof(Input%Node(ios)%Tau_Indx)
          if (allocated(Input%Node(ios)%H)) &
            num = num + 1d-6*sizeof(Input%Node(ios)%H)
          if (allocated(Input%Node(ios)%Var)) &
            num = num + 1d-6*sizeof(Input%Node(ios)%Var)
          if (allocated(Input%Node(ios)%Errors)) &
            num = num + 1d-6*sizeof(Input%Node(ios)%Errors)
          if (allocated(Input%Node(ios)%ebound)) &
            num = num + 1d-6*sizeof(Input%Node(ios)%ebound)
        end do
      end if

      ! Nodes_Flags
      if (allocated(Input%Nodes_Flags)) &
        num = num + 1d-6*sizeof(Input%Nodes_Flags)

      ! Nodes_Regul
      if (allocated(Input%Nodes_Regul)) &
        num = num + 1d-6*sizeof(Input%Nodes_Regul)

      ! Node_Type
      if (allocated(Input%Node_Type)) &
        num = num + 1d-6*sizeof(Input%Node_Type)

      ! Num_nodes
      if (allocated(Input%Num_nodes)) &
        num = num + 1d-6*sizeof(Input%Num_nodes)

      ! Indx_regul
      if (allocated(Input%Indx_regul)) &
        num = num + 1d-6*sizeof(Input%Indx_regul)

      ! Scal
      if (allocated(Input%Scal)) &
        num = num + 1d-6*sizeof(Input%Scal)

      ! Perturb
      if (allocated(Input%Perturb)) &
        num = num + 1d-6*sizeof(Input%Perturb)

      ! min_rel_Pert
      if (allocated(Input%min_rel_Pert)) &
        num = num + 1d-6*sizeof(Input%min_rel_Pert)

      ! Regul_weight
      if (allocated(Input%Regul_weight)) &
        num = num + 1d-6*sizeof(Input%Regul_weight)

      ! Atmo_strat_done
      if (allocated(Input%Atmo_strat_done)) &
        num = num + 1d-6*sizeof(Input%Atmo_strat_done)

      ! Weight
      if (allocated(Input%Weight)) &
        num = num + 1d-6*sizeof(Input%Weight)

      ! Atmo_strat
      if (allocated(Input%Atmo_strat)) &
        num = num + 1d-6*sizeof(Input%Atmo_strat)

      ! Weight_Factor
      if (allocated(Input%Weight_Factor)) &
        num = num + 1d-6*sizeof(Input%Weight_Factor)

      ! Sigma_Factor
      if (allocated(Input%Sigma_Factor)) &
        num = num + 1d-6*sizeof(Input%Sigma_Factor)

      end subroutine cram_input

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Geom variable\n
      !!  Geom(Geometry_class): Structure with geometric data\n
      !!           num(double): Memory count
      subroutine cram_geom(Geom,num)

      ! I/O

      type(Geometry_class), intent(in):: Geom
      double precision, intent(out):: num


      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(Geom)

      !
      ! Arrays
      !

      ! ithv
      if (allocated(Geom%ithv)) &
        num = num + 1d-6*sizeof(Geom%ithv)

      ! iphv
      if (allocated(Geom%iphv)) &
        num = num + 1d-6*sizeof(Geom%iphv)

      ! i_geom
      if (allocated(Geom%i_geom)) &
        num = num + 1d-6*sizeof(Geom%i_geom)

      ! i_scatt
      if (allocated(Geom%i_scatt)) &
        num = num + 1d-6*sizeof(Geom%i_scatt)

      ! V_mu
      if (allocated(Geom%V_mu)) &
        num = num + 1d-6*sizeof(Geom%V_mu)

      ! V_mux
      if (allocated(Geom%V_mux)) &
        num = num + 1d-6*sizeof(Geom%V_mux)

      ! V_muy
      if (allocated(Geom%V_muy)) &
        num = num + 1d-6*sizeof(Geom%V_muy)

      ! V_theta
      if (allocated(Geom%V_theta)) &
        num = num + 1d-6*sizeof(Geom%V_theta)

      ! V_phi
      if (allocated(Geom%V_phi)) &
        num = num + 1d-6*sizeof(Geom%V_phi)

      ! W_mu
      if (allocated(Geom%W_mu)) &
        num = num + 1d-6*sizeof(Geom%W_mu)

      ! W_mux
      if (allocated(Geom%W_mux)) &
        num = num + 1d-6*sizeof(Geom%W_mux)

      ! W_mux2
      if (allocated(Geom%W_mux2)) &
        num = num + 1d-6*sizeof(Geom%W_mux2)

      ! L_mu
      if (allocated(Geom%L_mu)) &
        num = num + 1d-6*sizeof(Geom%L_mu)

      ! L_theta
      if (allocated(Geom%L_theta)) &
        num = num + 1d-6*sizeof(Geom%L_theta)

      ! L_phi
      if (allocated(Geom%L_phi)) &
        num = num + 1d-6*sizeof(Geom%L_phi)

      ! V_muAA
      if (allocated(Geom%V_muAA)) &
        num = num + 1d-6*sizeof(Geom%V_muAA)

      ! V_siAA
      if (allocated(Geom%V_siAA)) &
        num = num + 1d-6*sizeof(Geom%V_siAA)

      ! W_muAA
      if (allocated(Geom%W_muAA)) &
        num = num + 1d-6*sizeof(Geom%W_muAA)

      ! V_gauss
      if (allocated(Geom%V_gauss)) &
        num = num + 1d-6*sizeof(Geom%V_gauss)

      ! V_mu_disk
      if (allocated(Geom%V_mu_disk)) &
        num = num + 1d-6*sizeof(Geom%V_mu_disk)

      ! V_CScatt
      if (allocated(Geom%V_CScatt)) &
        num = num + 1d-6*sizeof(Geom%V_CScatt)

      ! V_SScatt
      if (allocated(Geom%V_SScatt)) &
        num = num + 1d-6*sizeof(Geom%V_SScatt)

      ! TSL
      if (associated(Geom%TSL)) &
        num = num + 1d-6*sizeof(Geom%TSL)

      ! TBL
      if (associated(Geom%TBL)) &
        num = num + 1d-6*sizeof(Geom%TBL)

      ! TS
      if (associated(Geom%TS)) &
        num = num + 1d-6*sizeof(Geom%TS)

      ! TB
      if (associated(Geom%TB)) &
        num = num + 1d-6*sizeof(Geom%TB)

      end subroutine cram_geom

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in Bfield_class\n
      !!  Bfield(Bfield_class): Structure with magnetic field data\n
      !!           num(double): Memory count
      subroutine cram_bfield(Bfield,num)

      ! I/O

      type(Bfield_class), intent(in):: Bfield
      double precision, intent(out):: num

      ! Initialize
      num = 0d0

      !
      ! Arrays
      !

      ! Bstrength
      if (allocated(Bfield%Bstrength)) &
        num = num + 1d-6*sizeof(Bfield%Bstrength)

      ! Btheta
      if (allocated(Bfield%Btheta)) &
        num = num + 1d-6*sizeof(Bfield%Btheta)

      ! Bphi
      if (allocated(Bfield%Bphi)) &
        num = num + 1d-6*sizeof(Bfield%Bphi)

      ! Blos
      if (allocated(Bfield%Blos)) &
        num = num + 1d-6*sizeof(Bfield%Blos)

      ! Bpos
      if (allocated(Bfield%Bpos)) &
        num = num + 1d-6*sizeof(Bfield%Bpos)

      ! Azimuth
      if (allocated(Bfield%Azimuth)) &
        num = num + 1d-6*sizeof(Bfield%Azimuth)

      end subroutine cram_bfield

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Frec variable\n
      !!  Frec(Frequency_class): Structure with frequency data\n
      !!            num(double): Memory count
      subroutine cram_frec(Frec,num)

      ! I/O

      type(Frequency_class), intent(in):: Frec
      double precision, intent(out):: num

      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(Frec)

      !
      ! Arrays
      !

      ! Mpif0
      if (allocated(Frec%Mpif0)) &
        num = num + 1d-6*sizeof(Frec%Mpif0)

      ! Mpif1
      if (allocated(Frec%Mpif1)) &
        num = num + 1d-6*sizeof(Frec%Mpif1)

      ! Mlif0
      if (allocated(Frec%Mlif0)) &
        num = num + 1d-6*sizeof(Frec%Mlif0)

      ! Mlif1
      if (allocated(Frec%Mlif1)) &
        num = num + 1d-6*sizeof(Frec%Mlif1)

      ! IW_freq
      if (allocated(Frec%IW_freq)) &
        num = num + 1d-6*sizeof(Frec%IW_freq)

      ! Mntfreq
      if (allocated(Frec%Mntfreq)) &
        num = num + 1d-6*sizeof(Frec%Mntfreq)

      ! Mntfreqi
      if (allocated(Frec%Mntfreqi)) &
        num = num + 1d-6*sizeof(Frec%Mntfreqi)

      ! Mnpfreq
      if (allocated(Frec%Mnpfreq)) &
        num = num + 1d-6*sizeof(Frec%Mnpfreq)

      ! IW_freq_in
      if (allocated(Frec%IW_freq_in)) &
        num = num + 1d-6*sizeof(Frec%IW_freq_in)

      ! mapping
      if (allocated(Frec%mapping)) &
        num = num + 1d-6*sizeof(Frec%mapping)

      ! omega
      if (allocated(Frec%omega)) &
        num = num + 1d-6*sizeof(Frec%omega)

      ! W_freq
      if (allocated(Frec%W_freq)) &
        num = num + 1d-6*sizeof(Frec%W_freq)

      ! omega3
      if (allocated(Frec%omega3)) &
        num = num + 1d-6*sizeof(Frec%omega3)

      ! omega_ou
      if (allocated(Frec%omega_ou)) &
        num = num + 1d-6*sizeof(Frec%omega_ou)

      ! omega3_ou
      if (allocated(Frec%omega3_ou)) &
        num = num + 1d-6*sizeof(Frec%omega3_ou)

      ! exu
      if (associated(Frec%exu)) &
        num = num + 1d-6*sizeof(Frec%exu)

      end subroutine cram_frec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in the Flgsg variable\n
      !!  Flgsg(Fctsg_class): Structure with factorials, signs, and
      !!                      J-symbols\n
      !!         num(double): Memory count
      subroutine cram_flgsg(Flgsg,num)

      ! I/O

      type(Fctsg_class), intent(in):: Flgsg
      double precision, intent(out):: num


      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(Flgsg)

      !
      ! Arrays
      !

      ! flg
      if (allocated(Flgsg%flg)) &
        num = num + 1d-6*sizeof(Flgsg%flg)

      ! sg
      if (allocated(Flgsg%sg)) &
        num = num + 1d-6*sizeof(Flgsg%sg)

      ! J-symbols (approximated just to the number)
      ! Memoization is counter apart now
     !num = num + 8d-6*dble(nJs)

      end subroutine cram_flgsg

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in fudge_class\n
      !!  fudge(fudge_class): Structure with fudge data\n
      !!         num(double): Memory count
      subroutine cram_fudge(fudge,num)

      ! I/O

      type(fudge_class), intent(in):: fudge
      double precision, intent(out):: num

      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(fudge)

      !
      ! Arrays
      !

      ! fudge_v
      if (allocated(fudge%fudge_v)) &
        num = num + 1d-6*sizeof(fudge%fudge_v)

      end subroutine cram_fudge

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated in kurucz_class\n
      !! kurucz(kurucz_class): Structure with Kurucz line data\n
      !!          num(double): Memory count
      subroutine cram_kurucz(kurucz,num)

      ! I/O

      type(kurucz_class), intent(in):: kurucz
      double precision, intent(out):: num

      ! Local

      integer:: ios

      ! Initialize
      num = 0d0

      ! Non arrays
      num = num + 1d-6*sizeof(kurucz)

      !
      ! Arrays
      !

      ! tran
      if (allocated(kurucz%tran)) then
        do ios=lbound(kurucz%tran,1),ubound(kurucz%tran,1)
          num = num + 1d-6*sizeof(kurucz%tran(ios))
        end do
      end if

      end subroutine cram_kurucz

!#####################################################################
!#####################################################################
!#####################################################################

      !> Counts the data for the input frequency data\n
      !!   Red(Red_class): Structure with redistribution data\n
      !!      num(double): Memory count
      subroutine cram_red_frec(Red,num)

      ! I/O

      type(Red_class):: Red
      double precision, intent(out):: num

      ! Local

      integer:: indx,jndx

      ! Initialize
      num = 0d0

      !
      ! Arrays
      !

      ! Count indexing
      if (allocated(Red%izao)) &
        num = num + 1d-6*sizeof(Red%izao)

      ! If there is frequency data
      if (associated(Red%zao)) then

        ! For each index allocated
        do indx=1,Red%nzao

          ! Constants
          num = num + 1d-6*sizeof(Red%zao(indx))

          ! nf
          if (allocated(Red%zao(indx)%nf)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%nf)

          ! Mif0
          if (allocated(Red%zao(indx)%Mif0)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%Mif0)

          ! Mif1
          if (allocated(Red%zao(indx)%Mif1)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%Mif1)

          ! Rif0
          if (allocated(Red%zao(indx)%Rif0)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%Rif0)

          ! Rif1
          if (allocated(Red%zao(indx)%Rif1)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%Rif1)

          ! if0
          if (allocated(Red%zao(indx)%if0)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%if0)

          ! if1
          if (allocated(Red%zao(indx)%if1)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%if1)

          ! If input transition data
          if (associated(Red%zao(indx)%trani)) then

            ! For each input transition
            do jndx=1,size(Red%zao(indx)%trani)

              ! Constants
              num = num + 1d-6*sizeof(Red%zao(indx)%trani(jndx))

              ! Count frequency info
              if (allocated(Red%zao(indx)%trani(jndx)%omega)) &
                num = num + 1d-6*sizeof(Red%zao(indx)% &
                                        trani(jndx)%omega)
              if (allocated(Red%zao(indx)%trani(jndx)%W_freq)) &
                num = num + 1d-6*sizeof(Red%zao(indx)% &
                                        trani(jndx)%W_freq)
              if (allocated(Red%zao(indx)%trani(jndx)%mfreq)) &
                num = num + 1d-6*sizeof(Red%zao(indx)% &
                            trani(jndx)%mfreq)

            end do ! Input transitions

          end if ! Input transition data

          ! Emissivity
          if (allocated(Red%zao(indx)%eps20)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%eps20)
          if (allocated(Red%zao(indx)%eps21)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%eps21)
          if (allocated(Red%zao(indx)%eps22)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%eps22)
          if (allocated(Red%zao(indx)%eps23)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%eps23)
          if (allocated(Red%zao(indx)%rpf)) &
            num = num + 1d-6*sizeof(Red%zao(indx)%rpf)

        end do ! Indexes

      end if ! There is PRD data

      return

      end subroutine cram_red_frec

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count data in redistribution\n
      !!   Red(Red_class): Structure with redistribution data\n
      !!      num(double): Memory count
      subroutine cram_red_warr(Red,num)

      ! I/O

      type(Red_class):: Red
      double precision, intent(out):: num

      ! Local

      integer:: indx,jndx

      ! Initialize
      num = 0d0

      !
      ! Arrays
      !

      ! If there is redistriution data
      if (associated(Red%rzao)) then

        ! For each index allocated
        do indx=1,Red%nzao

          ! If input transition data
          if (associated(Red%rzao(indx)%trani)) then

            ! For each input transition
            do jndx=1,size(Red%rzao(indx)%trani)

              ! Constants
              num = num + 1d-6*sizeof(Red%Rzao(indx)%trani(jndx))

              ! Count PRD info
              if (allocated(Red%Rzao(indx)%trani(jndx)%Pwarr2)) &
                num = num + 1d-6*sizeof(Red%Rzao(indx)% &
                                        trani(jndx)%PWarr2)
              if (allocated(Red%Rzao(indx)%trani(jndx)%IWarr2)) &
                num = num + 1d-6*sizeof(Red%Rzao(indx)% &
                                        trani(jndx)%IWarr2)
              if (allocated(Red%Rzao(indx)%trani(jndx)%iPPRD)) &
                num = num + 1d-6*sizeof(Red%Rzao(indx)% &
                            trani(jndx)%iPPRD)

            end do ! Input transitions

          end if ! Input transition data

        end do ! Indexes

      end if ! There is redistribution data

      return

      end subroutine cram_red_warr

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count the data in normalizations for 1st order problem\n
      !!   Red(Red_class): Structure with redistribution data\n
      !!      num(double): Memory count
      subroutine cram_red_norm(Red,num)

      ! I/O

      type(Red_class):: Red
      double precision, intent(out):: num

      ! Local

      integer:: indx


      ! Initialize
      num = 0d0

      !
      ! Arrays
      !

      ! Count indexing
      if (allocated(Red%idzao)) &
        num = num + 1d-6*sizeof(Red%idzao)

      ! If there is Norm data
      if (associated(Red%dzao)) then

        ! For each index allocated
        do indx=1,Red%ndzao

          ! ADd content
          num = num + 1d-6*sizeof(Red%dzao(indx))

          ! Norms
          if (allocated(Red%dzao(indx)%Norm)) &
            num = num + 1d-6*sizeof(Red%dzao(indx)%Norm)
          if (allocated(Red%dzao(indx)%p)) &
            num = num + 1d-6*sizeof(Red%dzao(indx)%p)
          if (allocated(Red%dzao(indx)%cp)) &
            num = num + 1d-6*sizeof(Red%dzao(indx)%cp)

        end do ! Indexes

      end if ! There is Norm data

      return

      end subroutine cram_red_norm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count data in 1st order PRD norms\n
      !!   Red(Red_class): Structure with redistribution data\n
      !!      num(double): Memory count
      subroutine cram_red_1stord(Red,num)

      ! I/O

      type(Red_class):: Red
      double precision, intent(out):: num

      ! Local

      integer:: indx


      ! Initialize
      num = 0d0

      !
      ! Arrays
      !

      ! If there is Norm data
      if (associated(Red%pzao)) then

        ! For each index allocated
        do indx=1,Red%nzao

          ! ADd content
          num = num + 1d-6*sizeof(Red%pzao(indx))

          ! Norms
          if (allocated(Red%pzao(indx)%Norm)) &
            num = num + 1d-6*sizeof(Red%pzao(indx)%Norm)
          if (allocated(Red%pzao(indx)%p)) &
            num = num + 1d-6*sizeof(Red%pzao(indx)%p)
          if (allocated(Red%pzao(indx)%cp)) &
            num = num + 1d-6*sizeof(Red%pzao(indx)%cp)

        end do ! Indexes

      end if ! There is Norm data

      return

      end subroutine cram_red_1stord

!#####################################################################
!#####################################################################
!#####################################################################

      !> Resets the redistribution function structure\n
      !!   Red(Red_class): Structure with redistribution data\n
      !!      num(double): Memory count
      subroutine cram_red(Red,num)

      ! I/O

      type(Red_class):: Red
      double precision, intent(out):: num

      ! Local

      double precision:: num1


      ! Initialize
      num = 8d-6

      !
      ! Redistribution
      call cram_red_frec(Red,num1)
      num = num + num1
      call cram_red_warr(Red,num1)
      num = num + num1

      !
      ! Norm
      call cram_red_norm(Red,num1)
      num = num + num1
      call cram_red_1stord(Red,num1)
      num = num + num1

      return

      end subroutine cram_red

!#####################################################################
!#####################################################################
!#####################################################################

      !> Count memory allocated at hanle entry\n
      !!        Atom(Atom_class(:)): Structures with atomic data\n
      !! LTElines(LTEline_class(:)): Structures with LTE line data\n
      !!           Atmo(Atmo_class): Structure with atmospheric data\n
      !!       Bstrength(double(:)): Magnetic field strength\n
      !!              ndir(integer): Number of directions\n
      !!               pol(logical): If computing polarization\n
      !!                num(double): Estimated memory
      subroutine cram_estimate_norm(Atom,LTElines,Atmo, &
                                    Bstrength,ndir,pol,num)

      ! I/O

      type(Atom_class), dimension(:), allocatable, intent(in):: Atom
      type(LTEline_class), dimension(:), &
                           allocatable, intent(in):: LTElines
      type(Atmo_class), intent(in):: Atmo
      logical, intent(in):: pol
      integer, intent(in):: ndir
      double precision, dimension(:), intent(in):: Bstrength
      double precision, intent(out):: num

      ! Local

      logical:: lvel, field

      integer:: idir,iz,ia,jtran

      double precision:: vel


      ! Initialize
      lvel = .False.
      field = .False.

      num = 0d0

      ! Directions
      do idir=1,ndir

        ! Height
        do iz=Rz0,Rz1

          ! Check velocity if dynamic
          if (dyn) then
            vel = sqrt(Atmo%vx(iz)*Atmo%vx(iz) + &
                       Atmo%vy(iz)*Atmo%vy(iz) + &
                       Atmo%vz(iz)*Atmo%vz(iz))
            lvel = vel.gt.TINYVEL
          end if

          ! If no velocity, skip extra directions
          if (.not.lvel.and.idir.gt.1) cycle

          ! If polarized
          if (pol) then

            ! Magnetic?
            field = Bstrength(iz).gt.TINYB

            ! For each active atom
            do ia=1,nA

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! Skip absent
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! Add logical
                num = num + 4d-6

                ! Magnetic
                if (field) then

                  ! Add to ram
                  num = num + 8d-6*Atom(ia)%trano(jtran)%ncomB

                ! Non-magnetic
                else

                  ! Add to ram
                  num = num + 8d-6*Atom(ia)%trano(jtran)%ncomNB

                end if

              end do ! Transitions
            end do ! Active atoms

            ! Skip if no LTE lines
            if (.not.allocated(LTElines)) cycle

            ! For each LTE line
            do ia=1,size(LTElines)

              ! Skip absent
              if (LTElines(ia)%absent) cycle

              ! Skip too high z
              if (iz.lt.LTElines(ia)%Rz0) cycle

              ! Add logical
              num = num + 4d-6

            end do ! LTE lines

          ! No polarized
          else

            ! For each active atom
            do ia=1,nA

              ! For each transition
              do jtran=1,Atom(ia)%ntran

                ! Skip absent
                if (Atom(ia)%fflag(jtran)%absent) cycle

                ! Add logical
                num = num + 4d-6

                ! Add to ram
                num = num + 8d-6*Atom(ia)%nftran

              end do ! Transitions
            end do ! Active atoms

            ! Skip if no LTE lines
            if (.not.allocated(LTElines)) cycle

            ! For each LTE line
            do ia=1,size(LTElines)

              ! Skip absent
              if (LTElines(ia)%absent) cycle

              ! Skip too high z
              if (iz.lt.LTElines(ia)%Rz0) cycle

              ! Add logical
              num = num + 4d-6

            end do ! LTE lines

          end if ! Polarized

        end do ! Height
      end do ! Directions

      end subroutine cram_estimate_norm

!#####################################################################
!#####################################################################
!#####################################################################

      !> Compute the amount of RAM that is already reserved\n
      !!  MPID(MPI_class): Structure with MPI data
      double precision function cram_add(x)

      ! I/O

      integer, intent(in):: x

      integer:: y


      ! Calculate
      cram_add = PRAMc + VRAMc + WRAMc + RRAMc + &
                 BRAMc + MRAMc + TRAMc + ORAMc + &
                 FRAMc + ERAMc + SRAMc + DRAMc + &
                 DRAM2c

      ! Return
      return

      y = x

      end function cram_add

!#####################################################################
!#####################################################################
!#####################################################################

      end module cram_mod
