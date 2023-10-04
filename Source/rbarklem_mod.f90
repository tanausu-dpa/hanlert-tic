      !> Reading of atmospheric data
      module rbarklem_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC/HAO)
!     Roberto Casini (HAO)
!  Start:
!     07/12/2022
!  Last version:
!     08/07/2023 V3.0.2
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     08/07/2023:    V3.0.2 - Added initialization of Barklem
!                             parameters for LTE lines (TdPA)
!                           - Expect the energies in the model atom
!                             to be in cm^-1 (TdPA)
!                           - Added getBarklem_line (TdPA)
!
!     08/30/2022:    V3.0.1 - Bugfix: There were hard-coded indexes
!                             in the reading loops for default files.
!                             Changed to the appropriate variables
!                             now (TdPA)
!
!     07/12/2022:    V3.0.0 - First version (TdPA)
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
!  rBarklem:
!    Read the Barklem data (or initializes it)
!
!  getBarklem:
!    Get Barklem data of an atomic model
!
!  getBarklem_line:
!    Get Barklem data of an LTE line
!
!  barklem_inter:
!    Interpolation in the Barklem tables
!
!#####################################################################
!#####################################################################
!#####################################################################

      ! Use
      use aborted_mod
      use commons_mod
      use inter_mod
      use parameters_mod , only : ryd , mem
      use types_mod

      contains

!#####################################################################
!#####################################################################
!#####################################################################

      !> Reads a file with the atmospheric data.\n
      !!  Input(Input_class): Structure with settings data\n
      !!    Atom(Atom_class): Structure with the atomic data\n
      !!   Atomb(Atom_class): Structure with the atomic data for
      !!                      background opacities
      subroutine rBarklem(Input,Atom,Atomb)

      ! I/O
      type(Atom_class), dimension(:), intent(inout):: Atom
      type(Atom_class), dimension(:), intent(inout):: Atomb
      type(Input_class), intent(inout):: Input

      ! Local

      logical:: barklem,bin

      integer:: ia,itran,iterm,iterm1,i1,ios

      integer:: n1sp,n2sp,n1pd,n2pd,n1df,n2df
      double precision, dimension(:), allocatable:: x1sp,x2sp
      double precision, dimension(:), allocatable:: x1pd,x2pd
      double precision, dimension(:), allocatable:: x1df,x2df
      double precision, dimension(:,:), allocatable:: sp1,sp2
      double precision, dimension(:,:), allocatable:: pd1,pd2
      double precision, dimension(:,:), allocatable:: df1,df2

      !
      ! Check there is at least a barklem line
      !

      ! Initialize
      barklem = .False.

      ! For each active transition
      do ia=1,nA
        ! For each transition
        do itran=1,Atom(ia)%ntran
          ! If Barklem
          if (Atom(ia)%broad_type(itran).eq.0) then
            barklem = .True.
            exit
          end if
        end do
        if (barklem) exit
      end do

      ! If no Barklem yet
      if (.not.barklem) then

        ! For each passive transition
        do ia=1,nAb
          ! For each transition
          do itran=1,Atomb(ia)%ntran
            ! If Barklem
            if (Atomb(ia)%broad_type(itran).eq.0) then
              barklem = .True.
             exit
            end if
          end do
          if (barklem) exit
        end do

      end if ! No Barklem yet

      ! If no Barklem yet
      if (.not.barklem) then

        ! For each passive transition
        do ia=1,nLTEl
          ! If Barklem
          if (Input%LTEline(ia)%broad_type.eq.0) then
            barklem = .True.
            exit
          end if
        end do

      end if ! No Barklem yet

      ! If no Barklem, just leave
      if (.not.barklem) return

      !
      ! Read Barklem SP data
      !

      ! If there is file
      if (Input%bark_sp.ne.'NONE') then

        open (200,file=trim(Input%bark_sp), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        read(200,err=1100) n1sp
        read(200,err=1100) n2sp
        bin = .True.

      ! Hard-code
      else

        open(200,file=trim(Input%resource)//'spdata.dat',err=1000)

        n1sp = 21
        n2sp = 18
        bin = .False.

      end if

      allocate(x1sp(n1sp),x2sp(n2sp))
      allocate(sp1(n2sp,n1sp),sp2(n2sp,n1sp))

      ! If there is file
      if (bin) then

        read(200,err=1100) x1sp
        read(200,err=1100) x2sp
        read(200,err=1100) sp1
        read(200,err=1100) sp2

      ! Hard-code
      else

        ! Efective quantum numbers in the tables
        x1sp = (/ 1d0,1.1d0,1.2d0,1.3d0,1.4d0,1.5d0,1.6d0,1.7d0, &
                  1.8d0,1.9d0,2.0d0,2.1d0,2.2d0,2.3d0,2.4d0, &
                  2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0 /)
        x2sp = (/ 1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0, &
                  2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0, &
                  2.9d0,3.0d0 /)

        ! s-p data from Anstee and O'Mara 1995, MNRAS 276,859
        do i1=1,n1sp
          read(200,*,err=1100) sp1(:,i1)
        end do
        do i1=1,n1sp
          read(200,*,err=1100) sp2(:,i1)
        end do

      end if

      close(200)

      !
      ! Read Barklem PD data
      !

      ! If there is file
      if (Input%bark_pd.ne.'NONE') then

        open (200,file=trim(Input%bark_pd), &
              status='unknown', iostat=ios, err=2000, &
              access='stream', action='read', &
              form='unformatted')

        read(200,err=2100) n1pd
        read(200,err=2100) n2pd
        bin = .True.

      ! Hard-code
      else

        open(200,file=trim(Input%resource)//'pddata.dat',err=2000)

        n1pd = 18
        n2pd = 18
        bin = .False.

      end if

      allocate(x1pd(n1pd),x2pd(n2pd))
      allocate(pd1(n2pd,n1pd),pd2(n2pd,n1pd))

      ! If there is file
      if (bin) then

        read(200,err=2100) x1pd
        read(200,err=2100) x2pd
        read(200,err=2100) pd1
        read(200,err=2100) pd2

      ! Hard-code
      else

        ! Efective quantum numbers in the tables
        x1pd = (/ 1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0, &
                  2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0, &
                  2.9d0,3.0d0 /)
        x2pd = (/ 2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0, &
                  3.1d0,3.2d0,3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0, &
                  3.9d0,4.0d0 /)

        ! s-p data from Anstee and O'Mara 1995, MNRAS 276,859
        do i1=1,n1pd
          read(200,*,err=2100) pd1(:,i1)
        end do
        do i1=1,n1pd
          read(200,*,err=2100) pd2(:,i1)
        end do

      end if

      close(200)

      !
      ! Read Barklem DF data
      !

      ! If there is file
      if (Input%bark_df.ne.'NONE') then

        open (200,file=trim(Input%bark_df), &
              status='unknown', iostat=ios, err=3000, &
              access='stream', action='read', &
              form='unformatted')

        read(200,err=3100) n1df
        read(200,err=3100) n2df
        bin = .True.

      ! Hard-code
      else

        open(200,file=trim(Input%resource)//'dfdata.dat',err=3000)

        n1df = 18
        n2df = 18
        bin = .False.

      end if

      allocate(x1df(n1df),x2df(n2df))
      allocate(df1(n2df,n1df),df2(n2df,n1df))

      ! If there is file
      if (bin) then

        read(200,err=3100) x1df
        read(200,err=3100) x2df
        read(200,err=3100) df1
        read(200,err=3100) df2

      ! Hard-code
      else

        ! Efective quantum numbers in the tables
        x1df = (/ 2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0, &
                  3.1d0,3.2d0,3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0, &
                  3.9d0,4.0d0 /)
        x2df = (/ 3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0,3.9d0,4.0d0, &
                  4.1d0,4.2d0,4.3d0,4.4d0,4.5d0,4.6d0,4.7d0,4.8d0, &
                  4.9d0,5.0d0 /)

        ! s-p data from Anstee and O'Mara 1995, MNRAS 276,859
        do i1=1,n1df
          read(200,*,err=3100) df1(:,i1)
        end do
        do i1=1,n1df
          read(200,*,err=3100) df2(:,i1)
        end do

      end if

      close(200)

      ! For each active atom
      do ia=1,na

        ! Pair of terms
        do iterm=1,Atom(ia)%nMulti-1
          do iterm1=iterm+1,Atom(ia)%nMulti

            ! Transition
            itran = Atom(ia)%irad(iterm1,iterm)

            ! No transition, continue
            if (itran.le.0) cycle

            ! Check if doing Barklem
            call getBarklem(Atom(ia),itran,iterm,iterm1, &
                            n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                            n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                            n1df,n2df,x1df,x2df,df1,df2)

          end do
        end do
      end do

      ! Do with passives as well
      do ia=1,nab

        ! Pair of terms
        do iterm=1,Atomb(ia)%nMulti-1
          do iterm1=iterm+1,Atomb(ia)%nMulti

            ! Transition
            itran = Atomb(ia)%irad(iterm1,iterm)

            ! No transition, continue
            if (itran.le.0) cycle

            ! Check if doing Barklem
            call getBarklem(Atomb(ia),itran,iterm,iterm1, &
                            n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                            n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                            n1df,n2df,x1df,x2df,df1,df2)

          end do
        end do
      end do

      ! And LTE
      do ia=1,nLTEl

        ! Check if doing Barklem
        call getBarklem_line(Input%LTEline(ia), &
                             n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                             n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                             n1df,n2df,x1df,x2df,df1,df2)
      end do

      return

1000  umsg = 'Error opening Barklem sp file'
      call aborted
      return
1100  umsg = 'Error reading Barklem sp file'
      close(200)
      call aborted
      return
2000  umsg = 'Error opening Barklem pd file'
      call aborted
      return
2100  umsg = 'Error reading Barklem pd file'
      close(200)
      call aborted
      return
3000  umsg = 'Error opening Barklem df file'
      call aborted
      return
3100  umsg = 'Error reading Barklem df file'
      close(200)
      call aborted
      return

      end subroutine rBarklem

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get Barklem data for an atomic model\n
      !!  Atom(Atom_class): Structure with the atomic data\n
      !!    itran(integer): Transition index\n
      !!    iterm(integer): Lower term index\n
      !!   iterm1(integer): Upper term index\n
      !!     n1sp(integer): Dimension axis 1 sp\n
      !!     n2sp(integer): Dimension axis 1 sp\n
      !!   x1sp(dfloat(:)): Axis 1 sp\n
      !!   x2sp(dfloat(:)): Axis 2 sp\n
      !!  sp1(dfloat(:,:)): Table 1 sp\n
      !!  sp2(dfloat(:,:)): Table 2 sp\n
      !!     n1pd(integer): Dimension axis 1 pd\n
      !!     n2pd(integer): Dimension axis 1 pd\n
      !!   x1pd(dfloat(:)): Axis 1 pd\n
      !!   x2pd(dfloat(:)): Axis 2 pd\n
      !!  pd1(dfloat(:,:)): Table 1 pd\n
      !!  pd2(dfloat(:,:)): Table 2 pd\n
      !!     n1df(integer): Dimension axis 1 df\n
      !!     n2df(integer): Dimension axis 1 df\n
      !!   x1df(dfloat(:)): Axis 1 df\n
      !!   x2df(dfloat(:)): Axis 2 df\n
      !!  df1(dfloat(:,:)): Table 1 df\n
      !!  df2(dfloat(:,:)): Table 2 df
      subroutine getBarklem(Atom,itran,iterm,iterm1, &
                            n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                            n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                            n1df,n2df,x1df,x2df,df1,df2)

      ! I/O
      type(Atom_class), intent(inout):: Atom
      integer, intent(in):: itran,iterm,iterm1
      integer, intent(in):: n1sp,n2sp,n1pd,n2pd,n1df,n2df
      double precision, dimension(:), intent(in):: x1sp,x2sp
      double precision, dimension(:), intent(in):: x1pd,x2pd
      double precision, dimension(:), intent(in):: x1df,x2df
      double precision, dimension(:,:), intent(in):: sp1,sp2
      double precision, dimension(:,:), intent(in):: pd1,pd2
      double precision, dimension(:,:), intent(in):: df1,df2

      ! Local

      logical:: check
      integer:: itermc,i1,i2,t1,t2,msg,id,Z
      double precision:: aryd,alpha,sigma,neff1,neff2
      double precision, dimension(4):: args

      ! If not Barklem, return
      if (Atom%broad_type(itran).ne.0) return

      ! Correct rydberg energy for mass shift and calculate
      ! relative velocities for H and He
      aryd = ryd/(1d0 + mem/Atom%rmass)

      ! Take the arguments in the input
      args = Atom%broad_args(:,itran)

      !
      ! Find the first continuum
      !

      ! Initialize the index of the continuum so it can be used as
      ! flag
      itermc = -1

      ! Find the next continuum
      do i1=iterm1,Atom%nMulti

        if (Atom%stage(i1).gt.Atom%stage(iterm1)) then
          itermc = i1
          exit
        end if

      end do

      ! If we did not find the continuum
      if (itermc.lt.0.and.Atom%broad_type(itran).ne.2) then
        umsg = 'Could not find continuum in atom '// &
               Atom%Element//' and Van der Waals '// &
               'broadening not set to parametric.'
        call aborted
        return
      end if

      ! Next ion charge
      Z = Atom%stage(iterm)

      ! Initialize to use as flag
      id = -1
      check = .False.
      msg = 1

      ! Identify the type of transition
      if (nint(args(1)).eq.0.and.nint(args(3)).eq.1) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 0
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.0) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 0
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.2) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 1
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.1) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 1
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.3) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 2
      else if (nint(args(1)).eq.3.and.nint(args(3)).eq.2) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 2
      end if

      if (Z.gt.1) then
        id = -1
        msg = 2
      end if

      ! If we identified the transition as s-p, p-d or d-f
      if (id.ge.0) then

        msg = 3

        ! Calculate effective quantum numbers
        neff1 = Z*sqrt(aryd/(args(i1+1)*1d-5 - Atom%TRfreq(t1)))
        neff2 = Z*sqrt(aryd/(args(i2+1)*1d-5 - Atom%TRfreq(t2)))

        ! Get Barklem parameters from tables
        if (id.eq.0) then
          call barklem_inter(neff1,neff2,n1sp,n2sp, &
                             x1sp,x2sp,sp1,sp2, &
                             sigma,alpha,check)
        else if (id.eq.1) then
          call barklem_inter(neff1,neff2,n1pd,n2pd, &
                             x1pd,x2pd,pd1,pd2, &
                             sigma,alpha,check)
        else if (id.eq.2) then
          call barklem_inter(neff1,neff2,n1df,n2df, &
                             x1df,x2df,df1,df2, &
                             sigma,alpha,check)
        end if

      end if

      ! If sucessfull
      if (id.ge.0.and.check) then

        ! Store in args
        Atom%broad_args(1,itran) = sigma
        Atom%broad_args(2,itran) = alpha

      ! If failed to find the values in the table, or if the
      ! type of transition is not valid, do Unsold
      else

        ! Verbosity of error
        if (gpid.eq.0) then

          if (msg.eq.1) &
          write(umsg,'(A,i2,3A)') ' # Wrong parameters for '// &
                                  'Barklem broadening for '// &
                                  'transition ',itran,' of ', &
                                  Atom%Element,' atom, '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'
          if (msg.eq.2) &
          write(umsg,'(A,i2,3A)') ' # Barklem broadening only '// &
                                  'valid for neutral ions, '// &
                                  'transition ',itran,' of ', &
                                  Atom%Element,' atom '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'
          if (msg.eq.3) &
          write(umsg,'(A,i2,3A)') ' # Could not find values '// &
                                  'in the tables for the levels'// &
                                  'of transition ',itran,' of ', &
                                  Atom%Element,' atom, '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'
          call verbose

        end if

        Atom%broad_type(itran) = 1
        Atom%broad_args(:,itran) = (/ 1d0,0d0,1d0,0d0 /)

      end if ! Barklem inputs

      end subroutine getBarklem

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get Barklem data for a LTE line\n
      !! line(LTEline_class): Structure with the LTE line data\n
      !!       n1sp(integer): Dimension axis 1 sp\n
      !!       n2sp(integer): Dimension axis 1 sp\n
      !!     x1sp(dfloat(:)): Axis 1 sp\n
      !!     x2sp(dfloat(:)): Axis 2 sp\n
      !!    sp1(dfloat(:,:)): Table 1 sp\n
      !!    sp2(dfloat(:,:)): Table 2 sp\n
      !!       n1pd(integer): Dimension axis 1 pd\n
      !!       n2pd(integer): Dimension axis 1 pd\n
      !!     x1pd(dfloat(:)): Axis 1 pd\n
      !!     x2pd(dfloat(:)): Axis 2 pd\n
      !!    pd1(dfloat(:,:)): Table 1 pd\n
      !!    pd2(dfloat(:,:)): Table 2 pd\n
      !!       n1df(integer): Dimension axis 1 df\n
      !!       n2df(integer): Dimension axis 1 df\n
      !!     x1df(dfloat(:)): Axis 1 df\n
      !!     x2df(dfloat(:)): Axis 2 df\n
      !!    df1(dfloat(:,:)): Table 1 df\n
      !!    df2(dfloat(:,:)): Table 2 df
      subroutine getBarklem_line(line, &
                                 n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                                 n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                                 n1df,n2df,x1df,x2df,df1,df2)

      ! I/O
      type(LTEline_class), intent(inout):: line
      integer, intent(in):: n1sp,n2sp,n1pd,n2pd,n1df,n2df
      double precision, dimension(:), intent(in):: x1sp,x2sp
      double precision, dimension(:), intent(in):: x1pd,x2pd
      double precision, dimension(:), intent(in):: x1df,x2df
      double precision, dimension(:,:), intent(in):: sp1,sp2
      double precision, dimension(:,:), intent(in):: pd1,pd2
      double precision, dimension(:,:), intent(in):: df1,df2

      ! Local

      logical:: check
      integer:: i1,i2,msg,id,Z
      double precision:: e1,e2,aryd,alpha,sigma,neff1,neff2
      double precision, dimension(4):: args

      ! If not Barklem, return
      if (line%broad_type.ne.0) return

      ! Correct rydberg energy for mass shift and calculate
      ! relative velocities for H and He
      aryd = ryd/(1d0 + mem/line%rmass)

      ! Take the arguments in the input
      args = line%broad_args

      ! Next ion charge
      Z = line%stage

      ! Initialize to use as flag
      id = -1
      check = .False.
      msg = 1

      ! Identify the type of transition
      if (nint(args(1)).eq.0.and.nint(args(3)).eq.1) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 0
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.0) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 0
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.2) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 1
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.1) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 1
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.3) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 2
      else if (nint(args(1)).eq.3.and.nint(args(3)).eq.2) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 2
      end if

      if (Z.gt.1) then
        id = -1
        msg = 2
      end if

      ! If we identified the transition as s-p, p-d or d-f
      if (id.ge.0) then

        msg = 3

        ! Calculate effective quantum numbers
        neff1 = Z*sqrt(aryd/(args(i1+1)*1d-5 - e1))
        neff2 = Z*sqrt(aryd/(args(i2+1)*1d-5 - e2))

        ! Get Barklem parameters from tables
        if (id.eq.0) then
          call barklem_inter(neff1,neff2,n1sp,n2sp, &
                             x1sp,x2sp,sp1,sp2, &
                             sigma,alpha,check)
        else if (id.eq.1) then
          call barklem_inter(neff1,neff2,n1pd,n2pd, &
                             x1pd,x2pd,pd1,pd2, &
                             sigma,alpha,check)
        else if (id.eq.2) then
          call barklem_inter(neff1,neff2,n1df,n2df, &
                             x1df,x2df,df1,df2, &
                             sigma,alpha,check)
        end if

      end if

      ! If sucessfull
      if (id.ge.0.and.check) then

        ! Store in args
        line%broad_args(1) = sigma
        line%broad_args(2) = alpha

      ! If failed to find the values in the table, or if the
      ! type of transition is not valid, do Unsold
      else

        ! Verbosity of error
        if (gpid.eq.0) then

          if (msg.eq.1) &
          write(umsg,'(A)') ' # Wrong parameters for '// &
                            'Barklem broadening for '// &
                            'transition  of '// &
                            line%Element// &
                            ' atom, switch to Unsold without '// &
                            'any enhancement'
          if (msg.eq.2) &
          write(umsg,'(A)') ' # Barklem broadening only '// &
                            'valid for neutral ions, '// &
                            'transition of '// &
                            line%Element// &
                            ' atom switch to Unsold without '// &
                            'any enhancement'
          if (msg.eq.3) &
          write(umsg,'(A)') ' # Could not find values '// &
                            'in the tables for the levels'// &
                            'of transition of '// &
                            line%Element// &
                            ' atom, switch to Unsold without '// &
                            'any enhancement'
          call verbose

        end if

        line%broad_type = 1
        line%broad_args = (/ 1d0,0d0,1d0,0d0 /)

      end if ! Barklem inputs

      end subroutine getBarklem_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Gets coefficients by Barklem et al. given the effective
      !! principal quantum number\n
      !!             n1(dfloat): Effective quantum number level 1\n
      !!             n2(dfloat): Effective quantum number level 2\n
      !!           nn1(integer): Dimension axis 1\n
      !!           nn2(integer): Dimension axis 2\n
      !!         n1i(dfloat(:)): Axis 1\n
      !!         n2i(dfloat(:)): Axis 2\n
      !!      tab1(dfloat(:,:)): Table 1\n
      !!      tab2(dfloat(:,:)): Table 2\n
      !!              s(dfloat): Parameter to compute broadening\n
      !!              a(dfloat): Parameter to compute broadening\n
      !!            ch(logical): Bool to check success
      subroutine barklem_inter(n1,n2,nn1,nn2,n1i,n2i,tab1,tab2,s,a,ch)

      ! I/O

      logical, intent(out):: ch
      integer, intent(in):: nn1,nn2
      double precision, intent(in):: n1,n2
      double precision, dimension(:), intent(in):: n1i,n2i
      double precision, dimension(:,:), intent(in):: tab1,tab2
      double precision, intent(out):: s,a


      ! Flag failure
      ch = .False.

      ! If we are out of limits, return as failure
      if (n1.lt.n1i(1).or.n1.gt.n1i(nn1).or. &
          n2.lt.n2i(1).or.n2.gt.n2i(nn2)) then
        return
      end if

      ! Interpolate with splines the two coefficients
      call spline_2d(n1i,n2i,tab1,nn1,nn2,n1,n2,s)
      call spline_2d(n1i,n2i,tab2,nn1,nn2,n1,n2,a)

      ! Flag as success
      ch = .True.

      return

      end subroutine barklem_inter

!#####################################################################
!#####################################################################
!#####################################################################

      end module rbarklem_mod
