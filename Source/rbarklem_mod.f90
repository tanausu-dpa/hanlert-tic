      !> Reading Barklem broadening parameters
      module rbarklem_mod
!#####################################################################
!############################# HEADER ################################
!#####################################################################
!
!  Authors:
!     Tanaus\'u del Pino Alem\'an (IAC)
!  Start:
!     12/07/2022
!  Last version:
!     25/02/2025 V4.0.1
!
!#####################################################################
!#####################################################################
!
!  Changelog:
!
!     25/02/2025:    V4.0.1 - Bugfix: When checking the background
!                             atoms, the transition data from the
!                             active atom list was used instead (TdPA)
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
!  rBarklem
!    Read and initialize the Barklem broadening data
!
!  getBarklem
!    Setup the Barklem quantities to calculate the Van der Waals
!  contribution to the elastic broadening in a given transition from
!  an atomic model
!
!  getBarklem_line
!    Setup the Barklem quantities to calculate the Van der Waals
!  contribution to the elastic broadening in a given LTE transition
!
!  barklem_inter
!    Get the interpolated Barklem coefficients for a given pair of
!  effective quantum numbers
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
      !!    Input(Input_class): Structure with configuration data\n
      !!   Atom(Atom_class(:)): Structures with atomic data\n
      !!  Atomb(Atom_class(:)): Structures with atomic data for
      !!                        background atoms
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

      ! For each active atom
      do ia=1,nA

        ! For each transition
        do itran=1,Atom(ia)%ntran

          ! If Barklem
          if (Atom(ia)%broad_type(itran).eq.0) then

            ! Flag true and stop searching
            barklem = .True.
            exit

          end if ! Barklem broadening

        end do ! Transitions

        ! If found, stop search
        if (barklem) exit

      end do ! Active atoms

      ! If no Barklem found yet
      if (.not.barklem) then

        ! For each background atom
        do ia=1,nAb

          ! For each transition
          do itran=1,Atomb(ia)%ntran

            ! If Barklem
            if (Atomb(ia)%broad_type(itran).eq.0) then

              ! Flag found and stop search
              barklem = .True.
              exit

            end if ! Barklem broadening

          end do ! Transitions

          ! If found, stop search
          if (barklem) exit

        end do ! Backgrdouna atoms

      end if ! No Barklem found yet

      ! If no Barklem found yet
      if (.not.barklem) then

        ! For each LTE line
        do ia=1,nLTEl

          ! If Barklem
          if (Input%LTEline(ia)%broad_type.eq.0) then

            ! Flag found and stop
            barklem = .True.
            exit

          end if ! Barklem broadening

        end do ! LTE lines

      end if ! No Barklem found yet

      ! If no Barklem found, just leave
      if (.not.barklem) return


      !
      ! Read Barklem SP data
      !

      ! If there is a specified file
      if (Input%bark_sp.ne.'NONE') then

        ! Open file
        open (200,file=trim(Input%bark_sp), &
              status='unknown', iostat=ios, err=1000, &
              access='stream', action='read', &
              form='unformatted')

        ! Get dimensions
        read(200,err=1100) n1sp
        read(200,err=1100) n2sp

        ! Flag customized
        bin = .True.

      ! Hard-code
      else

        ! Open default
        open(200,file=trim(Input%resource)//'spdata.dat',err=1000)

        ! Set dimensions
        n1sp = 21
        n2sp = 18

        ! Flag default
        bin = .False.

      end if ! Specified file

      ! Allocate space for tabulation
      allocate(x1sp(n1sp),x2sp(n2sp))
      allocate(sp1(n2sp,n1sp),sp2(n2sp,n1sp))

      ! If there is a custom file
      if (bin) then

        ! Read tabulation
        read(200,err=1100) x1sp
        read(200,err=1100) x2sp
        read(200,err=1100) sp1
        read(200,err=1100) sp2

      ! Hard-coded default
      else

        ! Efective quantum numbers in the tables
        x1sp = (/ 1d0,1.1d0,1.2d0,1.3d0,1.4d0,1.5d0,1.6d0,1.7d0, &
                  1.8d0,1.9d0,2.0d0,2.1d0,2.2d0,2.3d0,2.4d0, &
                  2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0 /)
        x2sp = (/ 1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0, &
                  2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0, &
                  2.9d0,3.0d0 /)

        ! Read s-p data from Anstee and O'Mara 1995, MNRAS 276,859
        do i1=1,n1sp
          read(200,*,err=1100) sp1(:,i1)
        end do
        do i1=1,n1sp
          read(200,*,err=1100) sp2(:,i1)
        end do

      end if ! Custom or default file

      ! Close file
      close(200)


      !
      ! Read Barklem PD data
      !

      ! If there is a specified file
      if (Input%bark_pd.ne.'NONE') then

        ! Open file
        open (200,file=trim(Input%bark_pd), &
              status='unknown', iostat=ios, err=2000, &
              access='stream', action='read', &
              form='unformatted')

        ! Get dimensions
        read(200,err=2100) n1pd
        read(200,err=2100) n2pd

        ! Flag customized
        bin = .True.

      ! Hard-code
      else

        ! Open default
        open(200,file=trim(Input%resource)//'pddata.dat',err=2000)

        ! Set dimensions
        n1pd = 18
        n2pd = 18

        ! Flag default
        bin = .False.

      end if ! Specified file

      ! Allocate space for tabulation
      allocate(x1pd(n1pd),x2pd(n2pd))
      allocate(pd1(n2pd,n1pd),pd2(n2pd,n1pd))

      ! If there is a custom file
      if (bin) then

        ! Read tabulation
        read(200,err=2100) x1pd
        read(200,err=2100) x2pd
        read(200,err=2100) pd1
        read(200,err=2100) pd2

      ! Hard-coded default
      else

        ! Efective quantum numbers in the tables
        x1pd = (/ 1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0, &
                  2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0, &
                  2.9d0,3.0d0 /)
        x2pd = (/ 2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0, &
                  3.1d0,3.2d0,3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0, &
                  3.9d0,4.0d0 /)

        ! Read p-s data from Barklem and O'Mara 1997, MNRAS 290,102
        do i1=1,n1pd
          read(200,*,err=2100) pd1(:,i1)
        end do
        do i1=1,n1pd
          read(200,*,err=2100) pd2(:,i1)
        end do

      end if ! Custom or default file

      ! Close file
      close(200)


      !
      ! Read Barklem DF data
      !

      ! If there is a specified file
      if (Input%bark_df.ne.'NONE') then

        ! Open file
        open (200,file=trim(Input%bark_df), &
              status='unknown', iostat=ios, err=3000, &
              access='stream', action='read', &
              form='unformatted')

        ! Get dimensions
        read(200,err=3100) n1df
        read(200,err=3100) n2df

        ! Flag customized
        bin = .True.

      ! Hard-code
      else

        ! Open default
        open(200,file=trim(Input%resource)//'dfdata.dat',err=3000)

        ! Set dimensions
        n1df = 18
        n2df = 18

        ! Flag default
        bin = .False.

      end if ! Specified file

      ! Allocate space for tabulation
      allocate(x1df(n1df),x2df(n2df))
      allocate(df1(n2df,n1df),df2(n2df,n1df))

      ! If there is a custom file
      if (bin) then

        ! Read tabulation
        read(200,err=3100) x1df
        read(200,err=3100) x2df
        read(200,err=3100) df1
        read(200,err=3100) df2

      ! Hard-coded default
      else

        ! Efective quantum numbers in the tables
        x1df = (/ 2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0, &
                  3.1d0,3.2d0,3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0, &
                  3.9d0,4.0d0 /)
        x2df = (/ 3.3d0,3.4d0,3.5d0,3.6d0,3.7d0,3.8d0,3.9d0,4.0d0, &
                  4.1d0,4.2d0,4.3d0,4.4d0,4.5d0,4.6d0,4.7d0,4.8d0, &
                  4.9d0,5.0d0 /)

        ! Read d-f data from Anstee, O'Mara, and Ross 1998, MNRAS
        ! 296,1057
        do i1=1,n1df
          read(200,*,err=3100) df1(:,i1)
        end do
        do i1=1,n1df
          read(200,*,err=3100) df2(:,i1)
        end do

      end if ! Custom or default file

      ! Close file
      close(200)

      ! For each active atom
      do ia=1,na

        ! For each transition
        do itran=1,Atom(ia)%ntran

          ! Get terms
          iterm  = Atom(ia)%fst(itran)%iterml
          iterm1 = Atom(ia)%fst(itran)%itermu

          ! Check if doing Barklem
          call getBarklem(Atom(ia),itran,iterm,iterm1, &
                          n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                          n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                          n1df,n2df,x1df,x2df,df1,df2)

        end do ! Transitions
      end do ! Atoms

      ! For each background atom
      do ia=1,nab

        ! For each transition
        do itran=1,Atomb(ia)%ntran

          ! Get terms
          iterm  = Atomb(ia)%fst(itran)%iterml
          iterm1 = Atomb(ia)%fst(itran)%itermu

          ! Check if doing Barklem
          call getBarklem(Atomb(ia),itran,iterm,iterm1, &
                          n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                          n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                          n1df,n2df,x1df,x2df,df1,df2)

        end do ! Transitions
      end do ! Background atoms

      ! For each LTE line
      do ia=1,nLTEl

        ! Check if doing Barklem
        call getBarklem_line(Input%LTEline(ia), &
                             n1sp,n2sp,x1sp,x2sp,sp1,sp2, &
                             n1pd,n2pd,x1pd,x2pd,pd1,pd2, &
                             n1df,n2df,x1df,x2df,df1,df2)

      end do ! LTE lines

      ! Free
      deallocate(x1sp,x2sp,sp1,sp2)
      deallocate(x1pd,x2pd,pd1,pd2)
      deallocate(x1df,x2df,df1,df2)

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

      !> Setup the Barklem quantities to calculate the Van der Waals
      !! contribution to the elastic broadening in a given transition
      !! from an atomic model\n
      !!  Atom(Atom_class): Structure with atomic data\n
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

      ! If Van der Waals broadening is not Barklem, return
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

      !
      ! Find the next continuum
      !

      ! For all terms above the upper term
      do i1=iterm1,Atom%nMulti

        ! If stage larger than the upper term's
        if (Atom%stage(i1).gt.Atom%stage(iterm1)) then

          ! Assume we found the next continuum
          itermc = i1
          exit

        end if ! Stage larger than upper term's

      end do ! All terms above the upper term

      ! If we did not find the continuum
      if (itermc.lt.0.and.Atom%broad_type(itran).ne.2) then

        ! Issue error
        umsg = 'Could not find continuum in atom '// &
               Atom%Element//' and Van der Waals '// &
               'broadening not set to parametric.'
        call aborted
        return

      end if ! Continuum not found

      ! Next ion charge
      Z = Atom%stage(iterm)

      ! Initialize to use as flag
      id = -1
      check = .False.
      msg = 1

      !
      ! Identify the type of transition
      !

      ! s-p
      if (nint(args(1)).eq.0.and.nint(args(3)).eq.1) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 0
      ! p-s
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.0) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 0
      ! p-d
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.2) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 1
      ! d-p
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.1) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 1
      ! d-f
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.3) then
        t1 = iterm1
        i1 = 1
        t2 = iterm
        i2 = 3
        id = 2
      ! f-d
      else if (nint(args(1)).eq.3.and.nint(args(3)).eq.2) then
        t1 = iterm
        i1 = 3
        t2 = iterm1
        i2 = 1
        id = 2
      end if ! Type of transition

      ! No Barklem for ions
      if (Z.gt.1) then
        id = -1
        msg = 2
      end if

      ! If we identified the transition as s-p, p-d, or d-f
      if (id.ge.0) then

        ! Type of message
        msg = 3

        ! Calculate effective quantum numbers
        neff1 = Z*sqrt(aryd/(args(i1+1)*1d-5 - Atom%TRfreq(t1)))
        neff2 = Z*sqrt(aryd/(args(i2+1)*1d-5 - Atom%TRfreq(t2)))

        !
        ! Get Barklem parameters from tables
        !

        ! sp
        if (id.eq.0) then
          call barklem_inter(neff1,neff2,n1sp,n2sp, &
                             x1sp,x2sp,sp1,sp2, &
                             sigma,alpha,check)
        ! pd
        else if (id.eq.1) then
          call barklem_inter(neff1,neff2,n1pd,n2pd, &
                             x1pd,x2pd,pd1,pd2, &
                             sigma,alpha,check)
        ! df
        else if (id.eq.2) then
          call barklem_inter(neff1,neff2,n1df,n2df, &
                             x1df,x2df,df1,df2, &
                             sigma,alpha,check)
        end if ! electron orbital quantum numbers
      end if ! Correct transition identification

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

          ! Error type 1: Wrong parameters
          if (msg.eq.1) &
          write(umsg,'(A,i2,3A)') ' # Wrong parameters for '// &
                                  'Barklem broadening for '// &
                                  'transition ',itran,' of ', &
                                  Atom%Element,' atom, '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'

          ! Error type 2: Transition in an ion
          if (msg.eq.2) &
          write(umsg,'(A,i2,3A)') ' # Barklem broadening only '// &
                                  'valid for neutral ions, '// &
                                  'transition ',itran,' of ', &
                                  Atom%Element,' atom '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'

          ! Error type 3: Values out of the tabulation
          if (msg.eq.3) &
          write(umsg,'(A,i2,3A)') ' # Could not find values '// &
                                  'in the tables for the levels '// &
                                  'of transition ',itran,' of ', &
                                  Atom%Element,' atom, '// &
                                  'switch to Unsold without '// &
                                  'any enhancement'
          ! Verbose error
          call verbose

        end if ! There was an error

        ! Switch to Unsold
        Atom%broad_type(itran) = 1
        Atom%broad_args(:,itran) = (/ 1d0,0d0,1d0,0d0 /)

      end if ! Success or failure with setup of Barklem inputs

      end subroutine getBarklem

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get Barklem data for a LTE line\n
      !> Setup the Barklem quantities to calculate the Van der Waals
      !! contribution to the elastic broadening in a given LTE
      !! transition\n
      !! line(LTEline_class): Structure with LTE line data\n
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


      ! If Van der Waals broadening is not Barklem, return
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

      !
      ! Identify the type of transition
      !

      ! s-p
      if (nint(args(1)).eq.0.and.nint(args(3)).eq.1) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 0
      ! p-s
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.0) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 0
      ! p-d
      else if (nint(args(1)).eq.1.and.nint(args(3)).eq.2) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 1
      ! d-p
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.1) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 1
      ! d-f
      else if (nint(args(1)).eq.2.and.nint(args(3)).eq.3) then
        e1 = line%Eu
        i1 = 1
        e2 = line%El
        i2 = 3
        id = 2
      ! f-d
      else if (nint(args(1)).eq.3.and.nint(args(3)).eq.2) then
        e1 = line%El
        i1 = 3
        e2 = line%Eu
        i2 = 1
        id = 2
      end if ! Type of transition

      ! No Barklem for ions
      if (Z.gt.1) then
        id = -1
        msg = 2
      end if

      ! If we identified the transition as s-p, p-d or d-f
      if (id.ge.0) then

        ! Type of message
        msg = 3

        ! Calculate effective quantum numbers
        neff1 = Z*sqrt(aryd/(args(i1+1)*1d-5 - e1))
        neff2 = Z*sqrt(aryd/(args(i2+1)*1d-5 - e2))

        !
        ! Get Barklem parameters from tables
        !

        ! sp
        if (id.eq.0) then
          call barklem_inter(neff1,neff2,n1sp,n2sp, &
                             x1sp,x2sp,sp1,sp2, &
                             sigma,alpha,check)
        ! pd
        else if (id.eq.1) then
          call barklem_inter(neff1,neff2,n1pd,n2pd, &
                             x1pd,x2pd,pd1,pd2, &
                             sigma,alpha,check)
        ! df
        else if (id.eq.2) then
          call barklem_inter(neff1,neff2,n1df,n2df, &
                             x1df,x2df,df1,df2, &
                             sigma,alpha,check)
        end if ! electron orbital quantum numbers
      end if ! Correct transition identification

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

          ! Error type 1: Wrong parameters
          if (msg.eq.1) &
          write(umsg,'(A)') ' # Wrong parameters for '// &
                            'Barklem broadening for '// &
                            'transition  of '// &
                            line%Element// &
                            ' atom, switch to Unsold without '// &
                            'any enhancement'

          ! Error type 2: Transition in an ion
          if (msg.eq.2) &
          write(umsg,'(A)') ' # Barklem broadening only '// &
                            'valid for neutral ions, '// &
                            'transition of '// &
                            line%Element// &
                            ' atom switch to Unsold without '// &
                            'any enhancement'

          ! Error type 3: Values out of the tabulation
          if (msg.eq.3) &
          write(umsg,'(A)') ' # Could not find values '// &
                            'in the tables for the levels'// &
                            'of transition of '// &
                            line%Element// &
                            ' atom, switch to Unsold without '// &
                            'any enhancement'
          ! Verbose error
          call verbose

        end if ! There was an error

        ! Switch to Unsold
        line%broad_type = 1
        line%broad_args = (/ 1d0,0d0,1d0,0d0 /)

      end if ! Success or failure with setup of Barklem inputs

      end subroutine getBarklem_line

!#####################################################################
!#####################################################################
!#####################################################################

      !> Get the interpolated Barklem coefficients for a given pair of
      !! effective quantum numbers\n
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
