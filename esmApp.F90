module test_writer
  use ESMF
  implicit none
  private

  public :: write_singletile
  public :: write_multitile

contains
  subroutine write_singletile(decomp_dim1, decomp_dim2, fname)
    integer, intent(in) :: decomp_dim1
    integer, intent(in) :: decomp_dim2
    character(len=*), intent(in) :: fname

    type(ESMF_Grid) :: grid
    integer :: rc
    integer :: i, j
    integer :: lbnd(2), ubnd(2)
    real(ESMF_KIND_R8), pointer :: coordX(:), coordY(:)
    type(ESMF_ArraySpec) :: arraySpec
    type(ESMF_Field) :: field

    ! ------------------------------------------------------------------------
    ! Create grid
    ! ------------------------------------------------------------------------

    ! The grid creation here follows
    ! http://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/node5.html#SECTION05083100000000000000

    grid = ESMF_GridCreateNoPeriDim( &
         maxIndex = [10, 20], &
         regDecomp = [decomp_dim1, decomp_dim2], &
         coordSys = ESMF_COORDSYS_CART, &
         coordDep1 = [1], &
         coordDep2 = [2], &
         indexflag = ESMF_INDEX_GLOBAL, &
         rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_GridAddCoord(grid, staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    call ESMF_GridGetCoord(grid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do i = lbnd(1), ubnd(1)
       coordX(i) = i*10.0
    end do

    call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbnd, computationalUBound=ubnd, &
         farrayPtr=coordY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do j = lbnd(1), ubnd(1)
       coordY(j) = j*10.0
    end do

    ! ------------------------------------------------------------------------
    ! Create field
    ! ------------------------------------------------------------------------

    ! Set type and rank for ESMF arrayspec
    call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create field
    field = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='dummy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Fill field
    call ESMF_FieldFill(field, dataFillScheme='sincos', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ------------------------------------------------------------------------
    ! Write field
    ! ------------------------------------------------------------------------
    call ESMF_FieldWrite(field, fileName=fname, variableName='dummy', overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine write_singletile

  subroutine write_multitile(decomp_dim1, decomp_dim2, fname)
    ! decomp_dim1 and decomp_dim2 should be of size ntiles, and should have a value for
    ! each tile
    integer, intent(in) :: decomp_dim1(:)
    integer, intent(in) :: decomp_dim2(:)
    character(len=*), intent(in) :: fname

    type(ESMF_Grid) :: grid
    integer :: rc
    character(len=255), parameter :: mosaic_file = 'data/C96_mosaic.nc'
    character(len=255), parameter :: input_dir = 'data/'
    integer, parameter :: ntiles = 6
    integer :: decomptile(2,ntiles)
    integer :: n
    type(ESMF_ArraySpec) :: arraySpec
    type(ESMF_Field) :: field_x, field_y, field_x_copy, field_y_copy
    type(ESMF_Array) :: array_x
    type(ESMF_FieldBundle) :: fbundle
    type(ESMF_Decomp_Flag) :: decompflagPTile(2,ntiles)
    integer :: coordDimCount(2)
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:), coordPtr(:,:)
    integer :: ldeCount, lde

    ! Set decomposition
    decomptile(1,:) = decomp_dim1
    decomptile(2,:) = decomp_dim2
    do n = 1, ntiles
       decompflagPTile(:,n) = (/ ESMF_DECOMP_SYMMEDGEMAX, ESMF_DECOMP_SYMMEDGEMAX /)
    end do

    ! Create CS grid
    grid = ESMF_GridCreateMosaic(filename=trim(mosaic_file), &
         regDecompPTile=decomptile, tileFilePath=trim(input_dir), &
         decompflagPTile=decompflagPTile, &
         staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
         indexflag=ESMF_INDEX_GLOBAL, &
         name='fv3_grid', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Set type and rank for ESMF arrayspec
    call ESMF_ArraySpecSet(arraySpec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create fields
    field_x = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='x', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y_copy = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y_copy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Fill fields (this code based on implementation of ESMF_FieldFill, with dimCount = 2,
    ! coordDimCount = [2,2]): Fill the x field with the x coord and the y field with the y
    ! coord.
    call ESMF_FieldGet(field_x, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       ! x coord
       call ESMF_GridGetCoord(grid, coordDim=1, localDe=lde, farrayPtr=coordPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_FieldGet(field_x, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtr(:,:)
    end do
    call ESMF_FieldGet(field_y, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       ! y coord
       call ESMF_GridGetCoord(grid, coordDim=2, localDe=lde, farrayPtr=coordPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_FieldGet(field_y, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtr(:,:)
       call ESMF_FieldGet(field_y_copy, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtr(:,:)
    end do

    ! Create a copy of the x field that uses the same array, so we can test writing the
    ! same array twice from a single call.
    call ESMF_FieldGet(field_x, array=array_x)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_x_copy = ESMF_FieldCreate(grid, array_x, staggerloc=ESMF_STAGGERLOC_CENTER, &
         name='x_copy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create a FieldBundle with the two fields plus the copy. Note that the copy will
    ! reuse an existing IO decomposition.
    fbundle = ESMF_FieldBundleCreate(name="myfb", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_FieldBundleAdd(fbundle, [field_x, field_y, field_x_copy], rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Write fields
    call ESMF_FieldBundleWrite(fbundle, fileName=fname, overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Also write a new field (that just happens to have the same data as an existing
    ! field), to test that we can write a separate field to an existing file.
    call ESMF_FieldWrite(field_y_copy, fileName=fname, variableName='y_copy', overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine write_multitile

end module test_writer

program esmApp
  use ESMF
  use test_writer
  implicit none

  integer :: rc

  ! Init ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Note 10 DEs, in agreement with write_multitile
  call write_singletile( &
       decomp_dim1 = 2, &
       decomp_dim2 = 5, &
       fname = 'dummy_singletile.nc')

  ! Tiles 1, 3 and 5 each have one DE; tile 2 has 2 DEs along dimension 1; tile 4 has 3
  ! DEs along dimension 1; tile 6 has 2 DEs along dimension 2. The total processor count
  ! should be 10.
  call write_multitile( &
       decomp_dim1 = [1,2,1,3,1,1], &
       decomp_dim2 = [1,1,1,1,1,2], &
       fname = 'dummy_multitileA#.nc')

  ! Similar to the last version, but changing which tiles have multiple DEs; in
  ! particular, put multiple on the first tile
  call write_multitile( &
       decomp_dim1 = [3,1,2,1,1,1], &
       decomp_dim2 = [1,1,1,1,1,2], &
       fname = 'dummy_multitileB#.nc')

  ! Again similar, but swapping dim1 and dim2
  call write_multitile( &
       decomp_dim1 = [1,1,1,1,1,2], &
       decomp_dim2 = [3,1,2,1,1,1], &
       fname = 'dummy_multitileC#.nc')

  ! Finalize ESMF
  call ESMF_Finalize()

end program esmApp
