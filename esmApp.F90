module test_writer
  use ESMF
  implicit none
  private

  public :: write_singletile
  public :: write_and_read_multitile

contains
  subroutine write_singletile(decomp_dim1, decomp_dim2, fname)
    integer, intent(in) :: decomp_dim1
    integer, intent(in) :: decomp_dim2
    character(len=*), intent(in) :: fname

    type(ESMF_Grid) :: grid
    integer :: rc
    integer :: i, j, u1, u2
    integer :: lde, ldeCount
    integer :: lbndX(2), ubndX(2), lbndY(2), ubndY(2)
    real :: multiplier
    real(ESMF_KIND_R8), pointer :: coordX(:), coordY(:)
    real(ESMF_KIND_R8), pointer :: dataPtr4d(:,:,:,:)
    type(ESMF_ArraySpec) :: arraySpec
    type(ESMF_ArraySpec) :: arraySpec_w_ungridded
    type(ESMF_Field) :: field
    type(ESMF_Field) :: field_w_ungridded

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
         computationalLBound=lbndX, computationalUBound=ubndX, &
         farrayPtr=coordX, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do i = lbndX(1), ubndX(1)
       coordX(i) = i*10.0
    end do

    call ESMF_GridGetCoord(grid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         computationalLBound=lbndY, computationalUBound=ubndY, &
         farrayPtr=coordY, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    do j = lbndY(1), ubndY(1)
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

    ! ------------------------------------------------------------------------
    ! Create field with two ungridded dimensions
    ! ------------------------------------------------------------------------

    call ESMF_ArraySpecSet(arraySpec_w_ungridded, typekind=ESMF_TYPEKIND_R8, rank=4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create field
    field_w_ungridded = ESMF_FieldCreate(grid, arraySpec_w_ungridded, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='dummy_w_ungridded', &
         ungriddedLBound=[2,15], ungriddedUBound=[4,18], &
         ! 2nd and 4th dimensions are ungridded dimensions
         gridToFieldMap=[1,3], &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Fill field
    call ESMF_FieldGet(field_w_ungridded, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       call ESMF_FieldGet(field_w_ungridded, localDe=lde, farrayPtr=dataPtr4d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       do u1 = 2,4
          do u2 = 15,18
             do i = lbndX(1), ubndX(1)
                do j = lbndY(1), ubndY(1)
                   multiplier = 10.**(u2-15)
                   dataPtr4d(i,u1,j,u2) = u1*multiplier*(coordX(i) - coordY(j))
                end do
             end do
          end do
       end do
    end do

    ! ------------------------------------------------------------------------
    ! Write field with ungridded dimensions
    ! ------------------------------------------------------------------------
    call ESMF_FieldWrite(field_w_ungridded, fileName=fname, variableName='dummy_w_ungridded', overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine write_singletile

  subroutine write_and_read_multitile(decomp_dim1, decomp_dim2, fname)
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
    type(ESMF_ArraySpec) :: arraySpec_w_ungridded
    type(ESMF_Field) :: field_x, field_y, field_x_copy, field_y_copy
    type(ESMF_Field) :: field_w_ungridded
    type(ESMF_Field) :: field_x_read, field_y_read, field_x_copy_read, field_y_copy_read
    type(ESMF_Field) :: field_w_ungridded_read
    type(ESMF_Array) :: array_x
    type(ESMF_Array) :: array_x_read
    type(ESMF_FieldBundle) :: fbundle
    type(ESMF_FieldBundle) :: fbundle_read
    type(ESMF_Decomp_Flag) :: decompflagPTile(2,ntiles)
    integer :: coordDimCount(2)
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:), coordPtrX(:,:), coordPtrY(:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr4d(:,:,:,:)
    integer :: ldeCount, lde
    integer :: i, j, u1, u2
    real :: multiplier

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
    call ESMF_ArraySpecSet(arraySpec_w_ungridded, typekind=ESMF_TYPEKIND_R8, rank=4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create fields
    field_x = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='x', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_x_read = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='x', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y_read = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y_copy = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y_copy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_y_copy_read = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='y_copy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_w_ungridded = ESMF_FieldCreate(grid, arraySpec_w_ungridded, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='w_ungridded', &
         ungriddedLBound=[2,15], ungriddedUBound=[4,18], &
         ! 2nd and 4th dimensions are ungridded dimensions
         gridToFieldMap=[1,3], &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_w_ungridded_read = ESMF_FieldCreate(grid, arraySpec_w_ungridded, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='w_ungridded', &
         ungriddedLBound=[2,15], ungriddedUBound=[4,18], &
         ! 2nd and 4th dimensions are ungridded dimensions
         gridToFieldMap=[1,3], &
         rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Fill fields (this code based on implementation of ESMF_FieldFill, with dimCount = 2,
    ! coordDimCount = [2,2]): Fill the x field with the x coord and the y field with the y
    ! coord.
    call ESMF_FieldGet(field_x, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       ! x coord
       call ESMF_GridGetCoord(grid, coordDim=1, localDe=lde, farrayPtr=coordPtrX, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_FieldGet(field_x, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtrX(:,:)
    end do
    call ESMF_FieldGet(field_y, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       ! y coord
       call ESMF_GridGetCoord(grid, coordDim=2, localDe=lde, farrayPtr=coordPtrY, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_FieldGet(field_y, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtrY(:,:)
       call ESMF_FieldGet(field_y_copy, localDe=lde, farrayPtr=dataPtr, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       dataPtr(:,:) = coordPtrY(:,:)
    end do
    ! Fill field_w_ungridded in a more complex way
    call ESMF_FieldGet(field_w_ungridded, localDeCount=ldeCount)
    ! For now, ldeCount will always be 1, but handle generality
    do lde = 0, ldeCount-1
       call ESMF_GridGetCoord(grid, coordDim=1, localDe=lde, farrayPtr=coordPtrX, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_GridGetCoord(grid, coordDim=2, localDe=lde, farrayPtr=coordPtrY, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       call ESMF_FieldGet(field_w_ungridded, localDe=lde, farrayPtr=dataPtr4d, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
       do u1 = 2,4
          do u2 = 15,18
             do i = lbound(dataPtr4d, 1), ubound(dataPtr4d, 1)
                do j = lbound(dataPtr4d, 3), ubound(dataPtr4d, 3)
                   multiplier = 10.**(u2-15)
                   dataPtr4d(i,u1,j,u2) = u1*multiplier*(coordPtrX(i,j) - coordPtrY(i,j))
                end do
             end do
          end do
       end do
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
    call ESMF_FieldGet(field_x_read, array=array_x_read)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    field_x_copy_read = ESMF_FieldCreate(grid, array_x_read, staggerloc=ESMF_STAGGERLOC_CENTER, &
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
    fbundle_read = ESMF_FieldBundleCreate(name="myfb_read", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_FieldBundleAdd(fbundle_read, [field_x_read, field_y_read, field_x_copy_read], rc=rc)
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

    ! And write the field that has ungridded dimensions
    call ESMF_FieldWrite(field_w_ungridded, fileName=fname, variableName='w_ungridded', &
         overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Read fields
    call ESMF_FieldBundleRead(fbundle_read, fileName=fname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_FieldRead(field_y_copy_read, fileName=fname, variableName='y_copy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_FieldRead(field_w_ungridded_read, fileName=fname, variableName='w_ungridded', rc=rc)

  end subroutine write_and_read_multitile

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

  ! Note 10 DEs, in agreement with write_and_read_multitile
  call write_singletile( &
       decomp_dim1 = 2, &
       decomp_dim2 = 5, &
       fname = 'dummy_singletile.nc')

  ! Tiles 1, 3 and 5 each have one DE; tile 2 has 2 DEs along dimension 1; tile 4 has 3
  ! DEs along dimension 1; tile 6 has 2 DEs along dimension 2. The total processor count
  ! should be 10.
  call write_and_read_multitile( &
       decomp_dim1 = [1,2,1,3,1,1], &
       decomp_dim2 = [1,1,1,1,1,2], &
       fname = 'dummy_multitileA#.nc')

  ! Similar to the last version, but changing which tiles have multiple DEs; in
  ! particular, put multiple on the first tile
  call write_and_read_multitile( &
       decomp_dim1 = [3,1,2,1,1,1], &
       decomp_dim2 = [1,1,1,1,1,2], &
       fname = 'dummy_multitileB#.nc')

  ! Again similar, but swapping dim1 and dim2
  call write_and_read_multitile( &
       decomp_dim1 = [1,1,1,1,1,2], &
       decomp_dim2 = [3,1,2,1,1,1], &
       fname = 'dummy_multitileC#.nc')

  ! Finalize ESMF
  call ESMF_Finalize()

end program esmApp
