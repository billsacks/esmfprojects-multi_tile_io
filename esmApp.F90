module test_writer
  use ESMF
  implicit none
  private

  public :: write_singletile
  public :: write_multitile

contains
  subroutine write_singletile()
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

    ! Note 10 DEs, in agreement with write_multitile
    grid = ESMF_GridCreateNoPeriDim( &
         maxIndex = [10, 20], &
         regDecomp = [2, 5], &
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
    call ESMF_FieldWrite(field, fileName='dummy_singletile.nc', variableName='dummy', overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  end subroutine write_singletile

  subroutine write_multitile()
    type(ESMF_Grid) :: grid
    integer :: rc
    character(len=255), parameter :: mosaic_file = 'data/C96_mosaic.nc'
    character(len=255), parameter :: input_dir = 'data/'
    integer, parameter :: ntiles = 6
    integer :: decomptile(2,ntiles)
    integer :: n
    type(ESMF_ArraySpec) :: arraySpec
    type(ESMF_Field) :: field
    type(ESMF_Decomp_Flag) :: decompflagPTile(2,ntiles)

    ! Set decomposition
    ! Tiles 1, 3 and 5 each have one DE; tile 2 has 2 DEs along dimension 1; tile 4 has 3
    ! DEs along dimension 1; tile 6 has 2 DEs along dimension 2. The total processor count
    ! should be 10.
    decomptile(:,:) = 1
    decomptile(1,2) = 2
    decomptile(1,4) = 3
    decomptile(2,6) = 2
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

    ! Create field
    field = ESMF_FieldCreate(grid, arraySpec, staggerloc=ESMF_STAGGERLOC_CENTER, &
         indexflag=ESMF_INDEX_GLOBAL, name='dummy', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Fill field
    call ESMF_FieldFill(field, dataFillScheme='sincos', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Try to write field, this is failing !!!
    !call ESMF_FieldWrite(field, fileName='dummy.nc', variableName='dummy', overwrite=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

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

  call write_singletile()
  call write_multitile()

  ! Finalize ESMF
  call ESMF_Finalize()

end program esmApp
