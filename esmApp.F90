program esmApp
  use ESMF
  implicit none

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

  ! Init ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! Set decomposition
  ! The layout is fixes as 3,8 which means 3*8*ntiles = 144 core need to be used
  do n = 1, ntiles
     decomptile(1,n) = 3
     decomptile(2,n) = 8
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

  ! Finalize ESMF
  call ESMF_Finalize()

end program esmApp
