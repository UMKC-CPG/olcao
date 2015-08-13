module O_ReadDataSubs

   implicit none

   private
   public :: readLabel, readAndCheckLabel, readData

   integer :: readUnit
   integer :: writeUnit

   interface readData
      module procedure readDouble, read2Double, read3Double, readDoubleArray,&
            & readDoubleMatrix, readIntDouble, readInt, read3Int, readIntArray,&
            & readChar
   end interface readData

   contains


subroutine readLabel(readUnit,writeUnit)

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Local variables.
   character*70 :: foundLabel ! The label that is found by reading.
   character*6  :: formatString ! Used to print the label that was read.

   ! Read and trim the foundLabel that was read in.
   read (readUnit,fmt="(a70)") foundLabel
   foundLabel = trim(foundLabel)

   ! Create a format string for printing the foundLabel.
   write (formatString,fmt="(a2,i2.2,a1)") '(a',len_trim(foundLabel),')'
   write (writeUnit,fmt=formatString) foundLabel
   call flush (writeUnit)

end subroutine readLabel



subroutine readAndCheckLabel(readUnit,writeUnit,length,lookLabel)

   implicit none

   ! Passed parameters
   integer                :: length    ! Length of label.
   character (LEN=length) :: lookLabel ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Local variables.
   character (LEN=length) :: foundLabel   ! The label that is found by reading.
   character*6            :: formatString ! Used to print label that was read.

   ! Read and trim the foundLabel that was read in.
   read (readUnit,*) foundLabel

   ! Create a format string for printing the foundLabel.
   write (formatString,fmt="(a2,i2.2,a1)") '(a',length,')'
   write (writeUnit,fmt=formatString) foundLabel
   call flush (writeUnit)

   ! Compare the foundLabel and the lookLabel if they don't match, then stop.
   if (foundLabel /= lookLabel) then
      write (writeUnit,*) looklabel," label not found."
      stop
   endif

end subroutine readAndCheckLabel

subroutine readDouble (readUnit,writeUnit,doubleVar,length,lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   real (kind=double)     :: doubleVar ! The double that we wish to read.
   integer                :: length    ! Label length.
   character (LEN=length) :: lookLabel ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) doubleVar
   write (writeUnit,fmt="(e18.8)") doubleVar
   call flush (writeUnit)

end subroutine readDouble

subroutine read2Double (readUnit,writeUnit,doubleVar1,doubleVar2,length,&
      & lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   real (kind=double)     :: doubleVar1 ! The double that we wish to read.
   real (kind=double)     :: doubleVar2 ! The double that we wish to read.
   integer                :: length     ! Label length.
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameters.
   read (readUnit,*) doubleVar1, doubleVar2
   write (writeUnit,fmt="(2e28.8)") doubleVar1, doubleVar2
   call flush (writeUnit)

end subroutine read2Double

subroutine read3Double (readUnit,writeUnit,doubleVar1,doubleVar2,doubleVar3,&
      & length,lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   real (kind=double)     :: doubleVar1 ! The double that we wish to read.
   real (kind=double)     :: doubleVar2 ! The double that we wish to read.
   real (kind=double)     :: doubleVar3 ! The double that we wish to read.
   integer                :: length     ! Label length.
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameters.
   read (readUnit,*) doubleVar1, doubleVar2, doubleVar3
   write (writeUnit,fmt="(3e18.8)") doubleVar1, doubleVar2, doubleVar3
   call flush (writeUnit)

end subroutine read3Double

subroutine readDoubleArray (readUnit,writeUnit,numValues,doubleArray,length,&
      & lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer, intent (in)    :: numValues  ! Number of array values.
   real (kind=double), intent (out),&
         & dimension(numValues) :: doubleArray ! Array to read.
   integer, intent (in)    :: length     ! Label length
   character (LEN=length)  :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Local variables.
   character*9 :: formatString

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) doubleArray(1:numValues)
   write (formatString,fmt="(a1,i2.2,a6)") '(',min(numValues,4),'d18.8)'
   write (writeUnit,fmt=formatString) doubleArray(1:numValues)
   call flush (writeUnit)

end subroutine readDoubleArray


subroutine readDoubleMatrix (readUnit,writeUnit,numValues1,numValues2,&
      & doubleMatrix,length,lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer, intent (in)   :: numValues1 ! Number of matrix dim1 values.
   integer, intent (in)   :: numValues2 ! Number of matrix dim2 values.
   real (kind=double), intent (out), &
         & dimension(numValues1,numValues2) :: doubleMatrix ! Matrix to read.
   integer, intent (in)    :: length     ! Label length
   character (LEN=length)  :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Local variables.
   character*9 :: formatString
   integer :: value2

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) doubleMatrix(1:numValues1,1:numValues2)
   write (formatString,fmt="(a1,i2.2,a6)") '(',min(numValues1,4),'d18.8)'
   do value2 = 1, numValues2
      write (writeUnit,fmt=formatString) doubleMatrix(1:numValues1,value2)
   enddo
   call flush (writeUnit)

end subroutine readDoubleMatrix

subroutine readIntDouble (readUnit,writeUnit,intVar,doubleVar,length,lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer                :: intVar     ! The integer that we wish to read.
   real (kind=double)     :: doubleVar  ! The double that we wish to read.
   integer                :: length     ! Label length.
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameters.
   read (readUnit,*) intVar, doubleVar
   write (writeUnit,fmt="(i5,e28.8)") intVar, doubleVar
   call flush (writeUnit)

end subroutine readIntDouble

subroutine readInt (readUnit,writeUnit,integerVar,length,lookLabel)

   implicit none

   ! Passed parameters
   integer                :: integerVar ! The integer we wish to read.
   integer                :: length     ! Label length.
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) integerVar
   write (writeUnit,fmt="(i15)") integerVar
   call flush (writeUnit)

end subroutine readInt


subroutine read3Int (readUnit,writeUnit,intVar1,intVar2,intVar3,length,&
      & lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer                :: intVar1 ! The int that we wish to read.
   integer                :: intVar2 ! The int that we wish to read.
   integer                :: intVar3 ! The int that we wish to read.
   integer                :: length     ! Label length.
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameters.
   read (readUnit,*) intVar1, intVar2, intVar3
   write (writeUnit,fmt="(3i15)") intVar1, intVar2, intVar3
   call flush (writeUnit)

end subroutine read3Int


subroutine readIntArray (readUnit,writeUnit,numValues,intArray,length,lookLabel)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer                :: numValues ! Number of array values.
   integer, dimension(numValues) :: intArray ! Array to read.
   integer                :: length     ! Label length
   character (LEN=length) :: lookLabel  ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Local variables.
   character*6 :: formatString

   ! Check that the label is present and correct.
   if (length /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,length,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) intArray(1:numValues)
   write (formatString,fmt="(a1,i2.2,a3)") '(',numValues,'i4)'
   write (writeUnit,fmt=formatString) intArray(1:numValues)
   call flush (writeUnit)

end subroutine readIntArray



subroutine readChar (readUnit,writeUnit,varLength,charVar,labelLength,lookLabel)

   implicit none

   ! Passed parameters
   integer                     :: varLength   ! Length of variable to be read.
   character (LEN=varLength)   :: charVar     ! Character variable to be read.
   integer                     :: labelLength ! Label length.
   character (LEN=labelLength) :: lookLabel   ! Label that we are looking for.
   integer, intent(in)    :: readUnit  ! The unit number of the file from which
                                       ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Check that the label is present and correct.
   if (labelLength /= 0) then
      call readAndCheckLabel(readUnit,writeUnit,labelLength,lookLabel)
   endif

   ! Read and regurgitate the input parameter.
   read (readUnit,*) charVar
   write (writeUnit,*) charVar
   call flush (writeUnit)

end subroutine readChar

end module O_ReadDataSubs
