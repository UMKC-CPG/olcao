module O_WriteDataSubs

private
public :: writeLabel, writeData


interface writeData
   module procedure writeDoubleArray, writeInt, writeIntArray
end interface writeData

contains



subroutine writeLabel (outLabel,fileUnit)

   implicit none

   ! Define passed parameters
   character*(*) :: outLabel ! The label that is to be written
   integer :: fileUnit

   ! Local variables.
   character*6  :: formatString ! Used to print the label that is to be written.

   ! Create a format string for printing the outLabel.
   write (formatString,fmt="(a2,i2.2,a1)") '(a',len_trim(outLabel),')'
   write (fileUnit,fmt=formatString) outLabel
   call flush (fileUnit)

end subroutine writeLabel


subroutine writeDoubleArray (numValues,doubleArray,fileUnit)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer :: numValues ! Number of array values.
   real (kind=double), dimension(numValues) :: doubleArray ! Array to write.
   integer :: fileUnit ! File to write to.

   ! Local variables.
   integer :: i

   ! Write the array with a hard coded scheme of f18.8 numbers across.
   do i = 1, numValues
      write (fileUnit,ADVANCE="NO",fmt='(e18.8)') doubleArray(i)
      if (mod(i,4) == 0) write (fileUnit,*)
   enddo
   if (mod(numValues,4) /= 0) write (fileUnit,*)
   call flush (fileUnit)

end subroutine writeDoubleArray


subroutine writeInt (integerVar,fileUnit)

   ! Use necessary modules.
   use O_Kinds

   implicit none

   ! Passed parameters
   integer :: integerVar ! The integer we wish to write.
   integer :: fileUnit ! File to write to.

   ! Define local variables.
   integer :: fieldWidth
   character*12 :: formatString

   ! Compute the field width.
   fieldWidth = ceiling(log10(abs(real(integerVar,double))+1.0_double))

   ! Prepare the format string
   write (formatString,fmt='(a2,i9.9,a1)') '(i',fieldWidth,')'

   ! Write the integer.
   write (fileUnit,fmt=formatString) integerVar
!   write (fileUnit,fmt="(i15)") integerVar
   call flush (fileUnit)

end subroutine writeInt


subroutine writeIntArray (numValues,intArray,fileUnit)

   use O_Kinds

   implicit none

   ! Passed parameters
   integer :: numValues ! Number of array values.
   integer, dimension(numValues) :: intArray ! Array to write.
   integer :: fileUnit ! File to write to.

   ! Local variables.
   integer :: i

   ! Write the array with the hard coded scheme of i10 across.
   do i = 1, numValues
      write (fileUnit,ADVANCE="NO",fmt='(i8)') intArray(i)
      if (mod(i,10)==0) write (fileUnit,*)
   enddo
   if (mod(numValues,10) /= 0) write (fileUnit,*)
   call flush (fileUnit)

end subroutine writeIntArray


end module O_WriteDataSubs
