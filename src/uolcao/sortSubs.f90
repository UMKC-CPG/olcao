module O_SortSubs

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public
   private :: mergeLists

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine mergeSort (dataIn,dataOut,sortOrderOut,segmentBorders,dataLength)

   ! Import the precision variables.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:) :: dataIn
   real (kind=double), dimension (:) :: dataOut
   integer, dimension (:) :: sortOrderOut
   integer, dimension (:) :: segmentBorders
   integer :: dataLength

   ! Define the local variables
   integer :: i ! Loop index variable.
   integer :: initNumSegments
   integer :: currentNumSegments
   integer :: mergeListDirection ! 0 = merging from In to Out, 1 = Out to In.
   integer, allocatable, dimension (:) :: sortOrderIn

   ! Initialize the merge direction and segment info.
   mergeListDirection = 0
   initNumSegments    = size (segmentBorders) - 1
   currentNumSegments = initNumSegments

   ! Allocate space to hold the initial sorted order.
   allocate (sortOrderIn (dataLength))

   ! Initialize the sorting orders.
   do i = 1, dataLength
      sortOrderIn(i) = i
      sortOrderOut(i) = i
   enddo

   do i = 1, ceiling(log10(real(initNumSegments,double)) / log10(2.0_double))
      if (mergeListDirection == 0) then
         call mergeLists(dataIn,dataOut,sortOrderIn,sortOrderOut,&
               & segmentBorders, currentNumSegments)
      else
         call mergeLists(dataOut,dataIn,sortOrderOut,sortOrderIn,&
               & segmentBorders, currentNumSegments)
      endif

      mergeListDirection = mod(mergeListDirection + 1,2)

      currentNumSegments = currentNumSegments/2 + mod(currentNumSegments,2)
   enddo

   if (mergeListDirection == 0) then

      dataOut(:)      = dataIn(:)
      sortOrderOut(:) = sortOrderIn(:)
   endif

   ! Deallocate local arrays
   deallocate (sortOrderIn)

end subroutine mergeSort

subroutine mergeLists (fromData,toData,fromOrder,toOrder,segmentBorders,&
      & currentNumSegments)

   ! Import the kinds
   use O_Kinds

   ! Make sure no implicit variables are declared
   implicit none

   ! Define the passed dummy variables
   real (kind=double), dimension (:) :: fromData
   real (kind=double), dimension (:) :: toData
   integer, dimension (:) :: fromOrder
   integer, dimension (:) :: toOrder
   integer, dimension (:) :: segmentBorders
   integer :: currentNumSegments

   ! Define the local variables
   integer :: i ! Loop index variable
   integer :: toPos
   integer :: fromPos1
   integer :: fromPos2
   integer :: finalPos1
   integer :: finalPos2

   toPos = 1

   do i = 1, currentNumSegments/2  ! Integer division so 3/2 = 1 and 1/2 = 0.

      fromPos1  = segmentBorders(2*i-1)+1
      fromPos2  = segmentBorders(2*i)+1
      finalPos1 = fromPos2 - 1
      finalPos2 = segmentBorders(2*i+1)

      ! Check for the case where the fromPos1 and fromPos2 are the same (i.e.
      !   one of the two segments is empty.)
      if (fromPos1 == fromPos2) then
         do while (fromPos2 <= finalPos2)
            toData(toPos)  = fromData(fromPos2)
            toOrder(toPos) = fromOrder(fromPos2)
            fromPos2 = fromPos2 + 1
            toPos = toPos + 1
         enddo
      else

         do while (toPos <= finalPos2)

            if (fromData(fromPos1) < fromData(fromPos2)) then

               toData(toPos)  = fromData(fromPos1)
               toOrder(toPos) = fromOrder(fromPos1)
               toPos = toPos + 1

               fromPos1 = fromPos1 + 1
               if (fromPos1 > finalPos1) then

                  do while (fromPos2 <= finalPos2)

                     toData(toPos)  = fromData(fromPos2)
                     toOrder(toPos) = fromOrder(fromPos2)
                     fromPos2 = fromPos2 + 1
                     toPos = toPos + 1
                  enddo
               endif
            else

               toData(toPos)  = fromData(fromPos2)
               toOrder(toPos) = fromOrder(fromPos2)
               toPos = toPos + 1

               fromPos2 = fromPos2 + 1
               if (fromPos2 > finalPos2) then

                  do while (fromPos1 <= finalPos1)

                     toData(toPos)  = fromData(fromPos1)
                     toOrder(toPos) = fromOrder(fromPos1)
                     fromPos1 = fromPos1 + 1
                     toPos = toPos + 1
                  enddo
               endif
            endif
         enddo
      endif
      ! Adjust the segment border values
      segmentBorders(i+1) = finalPos2
   enddo

   ! Copy the last odd segment at the end.
   if (mod(currentNumSegments,2) .gt. 0) then
      do i = toPos, segmentBorders(currentNumSegments+1)
         toData(i)  = fromData(i)
         toOrder(i) = fromOrder(i)
      enddo
      segmentBorders(currentNumSegments/2+2) = &
            & segmentBorders(currentNumSegments+1)
   endif

end subroutine mergeLists

end module O_SortSubs
