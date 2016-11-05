module O_StringSubs

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! This function will get the leading alphabetical characters in a string and
!   returns them.  The assumption is that the string will be only 3 characters
!   long at the most.
function getLeadChars (label)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy parameters.
   character*70 :: label

   ! Define local and return variables.
   character*4 :: formatString
   character*3 :: getLeadChars
   integer     :: charIndex
   integer     :: firstNonAlphaPos
   integer     :: i

   ! Go through the string from the beginning and find the first character
   !   that is not within the A-Za-z range of the ASCII sequence.  (Note that
   !   a typical string will be of the form "Si1_1" so that the first non-
   !   alphabetical character will always be at an index > 1 (start counting
   !   from 1).)
   firstNonAlphaPos = 0 ! Initialize to avoid compiler warning.
   if (len(label) > 0) then
      do i = 1, len(label)
         charIndex = iachar(label(i:i))
         if ((charIndex < iachar('A')) .or. (charIndex > iachar('z'))) then
            firstNonAlphaPos = i
            exit
         elseif ((charIndex > iachar('Z')) .and. (charIndex < iachar('a'))) then
            firstNonAlphaPos = i
            exit
         endif
      enddo
      write (formatString,fmt="(a2,i1,a1)") '(a',firstNonAlphaPos-1,')'
   else
      stop 'A label has length 0'
   endif

   if (firstNonAlphaPos == 0) then
      stop "No a-z or A-Z characters found in label"
   endif

   write (getLeadChars,fmt=formatString) label
end function getLeadChars

end module O_StringSubs
