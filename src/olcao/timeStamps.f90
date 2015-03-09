module O_TimeStamps

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define local variables.
   character(len=8)  :: tempDate
   character(len=10) :: date
   character(len=10) :: tempTime
   character(len=12) :: time
   character(len=51) :: banner
   character(len=51), dimension(30) :: opLabels
   integer :: numOpCodes

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initOperationLabels

   ! Make sure that no funny variables are defined.
   implicit none

   ! Set the number of operation codes.
   numOpCodes = 30

   opLabels(1)  = '***************  Parse Input Files  ***************'
   opLabels(2)  = '***********  Get Implicit Information  ************'
   opLabels(3)  = '************  Lattice Initialization  *************'
   opLabels(4)  = '****  Init. Lattice Vector Finding Structures  ****'
   opLabels(5)  = '******  Get Exchange Correlation Mesh Counts  *****'
   opLabels(6)  = '***********  Initialize Setup HDF5 File  **********'
   opLabels(7)  = '*********  Make Exchange Correlation Mesh  ********'
   opLabels(8)  = '***************  Overlap Integrals  ***************'
   opLabels(9)  = '************  Kinetic Energy Integrals  ***********'
   opLabels(10) = '**********  Nuclear Potential Integrals  **********'
   opLabels(11) = '*********  Electronic Potential Integrals  ********'
   opLabels(12) = '**************  Momentum Integrals  ***************'
   opLabels(13) = '**************  Core Charge Density  **************'
   opLabels(14) = '*************  Electrostatic Matrices  ************'
   opLabels(15) = '****************  Secular Equation  ***************'
   opLabels(16) = '************  Populate Electron Levels  ***********'
   opLabels(17) = '*************  Valence Charge Density  ************'
   opLabels(18) = '***************  Make SCF Potential  **************'
   opLabels(19) = '************  PDOS/TDOS/LI Calculation  ***********'
   opLabels(20) = '*************  PSCF Integrals and MME  ************'
   opLabels(21) = '*****  Symmetric Band Structure Path KPoints  *****'
   opLabels(22) = '********  Bond Order and Effective Charge  ********'
   opLabels(23) = '*************  Transition Calculation  ************'
   opLabels(24) = '*************  Parse The Command Line  ************'
   opLabels(25) = '************  Evaluate WaveFn on Mesh  ************'
   opLabels(26) = '************  Three Center Bond Order  ************'

end subroutine initOperationLabels

subroutine timeStampStart (opCode)

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer :: opCode

   ! Make sure we are given a valid opCode.
   if ((opCode > numOpCodes) .or. (opCode < 1)) then
      write (20,*) "opCode range is (1,",numOpCodes,")"
      write (20,*) "Given opCode is:  ",opCode
      stop
   endif

   ! Print the leading banner followed by the appropriate label.
   write (20,*)
   write (20,fmt='(a51)') '***************************************************'
   write (20,fmt='(a51)') opLabels(opCode)

   ! Log the date and time
   call date_and_time(tempDate,tempTime)
   date = tempDate(1:4)//'/'//tempDate(5:6)//'/'//tempDate(7:)
   time = tempTime(1:2)//':'//tempTime(3:4)//':'//tempTime(5:)
   write (20,100) date,time
   write (20,*)
   call flush (20)

   100 format ('Date is: ',a10,' Time is: ',a12) ! Date and time output format

end subroutine timeStampStart

subroutine timeStampEnd (opCode)

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer :: opCode

   ! Make sure we are given a valid opCode.
   if ((opCode > numOpCodes) .or. (opCode < 1)) then
      write (20,*) "opCode range is between 0 and ",numOpCodes
      write (20,*) "Given opCode is:  ",opCode
      stop
   endif

   ! Log the date and time.
   write (20,*)
   call date_and_time(tempDate,tempTime)
   date = tempDate(1:4)//'/'//tempDate(5:6)//'/'//tempDate(7:)
   time = tempTime(1:2)//':'//tempTime(3:4)//':'//tempTime(5:)
   write (20,100) date,time

   ! Print the appropriate label followed by the trailing banner.
   write (20,fmt='(a51)') opLabels(opCode)
   write (20,fmt='(a51)') '***************************************************'
   write (20,*)
   call flush (20)

   100 format ('Date is: ',a10,' Time is: ',a12) ! Date and time output format

end subroutine timeStampEnd

end module O_TimeStamps
