module O_CommandLine

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Integer to track which number command line argument is to be read in next.
   integer :: nextArg

   ! Variables used by all OLCAO programs.
   integer :: basisCode ! 1=MB; 2=FB; 3=EB

   ! Variables used by a subset of OLCAO programs.  (Individually labeled.)
   integer :: excitedQN_n  ! main, band, dos, bond, optc, wave
   integer :: excitedQN_l  ! main, band, dos, bond, optc, wave

   ! Variables used by main only.
   integer :: doDOS  ! Tack on a DOS calculation (1) or not (0).
   integer :: doBond ! Tack on a bond order calculation (1) or not (0).

   ! Variables used by pscf only.
   integer :: doMOME ! Include computation of the momentum matrix elements.
   integer :: doSYBD ! Do (1) a symmetric band structure calculation along a
                     !   path in the BZ, or do a normal band structure
                     !   calculation (0).

   ! Variables used by optc only.
   integer :: stateSet   ! Which set of states to use for transition
                         !   probability computation.  0=Occupied to
                         !   unoccupied (normal optical properties); 1=ground
                         !   state core level occupied to excited state
                         !   conduction band unoccupied (ELNES/XANES);
                         !   2=sigma(E).
   integer :: serialXYZ  ! Perform the XYZ components in serial (1) or not (0).

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseSetupCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseSetupCommandLine



subroutine parseMainCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   call readExcitedQN

   ! Read a flag indicating that a DOS calculation should be tacked on to the
   !   end of the SCF iterations.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) doDOS
   write (20,*) "doDOS = ",doDOS

   ! Read a flag indicating that a BOND calculation should be tacked on to the
   !   end of the SCF iterations.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) doBond
   write (20,*) "doBond = ",doBond

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseMainCommandLine


subroutine parsePSCFCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   call readExcitedQN

   ! Read a flag indicating that the momentum matrix elements should also be
   !   calculated.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) doMOME
   write (20,*) "doMOME = ",doMOME

   ! Read a flag indicating whether we should do a band structure calculation
   !   using the given kpoints or a symmetric band calculation on a path
   !   between high symmetry kpoints in the BZ.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) doSYBD
   write (20,*) "doSYBD = ",doSYBD

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parsePSCFCommandLine


subroutine parseDOSCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.
   call readBasisCode

   call readExcitedQN

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseDOSCommandLine



subroutine parseBondCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   call readExcitedQN

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseBondCommandLine



subroutine parseOptcCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   ! Read the type of transition calculation to do (which defines what set of
   !   states will be used).
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) stateSet
   write (20,*) "stateSet = ",stateSet

   call readExcitedQN

   ! Read a flag indicating whether or not the xyz components should be
   !   computed in serial (1) or not (0).
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) serialXYZ
   write (20,*) "serialXYZ = ",serialXYZ

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseOptcCommandLine


subroutine parseWaveCommandLine

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCode

   call readExcitedQN

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseWaveCommandLine


subroutine readBasisCode

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Get the command line argument that defines the basis set to use.
   !   1=MB; 2=FB; 3=EB
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) basisCode

   ! Check to make sure that the basis code is a valid number.
   if ((basisCode < 1) .or. (basisCode > 3)) then
      write (20,*) "CLP 'basisCode' = ",basisCode
      write (20,*) "It should be > 0 and < 4."
      write (20,*) "1=MB; 2=FB; 3=EB"
      stop
   else
      write (20,*) "basisCode = ",basisCode
      write (20,*) "1=MB; 2=FB; 3=EB"
      write (20,*)
   endif

end subroutine readBasisCode

subroutine readExcitedQN

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Store the command line argument that defines the QN_n of which electron
   !   will be excited.  (0=ground state; 1=K; 2=L; 3=M; 4=N; ...)
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) excitedQN_n

   ! Store the command line argument that defines the QN_l of which electron
   !   will be excited.  (0=s; 1=p; 2=d; 3=f; ...)
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) excitedQN_l

   ! Check to make sure that QN_l < QN_n (or that QN_l == QN_n == 0 (gs)).
   if ((excitedQN_l < excitedQN_n) .or. &
         & ((excitedQN_n == 0) .and. (excitedQN_l==0))) then
      write (20,*) "excitedQN_n = ",excitedQN_n
      write (20,*) "excitedQN_l = ",excitedQN_l
      write (20,*)
   else
      write (20,*) "CLP 'excitedQN_n' = ",excitedQN_n
      write (20,*) "CLP 'excitedQN_l' = ",excitedQN_l
      write (20,*) "QN_l should be < QN_n"
      stop
   endif

end subroutine readExcitedQN


subroutine initCLP

   ! Make sure that nothing funny is declared.
   implicit none

   ! Make sure that the first command line parameter to be read is #1.
   nextArg = 1

   ! Initialize all possible command line parameters to zero.
   basisCode   = 0
   excitedQN_n = 0
   excitedQN_l = 0
   doDOS       = 0
   doBond      = 0
   doMOME      = 0
   doSYBD      = 0
   stateSet    = 0
   serialXYZ   = 0

end subroutine initCLP


end module O_CommandLine
