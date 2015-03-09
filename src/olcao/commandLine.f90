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
   type commandLineParameters
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

      ! Variables used by intg only.
      integer :: doMOME ! Include computation of the momentum matrix elements.

      ! Variables used by band only.
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
   end type commandLineParameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseSetupCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseSetupCommandLine



subroutine parseMainCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   call readExcitedQN(clp)

   ! Read a flag indicating that a DOS calculation should be tacked on to the
   !   end of the SCF iterations.
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%doDOS
   write (20,*) "clp%doDOS = ",clp%doDOS

   ! Read a flag indicating that a BOND calculation should be tacked on to the
   !   end of the SCF iterations.
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%doBond
   write (20,*) "clp%doBond = ",clp%doBond

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseMainCommandLine


subroutine parseIntgCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   ! Read a flag indicating that the momentum matrix elements should also be
   !   calculated.
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%doMOME
   write (20,*) "clp%doMOME = ",clp%doMOME

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseIntgCommandLine


subroutine parseBandCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   ! Read a flag indicating whether we should do a band structure calculation
   !   using the given kpoints or a symmetric band calculation on a path
   !   between high symmetry kpoints in the BZ.
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%doSYBD
   write (20,*) "clp%doSYBD = ",clp%doSYBD

   call readExcitedQN(clp)

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseBandCommandLine



subroutine parseDOSCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   call readExcitedQN(clp)

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseDOSCommandLine



subroutine parseBondCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   call readExcitedQN(clp)

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseBondCommandLine



subroutine parseOptcCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   ! Read the type of transition calculation to do (which defines what set of
   !   states will be used).
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%stateSet
   write (20,*) "clp%stateSet = ",clp%stateSet

   call readExcitedQN(clp)

   ! Read a flag indicating whether or not the xyz components should be
   !   computed in serial (1) or not (0).
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%serialXYZ
   write (20,*) "clp%serialXYZ = ",clp%serialXYZ

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseOptcCommandLine


subroutine parseWaveCommandLine(clp)

   ! Use necessary modules.
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Open the file that will be written to as output for this program.
   open(20,file='fort.20',status='unknown',form='formatted')

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP(clp)

   ! Begin Parsing the command line.

   call readBasisCode(clp)

   call readExcitedQN(clp)

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseWaveCommandLine


subroutine readBasisCode(clp)

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Get the command line argument that defines the basis set to use.
   !   1=MB; 2=FB; 3=EB
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%basisCode

   ! Check to make sure that the basis code is a valid number.
   if ((clp%basisCode < 1) .or. (clp%basisCode > 3)) then
      write (20,*) "CLP 'clp%basisCode' = ",clp%basisCode
      write (20,*) "It should be > 0 and < 4."
      write (20,*) "1=MB; 2=FB; 3=EB"
      stop
   else
      write (20,*) "clp%basisCode = ",clp%basisCode
      write (20,*) "1=MB; 2=FB; 3=EB"
      write (20,*)
   endif

end subroutine readBasisCode

subroutine readExcitedQN(clp)

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Store the command line argument that defines the QN_n of which electron
   !   will be excited.  (0=ground state; 1=K; 2=L; 3=M; 4=N; ...)
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%excitedQN_n

   ! Store the command line argument that defines the QN_l of which electron
   !   will be excited.  (0=s; 1=p; 2=d; 3=f; ...)
   call getarg(clp%nextArg,commandBuffer)
   clp%nextArg = clp%nextArg + 1
   read (commandBuffer,*) clp%excitedQN_l

   ! Check to make sure that QN_l < QN_n (or that QN_l == QN_n == 0 (gs)).
   if ((clp%excitedQN_l < clp%excitedQN_n) .or. &
         & ((clp%excitedQN_n == 0) .and. (clp%excitedQN_l==0))) then
      write (20,*) "clp%excitedQN_n = ",clp%excitedQN_n
      write (20,*) "clp%excitedQN_l = ",clp%excitedQN_l
      write (20,*)
   else
      write (20,*) "CLP 'clp%excitedQN_n' = ",clp%excitedQN_n
      write (20,*) "CLP 'clp%excitedQN_l' = ",clp%excitedQN_l
      write (20,*) "QN_l should be < QN_n"
      stop
   endif

end subroutine readExcitedQN


subroutine initCLP(clp)

   ! Make sure that nothing funny is declared.
   implicit none

   ! Define the passed parameters.
   type(commandLineParameters), intent(inout) :: clp

   ! Make sure that the first command line parameter to be read is #1.
   clp%nextArg = 1

   ! Initialize all possible command line parameters to zero.
   clp%basisCode   = 0
   clp%excitedQN_n = 0
   clp%excitedQN_l = 0
   clp%doDOS       = 0
   clp%doBond      = 0
   clp%doMOME      = 0
   clp%doSYBD      = 0
   clp%stateSet    = 0
   clp%serialXYZ   = 0

end subroutine initCLP


end module O_CommandLine
