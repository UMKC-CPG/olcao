module O_Input

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Define a type that holds various inputs needed in a typical run.
   type inputData
      ! Define shared control variables that are used by all (or most) programs.
      integer :: numStates    ! The total number of states in the system.
      integer :: numElectrons ! The number of electrons in the system.

      ! Define control variables that apply only to the main program or its
      !   population subroutine (used by many other programs).
      integer :: iterFlagTDOS ! Compute TDOS with each SCF iteration.
      real (kind=double) :: thermalSigma ! Thermal smearing in eV stored in a.u.
         ! This is in units of eV which can be converted to K by multiplying by
         !    11604.505.  Some frequently used values are 0.001eV ~= 12K,
         !    0.005eV ~= 58K, 0.025eV ~= 290K, 0.05eV ~= 580K.
      real (kind=double) :: maxThermalRange ! This is the maximum range (given
         ! in eV and converted to au) that the thermal smearing will extend
         ! over. Think of this as a cutoff.
      real (kind=double) :: fermiSearchLimit ! See populate.f90 for docs. This
         ! helps set the boundaries of the search for the fermi level.

      ! Note that for all sigma values (e.g. sigmaDOS, sigmaBOND, etc.) the
      !   expectation is that the user provides the FWHM in eV of the Gaussian
      !   that they want used for broadening. We have the equation
      !   FWHM = 2*sqrt(2*ln2)*sigma so that from the given value we can easily
      !   compute the sigma to be used in the broadening subroutine. The sigma
      !   value is equal to the standard deviation of g(r) = 1/(sigma sqrt(2Pi))
      !   exp(-r^2 / 2 sigma^2).

      ! Define control variables that apply only to the DOS program.
      integer :: detailCodePDOS
      real (kind=double) :: deltaDOS  ! Given in eV, stored in a.u.
      real (kind=double) :: sigmaDOS  ! Given in eV, stored in a.u.
      real (kind=double) :: eminDOS   ! Given in eV, stored in a.u.
      real (kind=double) :: emaxDOS   ! Given in eV, stored in a.u.

      ! Define control variables that apply only to the bond program.
      integer :: detailCodeBond ! 0=old-style (by types) 1=new-style (all atom)
      integer :: doBond3C ! 0=no 1=yes
      integer :: maxNumNeighbors ! This is used to limit memory allocations. It
            ! is the max # of bonding neighbors that a single atom can have.
      real (kind=double) :: maxBL ! Given in Angstroms, stored in a.u.
      real (kind=double) :: deltaBOND  ! Given in eV, stored in a.u.
      real (kind=double) :: sigmaBOND  ! Given in eV, stored in a.u.
      real (kind=double) :: eminBOND   ! Given in eV, stored in a.u.
      real (kind=double) :: emaxBOND   ! Given in eV, stored in a.u.

      ! Define control variables that apply to the PACS use of the optc program.
      integer :: excitedAtomPACS ! Atom number of the excited atom.
      integer :: firstInitStatePACS
      integer :: lastInitStatePACS
      real (kind=double) :: deltaPACS  ! Given in eV, stored in a.u.
      real (kind=double) :: sigmaPACS  ! Given in eV, stored in a.u.
      real (kind=double) :: onsetEnergySlackPACS ! Given in eV, stored in a.u.
      real (kind=double) :: energyWindowPACS     ! Given in eV, stored in a.u.
      real (kind=double) :: totalEnergyDiffPACS  ! Given in eV, stored in a.u.

      ! Define control variables for the usual use of the optc program.
      real (kind=double) :: deltaOPTC      ! Given in eV, stored in a.u.
      real (kind=double) :: sigmaOPTC      ! Given in eV, stored in a.u.
      real (kind=double) :: maxTransEnOPTC ! Given in eV, stored in a.u.
      real (kind=double) :: cutoffEnOPTC   ! Given in eV, stored in a.u.

      ! Define control variables that apply to the SIGE use of the optc program.
      real (kind=double) :: deltaSIGE      ! Given in eV, stored in a.u.
      real (kind=double) :: sigmaSIGE      ! Given in eV, stored in a.u.
      real (kind=double) :: maxTransEnSIGE ! Given in eV, stored in a.u.
      real (kind=double) :: cutoffEnSIGE   ! Given in eV, stored in a.u.

      ! Define control variables that apply only to the wave program.
      integer :: doRho ! Flag indicating whether or not to include state
            ! occupancy when making the plot.
      integer :: styleWAVE ! Flag requesting the form of the output.  1=3D
            ! OpenDX & 1D projections; 2=3D OpenDX only; 3=1D projections only.
      real (kind=double) :: eminWAVE ! Given in eV, stored in a.u.
      real (kind=double) :: emaxWAVE ! Given in eV, stored in a.u.
   end type inputData

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseInput(inDat,clp)

   ! Import necessary modules.
   use O_Kinds
   use O_CommandLine
   use O_Lattice
   use O_AtomicTypes
   use O_AtomicSites
   use O_PotTypes
   use O_PotSites
   use O_ReadDataSubs
   use O_KPoints
   use O_ExchangeCorrelation
   use O_TimeStamps

   ! Make sure no funny variables are defined.
   implicit none

   ! define Passed parameters.
   type(commandLineParameters) :: clp ! from O_CommandLine
   type(inputData), intent(inout) :: inDat

   ! Define the variables that hold the file units for reading and writing.
   integer :: readUnit, writeUnit
   

   ! Log the beginning of the input reading.
   call timeStampStart(1)

   ! Set the unit that will be written to as output for this program.  (This
   !   has already been opened when the command line was read.)
   writeUnit = 20
   ! Open the primary input file for reading.
   open(5,file='fort.5',status='old',form='formatted')
   readUnit = 5

   ! Read the title.
   call readTitle(readUnit,writeUnit)

   ! Read the atomic type (basis) information.
   call readAtomicTypes(readUnit,writeUnit,clp)

   ! Read the potential type information.
   call readPotTypes(readUnit,writeUnit)

   ! Read the exchange correlation mesh parameters.
   call readExchCorrMeshParameters(readUnit,writeUnit)

   ! Read shared control parameters:  
   call readSharedControl(readUnit,writeUnit,inDat,clp)

   ! Read main input control parameters.
   call readMainControl(readUnit,writeUnit,inDat)

   ! Read pdos input control parameters.
   call readDOSControl(readUnit,writeUnit,inDat)

   ! Read bond input control parameters.
   call readBondControl(readUnit,writeUnit,inDat)

   ! Read sybd input control parameters.
   call readSYBDControl(readUnit,writeUnit)

   ! Read PACS input control parameters.
   call readPACSControl(readUnit,writeUnit,inDat,clp)

   ! Read optc input control parameters.
   call readOptcControl(readUnit,writeUnit,inDat)

   ! Read sige input control parameters.
   call readSigeControl(readUnit,writeUnit,inDat)

   ! Read wave input control parameters.
   call readWaveControl(readUnit,writeUnit,inDat)

   ! Close the primary input file.
   close(5)



   ! Open the structure input file for reading.
   open(4,file='fort.4',status='old',form='formatted')
   readUnit = 4


   ! Read the unit cell vectors.
   call readRealCellVectors(readUnit,writeUnit) 

   ! Read the atomic site information.
   call readAtomicSites(readUnit,writeUnit)

   ! Read the potential site information.
   call readPotSites(readUnit,writeUnit)

   ! Close the structure input file.
   close(4)



   ! Open the kpoint mesh input file for reading.
   open(15,file='fort.15',status='old',form='formatted')
   readUnit = 15

   ! Read the kpoint mesh information.
   call readKPoints(readUnit,writeUnit)

   ! Close the kpoint mesh input file.
   close(15)

   ! Log the ending of the input reading.
   call timeStampEnd(1)

end subroutine parseInput


subroutine readTitle(readUnit,writeUnit)

   ! Import necessary modules.
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   character*80 :: titleLine

   ! Read the header for this input section.
   call readAndCheckLabel(readUnit,writeUnit,len('TITLE'),'TITLE')

   ! Read lines and regurgitate them until we find the END_TITLE tag.
   do while (.true.)
      read(5,*) titleLine
      write (20,fmt='(a80)') titleLine

      if (titleLine == 'END_TITLE') then
         exit
      endif
   enddo

end subroutine readTitle



subroutine readSharedControl(readUnit,writeUnit,inDat,clp)

   ! Import necessary modules.
   use O_Kinds
   use O_CommandLine
   use O_Lattice
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat
   type(commandLineParameters), intent(inout) :: clp

   ! Define local variables.
   real (kind=double) :: basisFnConvgTemp, electroConvgTemp ! Input parameters
         ! that define the thresholds used for basis function and electrostatic
         ! potential convergence.  (Cutoff distances for negligability). 
   real (kind=double), dimension(3) :: tempArray


   ! Read the header for this input section.
   call readAndCheckLabel(readUnit,writeUnit,len('SHARED_INPUT_DATA'),'SHARED_INPUT_DATA')


   ! Read the cutoff values for integrals.  (Used by setup and intg.)
   call readData(readUnit,writeUnit,basisFnConvgTemp,electroConvgTemp,&
         & len('BASISFUNCTION_AND_ELECTROSTATIC_CUTOFFS'),&
         & 'BASISFUNCTION_AND_ELECTROSTATIC_CUTOFFS')

   ! Check that the given values are valid.
   if (basisFnConvgTemp <= 0.0_double .or. basisFnConvgTemp > 1.0_double) then
      write (20,*) 'Wave Fn Convergence value is not in (0,1]'
      write (20,*) 'This is an error'
      stop
   endif
   if (electroConvgTemp <= 0.0_double .or. electroConvgTemp > 1.0_double) then
      write (20,*) 'Electrostatic Convergence value is not in (0,1]'
      write (20,*) 'This is an error'
      stop
   endif

   ! Immediately calculate and store the negative log of the above convergence
   !   values.  Please read the comments associated with this subroutine to
   !   understand the meaning of these values.
   call setCutoffThresh(basisFnConvgTemp,electroConvgTemp)

   write (20,FMT="(a,e12.5)")'Maximum basisFunction threshold = ',&
         & logBasisFnThresh
   write (20,FMT="(a,e12.5)")'Maximum electroStatic threshold = ',&
         & logElecThresh



   ! Read the number of states to use for the command line given basis code.
   call readData(readUnit,writeUnit,3,tempArray(:),len('NUM_STATES_TO_USE'),'NUM_STATES_TO_USE')
   inDat%numStates = tempArray(clp%basisCode)
   write (20,*) "Using ",inDat%numStates," states."
   call flush (20)


   ! Read the number of valence electrons in the system.
   call readData(readUnit,writeUnit,inDat%numElectrons,len('NUM_ELECTRONS'),'NUM_ELECTRONS')


   ! Read the thermal smearing parameters.
   call readData(readUnit,writeUnit,inDat%thermalSigma,inDat%maxThermalRange,&
         & len('THERMAL_SMEARING_SIGMA__SIGMA_CUTOFF'),&
         & 'THERMAL_SMEARING_SIGMA__SIGMA_CUTOFF')

   ! Apply convertion to atomic units.
   inDat%thermalSigma = inDat%thermalSigma / hartree

   ! Read the value that helps set the upper and lower limit on the search for
   !   the fermi level during the thermal smearing operation.
   call readData(readUnit,writeUnit,inDat%fermiSearchLimit,&
         & len('FERMI_LEVEL_SEARCH_LIMIT'),&
         & 'FERMI_LEVEL_SEARCH_LIMIT')

   ! Apply convertion to atomic units.
   inDat%fermiSearchLimit = inDat%fermiSearchLimit / hartree



end subroutine readSharedControl


subroutine readMainControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_Kinds
   use O_ReadDataSubs
   use O_PotTypes
   use O_Potential

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   ! Define local variables.
   integer :: i
   integer :: lastIt
   integer :: fbLevel
   integer :: XC_CODETemp
   integer :: numSpltTps
   integer :: typeID
   real (kind=double) :: rlxFact
   real (kind=double) :: cTest
   real (kind=double) :: spltAmt
   real (kind=double), allocatable, dimension (:) :: typeSpinSplitTemp


   call readAndCheckLabel(readUnit,writeUnit,len('MAIN_INPUT_DATA'),'MAIN_INPUT_DATA')

   call readData(readUnit,writeUnit,lastIt,len('LAST_ITERATION'),'LAST_ITERATION')

   call readData(readUnit,writeUnit,cTest,len('CONVERGENCE_TEST'),'CONVERGENCE_TEST')

   call readData(readUnit,writeUnit,XC_CODETemp,len('XC_CODE'),'XC_CODE')

   call readData(readUnit,writeUnit,fbLevel,len('FEEDBACK_LEVEL'),'FEEDBACK_LEVEL')

   call readData(readUnit,writeUnit,rlxFact,len('RELAXATION_FACTOR'),'RELAXATION_FACTOR')

   call readData(readUnit,writeUnit,inDat%iterFlagTDOS,len('EACH_ITER_FLAGS__TDOS'),&
         & 'EACH_ITER_FLAGS__TDOS')

   ! Read the default amount that all types will have their charge contribution
   !   split by in the case of a spin polarized calculation.
   call readData(readUnit,writeUnit,numSpltTps,spltAmt,len('NUM_SPLIT_TYPES__DEFAULT_SPLIT'),&
         & 'NUM_SPLIT_TYPES__DEFAULT_SPLIT')

   ! Assign the default spin splitting for all types in the system.
   allocate (typeSpinSplitTemp(numPotTypes))
   do i = 1, numPotTypes
      typeSpinSplitTemp(i) = spltAmt
   enddo

   ! Read specific splitting for specific types.
   call readAndCheckLabel(readUnit,writeUnit,len('TYPE_ID__SPIN_SPLIT_FACTOR'),&
         &'TYPE_ID__SPIN_SPLIT_FACTOR')
   do i = 1, numSpltTps
      call readData(readUnit,writeUnit,typeID,spltAmt,0,'')
      typeSpinSplitTemp(typeID) = spltAmt
   enddo

   call setPotControlParameters(fbLevel,lastIt,XC_CODETemp,rlxFact,cTest,&
         & typeSpinSplitTemp)

   deallocate (typeSpinSplitTemp)

end subroutine readMainControl


subroutine readDOSControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   call readAndCheckLabel(readUnit,writeUnit,len('DOS_INPUT_DATA'),&
                              'DOS_INPUT_DATA')

   ! Read all the input data.
   call readData(readUnit,writeUnit,inDat%deltaDOS,inDat%sigmaDOS,0,'')
   call readData(readUnit,writeUnit,inDat%eminDOS,inDat%emaxDOS,0,'')
   call readData(readUnit,writeUnit,inDat%detailCodePDOS,0,'')

   ! Apply necessary conversions to the units.  (eV -> hartree atomic units)
   inDat%deltaDOS = inDat%deltaDOS / hartree
   inDat%sigmaDOS = inDat%sigmaDOS / hartree / (2.0_double * sqrt(2.0_double *&
         & log(2.0_double)))
   inDat%eminDOS  = inDat%eminDOS  / hartree
   inDat%emaxDOS  = inDat%emaxDOS  / hartree

end subroutine readDOSControl


subroutine readBondControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   call readAndCheckLabel(readUnit,writeUnit,len('BOND_INPUT_DATA'),'BOND_INPUT_DATA')

   ! Read the data.
   call readData(readUnit,writeUnit,inDat%maxBL,0,'')
   call readData(readUnit,writeUnit,inDat%deltaBOND,inDat%sigmaBOND,0,'')
   call readData(readUnit,writeUnit,inDat%eminBOND,inDat%emaxBOND,0,'')
   call readData(readUnit,writeUnit,inDat%detailCodeBond,0,'')
   call readData(readUnit,writeUnit,inDat%doBond3C,0,'')
   call readData(readUnit,writeUnit,inDat%maxNumNeighbors,0,'')

   ! Apply necessary conversions to the units.  (Angstroms to bohrRadii a.u.
   !   and eV to hartree atomic units.)
   inDat%maxBL     = inDat%maxBL     / bohrRad
   inDat%deltaBOND = inDat%deltaBOND / hartree
   inDat%sigmaBOND = inDat%sigmaBOND / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))
   inDat%eminBOND  = inDat%eminBOND  / hartree
   inDat%emaxBOND  = inDat%emaxBOND  / hartree

end subroutine readBondControl


subroutine readSYBDControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_KPoints

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readSYBDKPoints(readUnit, writeUnit)

end subroutine readSYBDControl


subroutine readPACSControl(readUnit,writeUnit,inDat,clp)

   ! Use necessary modules
   use O_Kinds
   use O_ReadDataSubs
   use O_CommandLine
   use O_Constants

   ! Make sure no unwanted variables are declared.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat
   type(commandLineParameters), intent(inout) :: clp

   ! Define local variables.
   integer :: i
   integer :: numCorePACS
   real (kind=double), dimension(5) :: tempArray
   integer :: QN_n
   integer :: QN_l

   ! Initialize the first and last initial states.
   inDat%firstInitStatePACS = 0
   inDat%lastInitStatePACS  = 0


   ! Read the input control data for PACS calculations.
   call readAndCheckLabel(readUnit,writeUnit,len('PACS_INPUT_DATA'),'PACS_INPUT_DATA')

   call readData(readUnit,writeUnit,inDat%excitedAtomPACS,0,'')
   call readData(readUnit,writeUnit,inDat%deltaPACS,inDat%sigmaPACS,0,'')
   call readData(readUnit,writeUnit,inDat%onsetEnergySlackPACS,&
                  inDat%energyWindowPACS,0,'')
   call readData(readUnit,writeUnit,numCorePACS,0,'')
   do i = 1, numCorePACS
      call readData(readUnit,writeUnit,5,tempArray,0,'')
      QN_n = int(tempArray(1))
      QN_l = int(tempArray(2))
      if ((QN_n == clp%excitedQN_n) .and. (QN_l == clp%excitedQN_l)) then
         inDat%firstInitStatePACS = int(tempArray(3))
         inDat%lastInitStatePACS  = int(tempArray(4))
         inDat%totalEnergyDiffPACS = tempArray(5)
      endif
   enddo

   ! Apply necessary conversions of units.
   inDat%deltaPACS = inDat%deltaPACS / hartree
   inDat%sigmaPACS = inDat%sigmaPACS / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))
   inDat%onsetEnergySlackPACS = inDat%onsetEnergySlackPACS / hartree
   inDat%energyWindowPACS     = inDat%energyWindowPACS / hartree
   inDat%totalEnergyDiffPACS  = inDat%totalEnergyDiffPACS / hartree


end subroutine readPACSControl


subroutine readOptcControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   ! Read input OPTC control data.
   call readAndCheckLabel(readUnit,writeUnit,len('OPTC_INPUT_DATA'),'OPTC_INPUT_DATA')
   call readData(readUnit,writeUnit,inDat%cutoffEnOPTC,0,'')
   call readData(readUnit,writeUnit,inDat%maxTransEnOPTC,0,'')
   call readData(readUnit,writeUnit,inDat%deltaOPTC,0,'')
   call readData(readUnit,writeUnit,inDat%sigmaOPTC,0,'')

   ! Apply necessary conversions to a.u.
   inDat%cutoffEnOPTC   = inDat%cutoffEnOPTC / hartree
   inDat%maxTransEnOPTC = inDat%maxTransEnOPTC / hartree
   inDat%deltaOPTC = inDat%deltaOPTC / hartree
   inDat%sigmaOPTC = inDat%sigmaOPTC / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))

end subroutine readOptcControl


subroutine readSigeControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   ! Read the input sigma(E) control variables.
   call readAndCheckLabel(readUnit,writeUnit,len('SIGE_INPUT_DATA'),'SIGE_INPUT_DATA')
   call readData(readUnit,writeUnit,inDat%cutoffEnSIGE,0,'')
   call readData(readUnit,writeUnit,inDat%maxTransEnSIGE,0,'')
   call readData(readUnit,writeUnit,inDat%deltaSIGE,0,'')
   call readData(readUnit,writeUnit,inDat%sigmaSIGE,0,'')

   ! Apply necessary conversions to a.u.
   inDat%cutoffEnSIGE   = inDat%cutoffEnSIGE / hartree
   inDat%maxTransEnSIGE = inDat%maxTransEnSIGE / hartree
   inDat%deltaSIGE = inDat%deltaSIGE / hartree
   inDat%sigmaSIGE = inDat%sigmaSIGE / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))

end subroutine readSigeControl


subroutine readWaveControl(readUnit,writeUnit,inDat)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Lattice
   use O_Constants

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.
   type(inputData), intent(inout) :: inDat

   ! Read the input wave function / charge density plotting control variables.
   call readAndCheckLabel(readUnit,writeUnit,len('WAVE_INPUT_DATA'),'WAVE_INPUT_DATA')
   call readNumMeshPoints(readUnit,writeUnit)
   call readData(readUnit,writeUnit,inDat%eminWAVE,inDat%emaxWAVE,0,'')
   call readData(readUnit,writeUnit,inDat%doRho,0,'')
   call readData(readUnit,writeUnit,inDat%styleWAVE,0,'')

   ! Apply the necessary conversions of the data to atomic units.
   inDat%eminWAVE = inDat%eminWAVE / hartree
   inDat%emaxWAVE = inDat%emaxWAVE / hartree

end subroutine readWaveControl


end module O_Input
