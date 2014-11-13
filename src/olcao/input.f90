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

   ! Define control variables that apply only to the DOS program.
   integer :: detailCodePDOS
   real (kind=double) :: deltaDOS  ! Given in eV, stored in a.u.
   real (kind=double) :: sigmaDOS  ! Given in eV, stored in a.u.
   real (kind=double) :: eminDOS   ! Given in eV, stored in a.u.
   real (kind=double) :: emaxDOS   ! Given in eV, stored in a.u.

   ! Define control variables that apply only to the bond program.
   integer :: detailCodeBond
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

   ! Define control variables that apply to the usual use of the optc program.
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
   integer :: doRho ! Flag indicating whether or not to include state occupancy
         ! when making the plot.
   integer :: styleWAVE ! Flag requesting the form of the output.  1=3D OpenDX
         ! and 1D projections; 2=3D OpenDX only; 3=1D projections only.
   real (kind=double) :: eminWAVE ! Given in eV, stored in a.u.
   real (kind=double) :: emaxWAVE ! Given in eV, stored in a.u.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseInput

   ! Import necessary modules.
   use O_Kinds
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

   ! Log the beginning of the input reading.
   call timeStampStart(1)

   ! Set the unit that will be written to as output for this program.  (This
   !   has already been opened when the command line was read.)
   call setWriteUnit(20)

   ! Open the primary input file for reading.
   open(5,file='fort.5',status='old',form='formatted')
   call setReadUnit(5)

   ! Read the title.
   call readTitle

   ! Read the atomic type (basis) information.
   call readAtomicTypes

   ! Read the potential type information.
   call readPotTypes

   ! Read the exchange correlation mesh parameters.
   call readExchCorrMeshParameters

   ! Read shared control parameters:  
   call readSharedControl

   ! Read main input control parameters.
   call readMainControl

   ! Read pdos input control parameters.
   call readDOSControl

   ! Read bond input control parameters.
   call readBondControl

   ! Read sybd input control parameters.
   call readSYBDControl

   ! Read PACS input control parameters.
   call readPACSControl

   ! Read optc input control parameters.
   call readOptcControl

   ! Read sige input control parameters.
   call readSigeControl

   ! Read wave input control parameters.
   call readWaveControl

   ! Close the primary input file.
   close(5)



   ! Open the structure input file for reading.
   open(4,file='fort.4',status='old',form='formatted')
   call setReadUnit(4)

   ! Read the unit cell vectors.
   call readRealCellVectors

   ! Read the atomic site information.
   call readAtomicSites

   ! Read the potential site information.
   call readPotSites

   ! Close the structure input file.
   close(4)



   ! Open the kpoint mesh input file for reading.
   open(15,file='fort.15',status='old',form='formatted')
   call setReadUnit(15)

   ! Read the kpoint mesh information.
   call readKPoints

   ! Close the kpoint mesh input file.
   close(15)

   ! Log the ending of the input reading.
   call timeStampEnd(1)

end subroutine parseInput


subroutine readTitle

   ! Import necessary modules.
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   character*80 :: titleLine

   ! Read the header for this input section.
   call readAndCheckLabel(len('TITLE'),'TITLE')

   ! Read lines and regurgitate them until we find the END_TITLE tag.
   do while (.true.)
      read(5,*) titleLine
      write (20,fmt='(a80)') titleLine

      if (titleLine == 'END_TITLE') then
         exit
      endif
   enddo

end subroutine readTitle



subroutine readSharedControl

   ! Import necessary modules.
   use O_Kinds
   use O_CommandLine
   use O_Lattice
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   real (kind=double) :: basisFnConvgTemp, electroConvgTemp ! Input parameters
         ! that define the thresholds used for basis function and electrostatic
         ! potential convergence.  (Cutoff distances for negligability). 
   real (kind=double), dimension(3) :: tempArray


   ! Read the header for this input section.
   call readAndCheckLabel(len('SHARED_INPUT_DATA'),'SHARED_INPUT_DATA')


   ! Read the cutoff values for integrals.  (Used by setup and intg.)
   call readData(basisFnConvgTemp,electroConvgTemp,&
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
   call readData(3,tempArray(:),len('NUM_STATES_TO_USE'),'NUM_STATES_TO_USE')
   numStates = tempArray(basisCode)
write (20,*) "Using ",numStates," states."
call flush (20)


   ! Read the number of valence electrons in the system.
   call readData(numElectrons,len('NUM_ELECTRONS'),'NUM_ELECTRONS')


   ! Read the thermal smearing parameter.
   call readData(thermalSigma,len('THERMAL_SMEARING_SIGMA'),&
         & 'THERMAL_SMEARING_SIGMA')

   ! Apply convertions in to atomic units.
   thermalSigma = thermalSigma / hartree


end subroutine readSharedControl


subroutine readMainControl

   ! Use necessary modules
   use O_Kinds
   use O_ReadDataSubs
   use O_PotTypes
   use O_Potential

   implicit none

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


   call readAndCheckLabel(len('MAIN_INPUT_DATA'),'MAIN_INPUT_DATA')

   call readData(lastIt,len('LAST_ITERATION'),'LAST_ITERATION')

   call readData(cTest,len('CONVERGENCE_TEST'),'CONVERGENCE_TEST')

   call readData(XC_CODETemp,len('XC_CODE'),'XC_CODE')

   call readData(fbLevel,len('FEEDBACK_LEVEL'),'FEEDBACK_LEVEL')

   call readData(rlxFact,len('RELAXATION_FACTOR'),'RELAXATION_FACTOR')

   call readData(iterFlagTDOS,len('EACH_ITER_FLAGS__TDOS'),&
         & 'EACH_ITER_FLAGS__TDOS')

   ! Read the default amount that all types will have their charge contribution
   !   split by in the case of a spin polarized calculation.
   call readData(numSpltTps,spltAmt,len('NUM_SPLIT_TYPES__DEFAULT_SPLIT'),&
         & 'NUM_SPLIT_TYPES__DEFAULT_SPLIT')

   ! Assign the default spin splitting for all types in the system.
   allocate (typeSpinSplitTemp(numPotTypes))
   do i = 1, numPotTypes
      typeSpinSplitTemp(i) = spltAmt
   enddo

   ! Read specific splitting for specific types.
   call readAndCheckLabel(len('TYPE_ID__SPIN_SPLIT_FACTOR'),&
         &'TYPE_ID__SPIN_SPLIT_FACTOR')
   do i = 1, numSpltTps
      call readData (typeID,spltAmt,0,'')
      typeSpinSplitTemp(typeID) = spltAmt
   enddo

   call setPotControlParameters(fbLevel,lastIt,XC_CODETemp,rlxFact,cTest,&
         & typeSpinSplitTemp)

   deallocate (typeSpinSplitTemp)

end subroutine readMainControl


subroutine readDOSControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   call readAndCheckLabel(len('DOS_INPUT_DATA'),'DOS_INPUT_DATA')

   ! Read all the input data.
   call readData(deltaDOS,sigmaDOS,0,'')
   call readData(eminDOS,emaxDOS,0,'')
   call readData(detailCodePDOS,0,'')

   ! Apply necessary conversions to the units.  (eV -> hartree atomic units)
   deltaDOS = deltaDOS / hartree
   sigmaDOS = sigmaDOS / hartree
   eminDOS  = eminDOS  / hartree
   emaxDOS  = emaxDOS  / hartree

end subroutine readDOSControl


subroutine readBondControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   call readAndCheckLabel(len('BOND_INPUT_DATA'),'BOND_INPUT_DATA')

   ! Read the data.
   call readData(maxBL,0,'')
   call readData(deltaBOND,sigmaBOND,0,'')
   call readData(eminBOND,emaxBOND,0,'')
   call readData(detailCodeBond,0,'')

   ! Apply necessary conversions to the units.  (Angstroms to bohrRadii a.u.
   !   and eV to hartree atomic units.)
   maxBL     = maxBL     / bohrRad
   deltaBOND = deltaBOND / hartree
   sigmaBOND = sigmaBOND / hartree
   eminBOND  = eminBOND  / hartree
   emaxBOND  = emaxBOND  / hartree

end subroutine readBondControl


subroutine readSYBDControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_KPoints

   implicit none

   call readSYBDKPoints

end subroutine readSYBDControl


subroutine readPACSControl

   ! Use necessary modules
   use O_Kinds
   use O_ReadDataSubs
   use O_CommandLine
   use O_Constants

   ! Make sure no unwanted variables are declared.
   implicit none

   ! Define local variables.
   integer :: i
   integer :: numCorePACS
   real (kind=double), dimension(5) :: tempArray
   integer :: QN_n
   integer :: QN_l

   ! Initialize the first and last initial states.
   firstInitStatePACS = 0
   lastInitStatePACS  = 0


   ! Read the input control data for PACS calculations.
   call readAndCheckLabel(len('PACS_INPUT_DATA'),'PACS_INPUT_DATA')

   call readData(excitedAtomPACS,0,'')
   call readData(deltaPACS,sigmaPACS,0,'')
   call readData(onsetEnergySlackPACS,energyWindowPACS,0,'')
   call readData(numCorePACS,0,'')
   do i = 1, numCorePACS
      call readData(5,tempArray,0,'')
      QN_n = int(tempArray(1))
      QN_l = int(tempArray(2))
      if ((QN_n == excitedQN_n) .and. (QN_l == excitedQN_l)) then
         firstInitStatePACS = int(tempArray(3))
         lastInitStatePACS  = int(tempArray(4))
         totalEnergyDiffPACS = tempArray(5)
      endif
   enddo

   ! Apply necessary conversions of units.
   deltaPACS = deltaPACS / hartree
   sigmaPACS = sigmaPACS / hartree
   onsetEnergySlackPACS = onsetEnergySlackPACS / hartree
   energyWindowPACS     = energyWindowPACS / hartree
   totalEnergyDiffPACS  = totalEnergyDiffPACS / hartree


end subroutine readPACSControl


subroutine readOptcControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! Read input OPTC control data.
   call readAndCheckLabel(len('OPTC_INPUT_DATA'),'OPTC_INPUT_DATA')
   call readData(cutoffEnOPTC,0,'')
   call readData(maxTransEnOPTC,0,'')
   call readData(deltaOPTC,0,'')
   call readData(sigmaOPTC,0,'')

   ! Apply necessary conversions to a.u.
   cutoffEnOPTC   = cutoffEnOPTC / hartree
   maxTransEnOPTC = maxTransEnOPTC / hartree
   deltaOPTC = deltaOPTC / hartree
   sigmaOPTC = sigmaOPTC / hartree

end subroutine readOptcControl


subroutine readSigeControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Constants

   implicit none

   ! Read the input sigma(E) control variables.
   call readAndCheckLabel(len('SIGE_INPUT_DATA'),'SIGE_INPUT_DATA')
   call readData(cutoffEnSIGE,0,'')
   call readData(maxTransEnSIGE,0,'')
   call readData(deltaSIGE,0,'')
   call readData(sigmaSIGE,0,'')

   ! Apply necessary conversions to a.u.
   cutoffEnSIGE   = cutoffEnSIGE / hartree
   maxTransEnSIGE = maxTransEnSIGE / hartree
   deltaSIGE = deltaSIGE / hartree
   sigmaSIGE = sigmaSIGE / hartree

end subroutine readSigeControl


subroutine readWaveControl

   ! Use necessary modules
   use O_ReadDataSubs
   use O_Lattice
   use O_Constants

   implicit none

   ! Read the input wave function / charge density plotting control variables.
   call readAndCheckLabel(len('WAVE_INPUT_DATA'),'WAVE_INPUT_DATA')
   call readNumMeshPoints
   call readData(eminWAVE,emaxWAVE,0,'')
   call readData(doRho,0,'')
   call readData(styleWAVE,0,'')

   ! Apply the necessary conversions of the data to atomic units.
   eminWAVE = eminWAVE / hartree
   emaxWAVE = emaxWAVE / hartree

end subroutine readWaveControl


end module O_Input
