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
   integer :: outputCodeBondQ ! 0=new (all atom Machine readable) 1=old (types)
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseInput

   ! Import necessary modules.
   use O_Kinds
   use O_Lattice, only: readRealCellVectors
   use O_KPoints, only: readKPoints
   use O_PotTypes, only: readPotTypes
   use O_PotSites, only: readPotSites
   use O_AtomicTypes, only: readAtomicTypes
   use O_AtomicSites, only: readAtomicSites
   use O_ExchangeCorrelation, only: readExchCorrMeshParameters
   use O_ReadDataSubs
   use O_TimeStamps

   ! Make sure no funny variables are defined.
   implicit none

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
   call readAtomicTypes(readUnit,writeUnit)

   ! Read the potential type information.
   call readPotTypes(readUnit,writeUnit)

   ! Read the exchange correlation mesh parameters.
   call readExchCorrMeshParameters(readUnit,writeUnit)

   ! Read shared control parameters:  
   call readSharedControl(readUnit,writeUnit)

   ! Read main input control parameters.
   call readMainControl(readUnit,writeUnit)

   ! Read pdos input control parameters.
   call readDOSControl(readUnit,writeUnit)

   ! Read bond input control parameters.
   call readBondControl(readUnit,writeUnit)

   ! Read sybd input control parameters.
   call readSYBDControl(readUnit,writeUnit)

   ! Read PACS input control parameters.
   call readPACSControl(readUnit,writeUnit)

   ! Read optc input control parameters.
   call readOptcControl(readUnit,writeUnit)

   ! Read sige input control parameters.
   call readSigeControl(readUnit,writeUnit)

   ! Read wave input control parameters.
   call readWaveControl(readUnit,writeUnit)

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



subroutine readSharedControl(readUnit,writeUnit)

   ! Import necessary modules.
   use O_Kinds
   use O_Constants, only: hartree
   use O_CommandLine, only: basisCode
   use O_Lattice, only: logBasisFnThresh, logElecThresh, setCutoffThresh
   use O_ReadDataSubs

   ! Make sure no funny variables are defined.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   real (kind=double) :: basisFnConvgTemp, electroConvgTemp ! Input parameters
         ! that define the thresholds used for basis function and electrostatic
         ! potential convergence.  (Cutoff distances for negligability). 
   real (kind=double), dimension(3) :: tempArray


   ! Read the header for this input section.
   call readAndCheckLabel(readUnit,writeUnit,len('SHARED_INPUT_DATA'),&
         & 'SHARED_INPUT_DATA')


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
   call readData(readUnit,writeUnit,3,tempArray(:),len('NUM_STATES_TO_USE'),&
         & 'NUM_STATES_TO_USE')
   numStates = tempArray(basisCode)
   write (20,*) "Using ",numStates," states."
   call flush (20)


   ! Read the number of valence electrons in the system.
   call readData(readUnit,writeUnit,numElectrons,len('NUM_ELECTRONS'),&
      & 'NUM_ELECTRONS')


   ! Read the thermal smearing parameters.
   call readData(readUnit,writeUnit,thermalSigma,maxThermalRange,&
         & len('THERMAL_SMEARING_SIGMA__SIGMA_CUTOFF'),&
         & 'THERMAL_SMEARING_SIGMA__SIGMA_CUTOFF')

   ! Apply convertion to atomic units.
   thermalSigma = thermalSigma / hartree

   ! Read the value that helps set the upper and lower limit on the search for
   !   the fermi level during the thermal smearing operation.
   call readData(readUnit,writeUnit,fermiSearchLimit,&
         & len('FERMI_LEVEL_SEARCH_LIMIT'),&
         & 'FERMI_LEVEL_SEARCH_LIMIT')

   ! Apply convertion to atomic units.
   fermiSearchLimit = fermiSearchLimit / hartree



end subroutine readSharedControl


subroutine readMainControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Kinds
   use O_ReadDataSubs
   use O_Constants, only: hartree
   use O_PotTypes, only: numPotTypes
   use O_Potential, only: setPotControlParameters

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define local variables.
   integer :: i
   integer :: lastIt
   integer :: fbLevel
   integer :: xcCodeTemp
   integer :: plusUJFormTemp
   integer :: numPlusUJItemsTemp
   integer :: numSplitTypes
   integer :: itemID
   integer :: typeID
   integer, allocatable, dimension (:) :: plusUJItemIDTemp
   real (kind=double) :: rlxFact
   real (kind=double) :: cTest
   real (kind=double) :: hubbardU
   real (kind=double) :: hundJ
   real (kind=double) :: splitAmt
   real (kind=double), allocatable, dimension (:,:) :: plusUJItemValueTemp
   real (kind=double), allocatable, dimension (:) :: typeSpinSplitTemp


   call readAndCheckLabel(readUnit,writeUnit,len('MAIN_INPUT_DATA'),&
         & 'MAIN_INPUT_DATA')

   call readData(readUnit,writeUnit,lastIt,len('LAST_ITERATION'),&
         & 'LAST_ITERATION')

   call readData(readUnit,writeUnit,cTest,len('CONVERGENCE_TEST'),&
         & 'CONVERGENCE_TEST')

   call readData(readUnit,writeUnit,xcCodeTemp,len('XC_CODE'),'XC_CODE')

   call readData(readUnit,writeUnit,fbLevel,len('FEEDBACK_LEVEL'),&
         & 'FEEDBACK_LEVEL')

   call readData(readUnit,writeUnit,rlxFact,len('RELAXATION_FACTOR'),&
         & 'RELAXATION_FACTOR')

   call readData(readUnit,writeUnit,iterFlagTDOS,len('EACH_ITER_FLAGS__TDOS'),&
         & 'EACH_ITER_FLAGS__TDOS')

   ! At this point we need to read information related to the plusUJ terms.
   !   We will accept the data in three possible forms: by element, sequential
   !   type number, or atom. The determination is made on the first line of
   !   this section: (0, 1, 2) respectively. Depending on the value of the
   !   flag, the second number indicates the number of elements, types, or
   !   atoms for which UJ data is provided. The remaining data is therefore
   !   understood to be useful for individual atoms or for types or for
   !   elements. All of the data will be converted to the form of individual
   !   atoms internally within OLCAO, but in many cases it is assumed to be
   !   easier to specify the UJ terms for all atoms of a given element or
   !   sequential type number all at once instead of on a per-atom basis.

   ! Read the form in which the UJ terms will be given and the number of items.
   !   Presently, it is assumed that the UJ contribution will only be added to
   !   the orbital type (d or f) that is highest in energy, meaning that an
   !   atom with 4f orbitals cannot have a UJ contribution to the 3d orbitals.
   !   Note that the term "Item" is used here because depending on the Form,
   !   the input may be providing UJ values for all atoms of a given element,
   !   type, or even just for specific atoms.
   call readData(readUnit,writeUnit,plusUJFormTemp,numPlusUJItemsTemp,&
         & len('PLUSUJ_FORM__NUM_ITEMS'),'PLUSUJ_FORM__NUM_ITEMS')

   ! Allocate space to hold the U and J values for the designated number of
   !   atoms. Also allocate space to hold the type identifies.
   allocate (plusUJItemValueTemp(2,numPlusUJItemsTemp))
   allocate (plusUJItemIDTemp(numPlusUJItemsTemp))

   ! Read the U and J values for each of the specific ID'd items. For the
   !   element or sequential type form, the ID number given will correspond to
   !   the ID number given under the ATOM_TYPE_ID__SEQUENTIAL_NUMBER tag in the
   !   olcao.dat file (the first or last number of that line). For atoms, it
   !   corresponds to the atom number. Note that the U and J quantities must be
   !   provided in eV and that they will be converted to au.
   ! For reference, the ATOM_TYPE_ID__SEQUENTIAL_NUMBER tag provides a series
   !   of ID numbers. The first is the element ID, then the species ID within
   !   that element, and then the type ID within that species. The final number
   !   is the overall (system-wide) type ID.
   call readAndCheckLabel(readUnit,writeUnit,len('ID__U__J'),'ID__U__J')
   do i = 1, numPlusUJItemsTemp
      call readData(readUnit,writeUnit,itemID,hubbardU,hundJ,0,'')
      plusUJItemValueTemp(1,i) = hubbardU / hartree ! Convert from eV to au
      plusUJItemValueTemp(2,i) = hundJ / hartree ! Convert from eV to au
      plusUJItemIDTemp(i) = itemID
   enddo


   ! Read the default amount that all types will have their charge contribution
   !   split by in the case of a spin polarized calculation.
   call readData(readUnit,writeUnit,numSplitTypes,splitAmt,&
         & len('NUM_SPLIT_TYPES__DEFAULT_SPLIT'),&
         & 'NUM_SPLIT_TYPES__DEFAULT_SPLIT')

   ! Assign the default spin splitting for all types in the system.
   allocate (typeSpinSplitTemp(numPotTypes))
   typeSpinSplitTemp(:) = splitAmt

   ! Read specific splitting for specific types.
   call readAndCheckLabel(readUnit,writeUnit,&
         & len('TYPE_ID__SPIN_SPLIT_FACTOR'),'TYPE_ID__SPIN_SPLIT_FACTOR')
   do i = 1, numSplitTypes
      call readData(readUnit,writeUnit,typeID,splitAmt,0,'')
      typeSpinSplitTemp(typeID) = splitAmt
   enddo

   ! Now we pass the parameters we just read in to the Potential module so that
   !   it can initialize the values it stores.
   call setPotControlParameters(fbLevel,lastIt,xcCodeTemp,rlxFact,cTest,&
         & plusUJFormTemp,numPlusUJItemsTemp,plusUJItemIDTemp,&
         & plusUJItemValueTemp,typeSpinSplitTemp)

   deallocate (plusUJItemValueTemp)
   deallocate (plusUJItemIDTemp)
   deallocate (typeSpinSplitTemp)

end subroutine readMainControl


subroutine readDOSControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Constants, only: hartree
   use O_ReadDataSubs

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readAndCheckLabel(readUnit,writeUnit,len('DOS_INPUT_DATA'),&
                              'DOS_INPUT_DATA')

   ! Read all the input data.
   call readData(readUnit,writeUnit,deltaDOS,sigmaDOS,0,'')
   call readData(readUnit,writeUnit,eminDOS,emaxDOS,0,'')
   call readData(readUnit,writeUnit,detailCodePDOS,0,'')

   ! Apply necessary conversions to the units.  (eV -> hartree atomic units)
   deltaDOS = deltaDOS / hartree
   sigmaDOS = sigmaDOS / hartree / (2.0_double * sqrt(2.0_double *&
         & log(2.0_double)))
   eminDOS  = eminDOS  / hartree
   emaxDOS  = emaxDOS  / hartree

end subroutine readDOSControl


subroutine readBondControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Constants, only: hartree, bohrRad
   use O_ReadDataSubs

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readAndCheckLabel(readUnit,writeUnit,len('BOND_INPUT_DATA'),&
         & 'BOND_INPUT_DATA')

   ! Read the data.
   call readData(readUnit,writeUnit,maxBL,0,'')
   call readData(readUnit,writeUnit,deltaBOND,sigmaBOND,0,'')
   call readData(readUnit,writeUnit,eminBOND,emaxBOND,0,'')
   call readData(readUnit,writeUnit,outputCodeBondQ,0,'')
   call readData(readUnit,writeUnit,doBond3C,0,'')
   call readData(readUnit,writeUnit,maxNumNeighbors,0,'')

   ! Apply necessary conversions to the units.  (Angstroms to bohrRadii a.u.
   !   and eV to hartree atomic units.)
   maxBL     = maxBL     / bohrRad
   deltaBOND = deltaBOND / hartree
   sigmaBOND = sigmaBOND / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))
   eminBOND  = eminBOND  / hartree
   emaxBOND  = emaxBOND  / hartree

end subroutine readBondControl


subroutine readSYBDControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_ReadDataSubs
   use O_KPoints, only: readSYBDKPoints

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   call readSYBDKPoints(readUnit, writeUnit)

end subroutine readSYBDControl


subroutine readPACSControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Kinds
   use O_Constants, only: hartree
   use O_CommandLine, only: excitedQN_n, excitedQN_l
   use O_ReadDataSubs

   ! Make sure no unwanted variables are declared.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

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
   call readAndCheckLabel(readUnit,writeUnit,len('PACS_INPUT_DATA'),&
         & 'PACS_INPUT_DATA')

   call readData(readUnit,writeUnit,excitedAtomPACS,0,'')
   call readData(readUnit,writeUnit,deltaPACS,sigmaPACS,0,'')
   call readData(readUnit,writeUnit,onsetEnergySlackPACS,&
         & energyWindowPACS,0,'')
   call readData(readUnit,writeUnit,numCorePACS,0,'')
   do i = 1, numCorePACS
      call readData(readUnit,writeUnit,5,tempArray,0,'')
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
   sigmaPACS = sigmaPACS / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))
   onsetEnergySlackPACS = onsetEnergySlackPACS / hartree
   energyWindowPACS     = energyWindowPACS / hartree
   totalEnergyDiffPACS  = totalEnergyDiffPACS / hartree


end subroutine readPACSControl


subroutine readOptcControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Constants, only: hartree
   use O_ReadDataSubs

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Read input OPTC control data.
   call readAndCheckLabel(readUnit,writeUnit,len('OPTC_INPUT_DATA'),&
         & 'OPTC_INPUT_DATA')
   call readData(readUnit,writeUnit,cutoffEnOPTC,0,'')
   call readData(readUnit,writeUnit,maxTransEnOPTC,0,'')
   call readData(readUnit,writeUnit,deltaOPTC,0,'')
   call readData(readUnit,writeUnit,sigmaOPTC,0,'')

   ! Apply necessary conversions to a.u.
   cutoffEnOPTC   = cutoffEnOPTC / hartree
   maxTransEnOPTC = maxTransEnOPTC / hartree
   deltaOPTC = deltaOPTC / hartree
   sigmaOPTC = sigmaOPTC / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))

end subroutine readOptcControl


subroutine readSigeControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Constants, only: hartree
   use O_ReadDataSubs

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Read the input sigma(E) control variables.
   call readAndCheckLabel(readUnit,writeUnit,len('SIGE_INPUT_DATA'),&
         & 'SIGE_INPUT_DATA')
   call readData(readUnit,writeUnit,cutoffEnSIGE,0,'')
   call readData(readUnit,writeUnit,maxTransEnSIGE,0,'')
   call readData(readUnit,writeUnit,deltaSIGE,0,'')
   call readData(readUnit,writeUnit,sigmaSIGE,0,'')

   ! Apply necessary conversions to a.u.
   cutoffEnSIGE   = cutoffEnSIGE / hartree
   maxTransEnSIGE = maxTransEnSIGE / hartree
   deltaSIGE = deltaSIGE / hartree
   sigmaSIGE = sigmaSIGE / hartree /(2.0_double*sqrt(2.0_double * &
         & log(2.0_double)))

end subroutine readSigeControl


subroutine readWaveControl(readUnit,writeUnit)

   ! Use necessary modules
   use O_Lattice
   use O_Constants, only: hartree
   use O_ReadDataSubs

   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Read the input wave function / charge density plotting control variables.
   call readAndCheckLabel(readUnit,writeUnit,len('WAVE_INPUT_DATA'),&
         & 'WAVE_INPUT_DATA')
   call readNumMeshPoints(readUnit,writeUnit)
   call readData(readUnit,writeUnit,eminWAVE,emaxWAVE,0,'')
   call readData(readUnit,writeUnit,doRho,0,'')
   call readData(readUnit,writeUnit,styleWAVE,0,'')

   ! Apply the necessary conversions of the data to atomic units.
   eminWAVE = eminWAVE / hartree
   emaxWAVE = emaxWAVE / hartree

end subroutine readWaveControl


end module O_Input
