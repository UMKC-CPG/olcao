# This python script is part of the OLCAO python module
# Author: James currie <jecyrd@mail.umkc.edu>

import os
import sys

############# Data structures for organizational purposes. ###############
# Information important to the OLCAO Environment
class constants():
  def __init__(self):
    self.eVhartree = 27.21138505

class olcaoEnv():
  def __init__(self):
    self.obranch = os.environ.get('BRANCH')
    self.odir = os.environ.get('OLCAO_DIR')
    self.obin = os.environ.get('OLCAO_BIN')
    self.odata = os.environ.get('OLCAO_DATA')
    self.otemp = os.environ.get('OLCAO_TEMP')
    self.otar = os.environ.get('OLCAO_TAR')
    self.otouch = os.environ.get('OLCAO_TOUCH')
    self.odouble = os.environ.get('OLCAO_DOUBLE')
    self.opbs = os.environ.get('OLCAO_PBS')
    self.olsf = os.environ.get('OLCAO_LSF')
    self.obash = os.environ.get('OLCAO_BASH')
    self.ovaspdir = os.environ.get('VASPPOT_DIR')

    if (self.obin == None):
      print 'Did you forget to source the ".olcao" file?'
      sys.exit()

  def reset(self):
    self.__init__()

# Execution Control
class exeCtrl():
  def __init__(self):
    self.niceValue = 0
    self.test = 0
    self.edge = 'gs' 
    self.QN_n = None
    self.QN_l = None
    self.altSCF = 0
    self.altPSCF = 0
    self.spinPol = 0
    self.serialxyz = 0
  def reset(self):
    self.__init__()

# Information about which job to run
class jobInfo():
  def __init__(self):
    self.dosRun =     0
    self.pacsRun =    0
    self.sybdRun =    0
    self.bondRun =    0
    self.optcRun =    0
    self.sigeRun =    0
    self.waveRun =    0
    self.dosRunSCF =  0
    self.bondRunSCF = 0
    self.waveRunSCF = 0
  def reset(self):
    self.__init__()

# Basis Information
class basisInfo():
  def __init__(self):
    self.scfBasis      = 'fb'
    self.scfBasisSet   = 0
    self.scfBasisCode  = None
    self.pscfBasis     = 'no'
    self.pscfBasisSet  = 0
    self.pscfBasisCode = None
  def reset(self):
    self.__init__()

# Complete file and directory names
class filesAndDirectories():
  def __init__(self):
    self.OLCAOlock    = 'OLCAOlock'
    self.OLCAOkill    = 'OLCAOkill'
    self.runtime      = 'runtime'
    self.intermediate = 'intermediate'
    self.energy       = 'enrg'
    self.iteration    = 'iter'
    self.kp_scf       = 'kp-scf'
    self.kp_pscf      = 'kp-pscf'
    self.kp_alt       = 'kp-alt'
    self.moment       = 'mgmo'
    self.initPot      = 'scfV'
    self.structure    = 'structure'
    self.olcao        = 'olcao'
    self.atomPos      = 'atomPos'
    self.lattice      = 'lattice'
    
    self.proj_home    = None
    self.inputs       = None
  def reset(self):
    self.__init__()

# File nature components
class fileNatureComps():
  def __init__(self):
    self.setup = 'setup'
    self.main  = 'main'
    self.intg  = 'intg'
    self.band  = 'band'
    self.sybd  = 'sybd'
    self.dos   = 'dos'
    self.bond  = 'bond'
    self.optc  = 'optc'
    self.pacs  = 'pacs'
    self.sige  = 'sige'
    self.wave  = 'wave'
  def reset(self):
    self.__init__()

# Auxiliary suffix components
class auxSuffixComps():
  def __init__(self):
    self.loci       = '.loci'
    self.tdos       = '.t'
    self.pdos       = '.p'
    self.elf        = '.elf'
    self.epsl       = '.eps1'
    self.epsli      = '.eps1i'
    self.eps2       = '.eps2'
    self.cond       = '.cond'
    self.ref        = '.ref'
    self.pot        = '.pot'
    self.valeRho    = '.rhoV'
    self.up         = '.up'
    self.dn         = '.dn'
    self.upPdn      = '.up+dn'
    self.upMdn      = '.up-dn'
    self.minusNeut  = '-N'
    self.profile    = '.prof'
  def reset(self):
    self.__init__()

# File Name Extensions
class fnExtensions():
  def __init__(self):
    self.plot = '.plot'
    self.raw  = '.raw'
    self.out  = '.out'
    self.dat  = '.dat'
    self.alt  = '.alt'
    self.dx   = '.dx'
    self.hdf  = '.hdf5'
  def reset(self):
    self.__init__()

# File Name Extensions
class execNames():
  def __init__(self):
    self.exe = ['' for x in xrange(8)]
    self.exe[0] = 'OLCAOsetup'
    self.exe[1] = 'OLCAOmain'
    self.exe[2] = 'OLCAOintg'
    self.exe[3] = 'OLCAOband'
    self.exe[4] = 'OLCAOdos'
    self.exe[5] = 'OLCAObond'
    self.exe[6] = 'OLCAOoptc'
    self.exe[7] = 'OLCAOwave'
    self.exeMechanism = None
  def reset(self):
    self.__init__()

# All encompassing object
class olcaoDS():
  def __init__(self):
    # Setup Structures
    self.env = olcaoEnv()
    self.control = exeCtrl()
    self.jobInfo = jobInfo()
    self.basisInfo = basisInfo()
    self.ioNames = filesAndDirectories()
    self.fileNature = fileNatureComps()
    self.auxSuffix = auxSuffixComps()
    self.fExtensions = fnExtensions()
    self.exeNames = execNames()
    self.consts = constants()

#    # Get current dir to go back to later, and cd to jobdir
#    # This will be usefull for parallel version
#    os.chdir(jobdir)
#
#    # Parse Command Line
#    self.parseCLP(clp)
#
#    # Setup Executable Environment
#    self.getTempDir()
#    self.initDirectories()
#    self.initExes()
#    self.setEdgeCode()
#
#    # Finally set the appropriate spinPol value for the
#    # XC_Code parameter in the olcao.dat
#    self.XC_CodeSpin()

  # Subroutine to parse all the command line options
  # clp is list holding the command line parameters
  def parseCLP(self,clp):
    numParams = len(clp)

    # Initialize a counter
    counter = 0
    
    # Loop over and handle command line parameters
    while (counter < numParams):
      if (clp[counter] == '-scf'):
        self.basisInfo.scfBasisSet=1
        counter += 1
        self.basisInfo.scfBasis = clp[counter].lower()

      elif (clp[counter] == '-pscf'):
        self.basisInfo.pscfBasisSet=1
        counter += 1
        self.basisInfo.pscfBasis = clp[counter].lower()

      elif (clp[counter] == '-scfaltkp'):
        self.control.altSCF = 1

      elif (clp[counter] == '-pscfaltkp'):
        self.control.altPSCF = 1

      elif (clp[counter] == '-serialxyz'):
        self.control.serialxyz = 1

      elif (clp[counter] == '-pacs'):
        self.jobInfo.pacsRun = 1
        counter += 1
        self.control.edge = clp[counter]
        if (self.basisInfo.pscfBasisSet == 0):
          self.basisInfo.pscfBasis = 'eb'
          self.basisInfo.pscfBasisSet = 1

      elif (clp[counter] == '-sybd'):
        self.jobInfo.sybdRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.pscfBasisSet == 0):
          self.basisInfo.pscfBasis = 'fb'
          self.basisInfo.pscfBasisSet = 1

      elif (clp[counter] == '-dos'):
        self.jobInfo.dosRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        print "Basis Set: ",self.basisInfo.pscfBasisSet
        if (self.basisInfo.pscfBasisSet == 0):
          self.basisInfo.pscfBasis = 'fb'
          self.basisInfo.pscfBasisSet = 1

      elif (clp[counter] == '-scfdos'):
        self.jobInfo.dosRunSCF = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]

      elif (clp[counter] == '-bond'):
        self.jobInfo.bondRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.pscfBasisSet == 0):
          self.basisInfo.pscfBasis = 'mb'
          self.basisInfo.pscfBasisSet = 1

      elif (clp[counter] == '-scfbond'):
        self.jobInfo.bondRunSCF = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.scfBasisSet == 0):
          self.basisInfo.scfBasis = 'mb'
          self.basisInfo.scfBasisSet = 1

      elif (clp[counter] == '-optc'):
        self.jobInfo.optcRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.scfBasisSet == 0):
          self.basisInfo.scfBasis = 'eb'
          self.basisInfo.scfBasisSet = 1
    
      elif (clp[counter] == '-sige'):
        self.jobInfo.sigeRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.scfBasisSet == 0):
          self.basisInfo.scfBasis = 'fb'
          self.basisInfo.scfBasisSet = 1
      
      elif (clp[counter] == '-wave'):
        self.jobInfo.waveRun = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]
        if (self.basisInfo.scfBasisSet == 0):
          self.basisInfo.scfBasis = 'fb'
          self.basisInfo.scfBasisSet = 1
      
      elif (clp[counter] == '-wavescf'):
        self.jobInfo.waveRunSCF = 1
        if ((counter == numParams-1) or ('-' in clp[counter+1])):
            self.control.edge = 'gs'
        else:
          counter += 1
          self.control.edge = clp[counter]

      elif (clp[counter] == "-nice"):
        counter += 1
        self.control.niceValue = clp[counter]

      elif (clp[counter] == "-test"):
        counter += 1
        self.control.test = clp[counter]

      else:
        print "UNKNOWN COMMAND LINE PARAMETER ", clp[counter], \
            " ABORTING.\n"
        sys.exit()
      # Bottom of loop increment counter
      counter += 1
  
  def getTempDir(self):
    # Get and separate the current directory. 
    # Then pop off the first value of the list because it'll 
    # just be an empty string
    directories = os.getcwd().split('/')
    directories.pop(0)
  
    # Find the location in the list that contains the user name
    unameloc = directories.index(os.getenv('USER'))
  
    # Append the remaining directories to the temp directory
    for x in range(unameloc+1,len(directories)):
      self.env.otemp += '/' + directories[x]
  
  # Initialize the directory names and directory structures
  def initDirectories(self):
    # Set project home and inputs directory paths
    self.ioNames.proj_home = os.getcwd()
    self.ioNames.inputs = self.ioNames.proj_home + '/inputs'
  
    # Create the temp directory if it does not already exist
    if (not os.path.exists(self.env.otemp)):
        systring = self.env.obin + '/mkdirhier ' + self.env.otemp
        os.system(systring)
  
    # If a link to the temp directory does not yet exist, then create one.
    # If a link already exists, but points to the wrong place, then rename
    # the existing link and create a new one.
    if (not os.path.exists(self.ioNames.intermediate)):
      systring = 'ln -s ' + self.env.otemp + ' ' + self.ioNames.intermediate
      os.system(systring)
    else:
      intermediateLoc = os.path.realpath(self.ioNames.intermediate)
      print intermediateLoc
      print self.env.otemp
      if (intermediateLoc != self.env.otemp):
        # Move intermediate directory
        systring = 'mv ' + self.ioNames.intermediate 
        systring += ' ' + self.ioNames.intermediate + 'FIXME'
        os.system(systring)
        # Create new intermediate directory
        systring = 'ln -s ' + self.env.otemp + ' ' 
        systring += self.ioNames.intermediate
        os.system(systring)

  def checkGammaKP(self,fn):
    kpfile = open(fn,'r')
    lines = [x for x in kpfile.readlines()]
    kpfile.close()

    numkp = int(lines[1])
    splist = lines[-1].split()
    kpvals = [float(splist[-3]),float(splist[-2]),float(splist[-1])]
    
    if (numkp==1):
      if ((kpvals[0] == 0) and (kpvals[1] == 0) and (kpvals[2] == 0)):
        return 1

    # Else either of those conditions above are false
    return 0


  def initExes(self):
    # Check each kpoint input file to see if the gamma kp is needed
    # anywhere in the calculation
    fstring = self.ioNames.inputs + '/' + self.ioNames.kp_scf
    fstring += self.fExtensions.dat
    altGamma = self.checkGammaKP(fstring)
   
    # If the alternate kpoint set is used for either the scf or
    # pscf then use the altGamma flag for the scf or pscf.
    # Else use information in file to set flag.
    if (self.control.altSCF == 1):
      scfGamma = altGamma
    else:
      fstring = self.ioNames.inputs + '/' + self.ioNames.kp_scf
      fstring += self.fExtensions.dat
      scfGamma = self.checkGammaKP(fstring)
    if (self.control.altPSCF == 1):
      pscfGamma = altGamma
    else:
      fstring = self.ioNames.inputs + '/' + self.ioNames.kp_scf
      fstring += self.fExtensions.dat
      pscfGamma = self.checkGammaKP(fstring)

    # Determine if the calculation is going to use the gamma kpoint
    # and prepend "g" onto the list of executables
    if (scfGamma == 1):
      # Setup and Main comprise the SCF Portion
      self.exeNames.exe[0] = 'g' + self.exeNames.exe[0]
      self.exeNames.exe[1] = 'g' + self.exeNames.exe[1]
    if (scfGamma == 1):
      # All other executables comprise the PSCF portion
      for x in range(2,len(self.exeNames.exe)):
        self.exeNames.exe[x] = 'g' + self.exeNames.exe[x]


    # Determine if the test exectuable should be used. If so, then
    # append the test string
    if (self.control.test == 1):
      # All other executables comprise the PSCF portion
      for x in xrange(len(self.exeNames.exe)):
        self.exeNames.exe[x] = 'test' + self.exeNames.exe[x]

    # Determine the command lin execution mechanism
    if ((self.control.niceValue > 19) or (self.control.niceValue < -20)):
      self.exeNames.exeMechanism = "time"
    else:
      self.exeNames.exeMechanism = "nice -n " \
                                  + str(self.control.niceValue) \
                                  + " time"

  def setEdgeCode(self):
    if (self.control.edge == 'gs'):
      self.control.QN_n = '0'
      self.control.QN_l = '0'
    elif (self.control.edge == '1s'):
      self.control.QN_n = '1'
      self.control.QN_l = '0'
    elif (self.control.edge == '2s'):
      self.control.QN_n = '2'
      self.control.QN_l = '0'
    elif (self.control.edge == '2p'):
      self.control.QN_n = '2'
      self.control.QN_l = '1'
    elif (self.control.edge == '3s'):
      self.control.QN_n = '3'
      self.control.QN_l = '0'
    elif (self.control.edge == '3p'):
      self.control.QN_n = '3'
      self.control.QN_l = '1'
    elif (self.control.edge == '3d'):
      self.control.QN_n = '3'
      self.control.QN_l = '2'
    elif (self.control.edge == '4s'):
      self.control.QN_n = '4'
      self.control.QN_l = '0'
    elif (self.control.edge == '4p'):
      self.control.QN_n = '4'
      self.control.QN_l = '1'
    elif (self.control.edge == '4d'):
      self.control.QN_n = '4'
      self.control.QN_l = '2'
    elif (self.control.edge == '4f'):
      self.control.QN_n = '4'
      self.control.QN_l = '3'
    elif (self.control.edge == '5s'):
      self.control.QN_n = '5'
      self.control.QN_l = '0'
    elif (self.control.edge == '5p'):
      self.control.QN_n = '5'
      self.control.QN_l = '1'
    elif (self.control.edge == '5d'):
      self.control.QN_n = '5'
      self.control.QN_l = '2'
    elif (self.control.edge == '6s'):
      self.control.QN_n = '6'
      self.control.QN_l = '0'
    elif (self.control.edge == '6p'):
      self.control.QN_n = '6'
      self.control.QN_l = '1'
    elif (self.control.edge == '7s'):
      self.control.QN_n = '7'
      self.control.QN_l = '0'

    # Determine the number code for which basis set to use for scf 
    if (self.basisInfo.scfBasis == 'mb'):
      self.basisInfo.scfBasisCode = '1'
    elif (self.basisInfo.scfBasis == 'fb'):
      self.basisInfo.scfBasisCode = '2'
    elif (self.basisInfo.scfBasis == 'eb'):
      self.basisInfo.scfBasisCode = '3'
    
    # Determine the number code for which basis set to use for pscf 
    if (self.basisInfo.pscfBasis == 'mb'):
      self.basisInfo.pscfBasisCode = '1'
    elif (self.basisInfo.pscfBasis == 'fb'):
      self.basisInfo.pscfBasisCode = '2'
    elif (self.basisInfo.pscfBasis == 'eb'):
      self.basisInfo.pscfBasisCode = '3'

  # This subroutine has the purpose of figuring out what the value for
  # self.control.spinPol should be based on the XC_CODE given in
  # the olcao.dat file
  def XC_CodeSpin(self):
    # Read in the olcao.dat file and find the xc_code line
    fn = self.ioNames.inputs+'/' \
        + self.ioNames.olcao \
        + self.fExtensions.dat
    olcaodat = open(fn,'r')
    lines = [x for x in olcaodat.readlines()]

    found = 0
    for x in xrange(len(lines)):
      if 'XC_CODE' in lines[x]:
        found = x+1
        break

    if (found == 0):
      print "Could not find XC_CODE defaulting to 100\n" \
            + "Spinpol defaulting to 1\n"
      self.control.spinPol = 1

    else:
      # Else read in the xc_code.dat file and find the appropriate
      # spinPol value
      xc_code = int(lines[found])

      fn = self.env.odata + '/xc_code.dat'
      xcfile = open(fn,'r')
      lines = [x for x in xcfile.readlines()]
      xcfile.close()
     
      for x in xrange(len(lines)):
        # First split the values
        splist = lines[x].split()
        if ((len(splist)>1) and (xc_code == splist[1])):
          self.control.spinPol = splist[2]
          break


