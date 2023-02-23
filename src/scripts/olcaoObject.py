# This python script is part of the OLCAO python module
# Author: James Currie <jecyrd@mail.umkc.edu>

import os
import sys
import time
from olcaoStructures import *

class olcaoJob(olcaoDS):
  def __init__(self,jobdir,clp):
    olcaoDS.__init__(self)
    self.jobdir = jobdir
    self.clp = clp
    self.origDir = os.getcwd() 

    # Empty object later to be used for the runtime file.
    self.runtime = None

  # This method removes the lock file, signally the end of the run.
  # Also prints to the runtime file that the program sequence is complete.
  def cleanUp(self):
    syscmd = 'rm ' + self.env.otemp + '/' + self.ioNames.OLCAOlock
    os.system(syscmd)

  # This method creates the runtime file object
  def openRuntime(self):
    self.runtime = open('runtime ','w')

  # This method closes the runtime file object
  def closeRuntime(self):
    self.runtime.close()

  # This subroutine makes a lock file to signal to other OLCAO jobs
  # that an olcao job is currently running in this directory.
  def createLockFile(self):
    lockfile = self.env.otemp + '/' + self.ioNames.OLCAOlock 
    if (os.path.isfile(lockfile)):
      print "Lock File found in " + self.env.otemp + '\n' \
            + 'Is another olcao script runnning?\n' \
            + 'Did an olcao script die badly?\n'
      sys.exit()
    else:
      os.system("touch " + lockfile)

  # This subroutine initializes the olcaoDataStructure member to this
  # class. It is neveccary to run in order to pull in all the
  # environment information.
  def setupOLCAOenv(self):
    # Parse Command Line
    self.parseCLP(self.clp)

    # Setup Executable Environment
    self.getTempDir()
    self.initDirectories()
    self.initExes()
    self.setEdgeCode()

    # Finally set the appropriate spinPol value for the XC_Code
    # parameter in the olcao.dat file.
    self.XC_CodeSpin()

  # This method is for running the various steps specified in the
  # various number of inputs
  def runOLCAO(self):
    # Calculate the setup portion of the program. This is independent
    # of the edge.
    if (self.basisInfo.scfBasis != "no"):
      self.runTask(self.fileNature.setup, self.basisInfo.scfBasis,"")

    if (self.jobInfo.pacsRun == 1):
      # We are running a pacs calculation.
      self.runTask(self.fileNature.main, self.basisInfo.scfBasis,"gs")
      self.runTask(self.fileNature.intg, self.basisInfo.pscfBasis,"gs")
      self.runTask(self.fileNature.band, self.basisInfo.pscfBasis,"gs")

      if (self.control.edge != "gs"):
        self.runTask(self.fileNature.main, self.basisInfo.scfBasis, \
            self.control.edge)

        # Now that both main calculations are done for the PACS calculation
        # We can update the pacs input file with the total energy
        # difference and the energy range.
        self.updatePACS(self.control.edge, self.basisInfo.scfBasis, \
            self.basisInfo.pscfBasis)
      
        self.runTask(self.fileNature.intg, self.basisInfo.pscfBasis, \
            self.control.edge)
        self.runTask(self.fileNature.band, self.basisInfo.pscfBasis, \
            self.control.edge)
        self.runTask(self.fileNature.pacs, self.basisInfo.pscfBasis, \
            self.control.edge)
    else:
      # We are doing some other caluclation
      # i.e. PDOS, SYBD, BOND, OPTC, SIGE, WAVE
      if (self.basisInfo.scfBasis != "no"):
        self.runTask(self.fileNature.main, self.basisInfo.scfBasis, \
            self.control.edge)
      if (self.basisInfo.pscfBasis != "no"):
        self.runTask(self.fileNature.intg, self.basisInfo.pscfBasis, \
            self.control.edge)
        if (self.jobInfo.sybdRun == 1):
          self.runTask(self.fileNature.sybd, self.basisInfo.pscfBasis, \
             self.control.edge)
        else:
          self.runTask(self.fileNature.band, self.basisInfo.pscfBasis, \
             self.control.edge)
          if (self.jobInfo.dosRun ==1):
            self.runTask(self.fileNature.dos, self.basisInfo.pscfBasis, \
               self.control.edge)
          elif (self.jobInfo.bondRun ==1):
            self.runTask(self.fileNature.bond, self.basisInfo.pscfBasis, \
               self.control.edge)
          elif (self.jobInfo.optcRun ==1):
            self.runTask(self.fileNature.optc, self.basisInfo.pscfBasis, \
               self.control.edge)
          elif (self.jobInfo.sigeRun ==1):
            self.runTask(self.fileNature.sige, self.basisInfo.pscfBasis, \
               self.control.edge)
          elif (self.jobInfo.waveRun ==1):
            self.runTask(self.fileNature.wave, self.basisInfo.pscfBasis, \
               self.control.edge)

  # This subroutine sets up the olcao environment, creates a lock file
  # and then runs the olcao job. If different steps are wanted for
  # whatever reason use setupOLCAOenv, runOLCAO, and createLockFile, 
  # methods to this object
  def run(self):
    # First setup the OLCAO environment
    self.setupOLCAOenv()

    # First create a lock file
    self.createLockFile()

    # Create the runtime file object
    self.openRuntime()

    # Finally run the OLCAO job
    self.runOLCAO()

    # Close the runtime file object
    self.closeRuntime()

    # Clean up
    self.cleanUp()

  # A generic run task subroutine. This subroutine will run a requested
  # task according to the input code it is given.
  def runTask(self, taskName, taskBasis, taskEdge):
    taskBasis = '-' + taskBasis # Prepend a hyphen
    taskEdge += '_' # Append an underscore

    # If running setup, correct the taskEdge since the edge doesn't matter
    if (taskName == self.fileNature.setup):
      taskEdge = ''

    # Change the working directory to the intermediate location
    os.chdir(self.env.otemp)

    # Mark the date and time of the begining of this task in runtime
    runstr = taskEdge + taskName + taskBasis
    self.runtime.write(runstr + '\n')
    runstr = time.strftime("%d-%m-%Y)") + '\n'
    runstr += time.strftime("%H:%M:%S")
    self.runtime.write(runstr + '\n')

    # Identify the task kpoint file to use.
    taskKP = self.getTaskKP(taskName)

    # Determine if the alternate KPoints were used and if so assign
    # a value to the alt string so that the outputs etc created with it
    # can be kept separate.
    if ('kp-alt' in taskKP):
      alt = ".alt"
    else:
      alt = ""
    
    # Determine if the output for this task already exists.
    checkfile = self.ioNames.proj_home + '/' \
        + taskEdge + taskName + taskBasis + alt + self.fExtensions.out
    if (not os.path.isfile(checkfile)):
      # Mark the lock file with the current exectuable.
      sysstr = "echo " + taskName + taskBasis + " | cat > " \
          + self.env.otemp + '/' + self.ioNames.OLCAOlock
      os.system(sysstr)

      # Check to see if the necessary input files are in the proj_home
      # directory. If they are not, then try to copy them from the input
      # directory into the proj_home directory. This is done so that a 
      # "pure" version of the input exists in the inputs directory while a
      # "user modifiable" version of the input exists in the proj_home
      # directory
      self.checkInput(taskName,taskBasis,taskEdge,taskKP)

      # Copy the task input files from the proj_home directory to temp
      self.copyInput(taskName,taskBasis,taskEdge,taskKP,alt)

      # Create the proper command line parameters (CLP) for this task.
      taskCLP = self.determineCLP(taskName,taskEdge)
      #print "Task Command Line: " + taskName + ' ' + taskCLP + '\n'

      # Execute the requested task program and save the runtime data.
      self.executeProgram(taskName, taskCLP)

      # Manage the output files
      self.manageOutput(taskName,taskBasis,taskEdge,taskKP,alt)

    else:
      runstr = '\n' + taskEdge + taskName + taskBasis \
              + self.fExtensions.out + 'file already exists.\n' \
              + "Skipping " + taskEdge + taskName + taskBasis \
              + 'execution.\n\n'
      self.runtime.write(runstr)
    if (int(self.env.otouch) == 1):
      os.system('touch ' + self.env.otemp + '/*')

    # Mark the date and time of the ending of this task.
    runstr = 'End: ' + time.strftime('%H:%M:%S]n') 
    self.runtime.write(runstr + '\n\n')

    # Check to see if OLCAOkill file exists. If it does, kill execution
    # of the script.
    if (os.path.isfile(self.ioNames.OLCAOkill)):
      sys.exit()

  # Obtainn the name of the kpoint file to use for the current task.
  def getTaskKP(self, taskName):
    if ((taskName == self.fileNature.setup) or \
              (taskName == self.fileNature.main)):
      if (self.control.altSCF == 0):
        taskKP = self.ioNames.kp_scf
      else:
        taskKP = self.ioNames.kp_alt
    else:
      if (self.control.altPSCF == 0):
        taskKP = self.ioNames.kp_pscf
      else:
        taskKP = self.ioNames.kp_alt

    return taskKP

  def checkInput(self,taskName,taskBasis,taskEdge,taskKP):
    # Check for the primary "olcao" input file, and the system structure.
    datcheck = self.ioNames.proj_home + '/'
    datcheck += self.ioNames.olcao + self.fExtensions.dat

    structcheck = self.ioNames.proj_home + '/'
    structcheck += self.ioNames.structure + self.fExtensions.dat
    if (not os.path.isfile(datcheck)):
      syscall = 'cp ' + self.ioNames.inputs + '/'
      syscall += self.ioNames.olcao + self.fExtensions.dat
      syscall += ' ' + self.ioNames.proj_home + '/'

      os.system(syscall)
    if (not os.path.isfile(structcheck)):
      syscall = 'cp ' + self.ioNames.inputs + '/'
      syscall += self.ioNames.structure + self.fExtensions.dat
      syscall += ' ' + self.ioNames.proj_home + '/'

      os.system(syscall)

    # Check for the set of kpoints needed for this task.
    kpcheck = self.ioNames.proj_home + '/' + taskKP + self.fExtensions.dat
    if (not os.path.isfile(kpcheck)):
      syscall = 'cp ' + self.ioNames.inputs + '/' + taskKP 
      syscall += self.fExtensions.dat + ' ' + self.ioNames.proj_home + '/'

      os.system(syscall)

    # Check for a potential coefficient file for main.
    if (taskName == self.fileNature.main):
      potcheck = self.ioNames.proj_home + '/'
      potcheck += self.ioNames.initPot + self.fExtensions.dat
      if (not os.path.isfile(potcheck)):
        syscall = 'cp ' + self.ioNames.inputs + '/' + self.ioNames.initPot
        syscall += self.fExtensions.dat + ' ' + self.ioNames.proj_home + '/'

        os.system(syscall)

  # This method will copy all of the data files for each type of task
  def copyInput(self,taskName,taskBasis,taskEdge,taskKP,alt):
    # Define the temporary basis tag
    tempBasis = '-temp'

    # Make some local variables because the syntax can be cumbersome
    proj_home = self.ioNames.proj_home
    dat = self.fExtensions.dat
    olcaodat = self.ioNames.olcao + dat
    structuredat = self.ioNames.structure + dat

    os.system('cp ' + proj_home + '/' + olcaodat + ' fort.5')
    os.system('cp ' + proj_home + '/' + structuredat + ' fort.4')
    os.system('cp ' + proj_home + '/' + taskKP + dat + ' fort.15')

    # Get the appropriate hdf5 stored intermediate data.
    if (taskName == self.fileNature.main):
      sysstr = 'mv ' + self.fileNature.setup + taskBasis
      sysstr += alt + self.fExtensions.hdf +' ' 
      sysstr += self.fileNature.setup + tempBasis + self.fExtensions.hdf
      os.system(sysstr)
      self.getPotCoeffs(taskName,taskBasis,taskEdge)
    elif (taskName == self.fileNature.intg):
      self.getPotCoeffs(taskName,taskBasis,taskEdge)

    elif(taskName==self.fileNature.band or taskName==self.fileNature.sybd):
      sysstr = 'mv ' +  taskEdge + self.fileNature.intg + taskBasis 
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.fileNature.intg + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

    elif(taskName==self.fileNature.dos or taskName==self.fileNature.bond):
      sysstr = 'mv ' +  taskEdge + self.fileNature.band + taskBasis
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.fileNature.band + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

    elif(taskName==self.fileNature.optc or taskName==self.fileNature.sige):
      sysstr = 'mv ' + taskEdge + self.fileNature.band + taskBasis
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.filenature.band + taskBasis + self.fExtensions.hdf
      os.system(sysstr)

      sysstr = 'mv ' +  taskEdge + self.fileNature.intg + taskBasis 
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.fileNature.intg + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

    elif (taskName == self.fileNature.pacs):
      sysstr = 'mv gs_' + selfyfileNature.intg + taskBasis 
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.fileNature.intg + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

      sysstr = 'mv gs_' + selfyfileNature.band + taskBasis 
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.fileNature.band + tempBasis + self.fExtensions.hdf
      os.system(sysstr)
      
      sysstr = 'mv ' + taskEdge + self.fileNature.band + taskBasis
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.filenature.band + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

    elif (taskName == self.fileNature.wave):
      sysstr = 'mv ' + taskEdge + self.fileNature.band + taskBasis
      sysstr += alt + self.fExtensions.hdf + ' '
      sysstr += self.filenature.band + tempBasis + self.fExtensions.hdf
      os.system(sysstr)

      self.getPotCoeffs(taskName,taskBasis,taskEdge)

  def getPotCoeffs(self, taskName, taskBasis, taskEdge):
    
    # Establish the order of alternate basis sets to check for a set of
    # potential coefficients
    altBasis = ['-eb','-fb','-eb']

    # Setup does not need any set of potential coefficients
    if (taskName == self.fileNature.setup):
      return
    
    thing = taskEdge + self.ioNames.initPot + taskBasis + \
        self.fExtensions.dat

    # Check for the potential from the same edge and basis of the current
    # task.
    if (os.path.isfile(self.ioNames.proj_home+'/'+thing)):
      self.copy(self.ioNames.proj_home+'/'+thing,'fort.8')
      self.runtime.write("Using " + thing + '\n')
      return

    # Check each of the alternate basis sets for the current edge.
    for basis in altBasis:
      if (os.path.isfile(self.ioNames.proj_home+'/'+thing)):
        self.copy(self.ioNames.proj_home+'/'+thing,'fort.8')
        self.runtime.write("Using " + thing + '\n')
        return

    # Check the ground state of the current basis.
    thing = 'gs_' + self.ioNames.initPot + taskBasis + self.fExtensions.dat
    if (os.path.isfile(self.ioNames.proj_home+'/'+thing)):
      self.copy(self.ioNames.proj_home+'/'+thing,'fort.8')
      self.runtime.write("Using " + thing + '\n')
      return

    # Check each of the alternate basis sets for the ground state
    for basis in altBasis:
      if (os.path.isfile(self.ioNames.proj_home+'/'+thing)):
        self.copy(self.ioNames.proj_home+'/'+thing,'fort.8')
        self.runtime.write("Using " + thing + '\n')
        return

    # Check for the initial non-scf potential int he project home direct.
    thing = self.ioNames.initPot + self.fExtensions.dat 
    if (os.path.isfile(self.ioNames.proj_home+'/'+thing)):
      self.copy(self.ioNames.proj_home+'/'+thing,'fort.8')
      self.runtime.write("Using " + thing + '\n')
      return

    # Check for the initial non-scf potential in the inputs directory
    if (os.path.isfile(self.ioNames.inputs+'/'+thing)):
      self.copy(self.ioNames.inputs+'/'+thing,'fort.8')
      self.runtime.write("Using " + thing + '\n')
      return

  # This subroutine will determine what the command line parameters for
  # the task execution should be.
  def determineCLP(self,taskName,taskEdge):
    # Initialize empty string to return
    taskCLP = ''

    if (taskName == self.fileNature.setup):
      taskCLP = self.basisInfo.scfBasisCode

    elif (taskName == self.fileNature.main):
      # The first CLP indicates which basis to use.
      taskCLP = self.basisInfo.scfBasisCode + " "

      # The next CLPs signal the edge that should be computed for. Note
      # that in the case of the ground state run we simply make an explicit
      # request for no excitation state.
      if (taskEdge == "gs_"):
        taskCLP += "0 0 "
      else:
        taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

      # The next CLP indicates if a dos calculation should also be done.
      taskCLP += str(self.jobInfo.dosRunSCF) + " "

      # The next CLP indicates if a bond calculation should also be done.
      taskCLP += str(self.jobInfo.bondRunSCF) + " "
    
    elif (taskName == self.fileNature.intg):
      # The first CLP indicates which basis to use.
      taskCLP = self.basisInfo.pscfBasisCode + " "

      # This CLP indicates whether the momentum matrix should be computed.
      if (((self.jobInfo.pacsRun == 1) and ('gs' in taskEdge)) or \
            (self.jobInfo.optcRun == 1) or (self.jobInfo.sigeRun)):
        taskCLP += "1"
      else:
        taskCLP += "0"

    elif (taskName == self.fileNature.band):
      # The first CLP indicates which basis to use.
      taskCLP = self.basisInfo.pscfBasisCode + " 0 "

      # The next CLPs signal the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

    elif (taskName == self.fileNature.sybd):
      # Do a sybd style calculation.
      taskCLP = self.basisInfo.pscfBasisCode + " 1 "

      #The next CLP signals the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

    elif (taskName == self.fileNature.dos):
      # The first CLP indicates what basis to use
      taskCLP = self.basisInfo.pscfBasisCode + " "

      #The next CLP signals the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

    elif (taskName == self.fileNature.bond):
      # The first CLP indicates what basis to use
      taskCLP = self.basisInfo.pscfBasisCode + " "

      #The next CLP signals the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

    elif (taskName == self.fileNature.optc):
      # The first CLP indicates what basis to use
      taskCLP = self.basisInfo.pscfBasisCode + " "

      # The second specifies that we are doing a normal optical properties
      # calculation (dielectric function).
      taskCLP += "0 "

      # The next CLPs signal the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

      # The final value indicates whether or not to computer the x,y,z
      # components in serial or all together (the default).
      taskCLP += self.control.serialxyz + " "

    elif (taskName == self.fileNature.pacs):
      # The first CLP indicates what basis to use
      taskCLP = self.basisInfo.pscfBasisCode + " "

      # The second specifies that we are doing a XANES/ELNES calculation.
      taskCLP += "1 "

      # The next CLPs signal the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "

      # The final value indicates whether or not to compute the x,y,z
      # components inserial or all together (the default).
      taskCLP += self.control.serialxyz + " "

    elif (taskName == self.fileNature.wave):
      # The first CLP indicates what basis to use
      taskCLP = self.basisInfo.pscfBasisCode + " "

      # The next CLPs signal the edge that should be computed for.
      taskCLP += self.control.QN_n + " " + self.control.QN_l + " "
   
    return taskCLP

  def executeProgram(self, taskName, taskCLP):
    # Need some local vars initialized
    executable = ''
    output = ''
    values = []

    # Get the name of the executable
    if (taskName == self.fileNature.setup):
      executable = self.exeNames.exe[0]
    elif (taskName == self.fileNature.main):
      executable = self.exeNames.exe[1]
    elif (taskName == self.fileNature.intg):
      executable = self.exeNames.exe[2]
    elif (taskName == self.fileNature.band):
      executable = self.exeNames.exe[3]
    elif (taskName == self.fileNature.sybd):
      executable = self.exeNames.exe[3]
      if (executable[0] == 'g'):
        executable = executable[1:]
    elif (taskName == self.fileNature.dos):
      executable = self.exeNames.exe[4]
    elif (taskName == self.fileNature.bond):
      executable = self.exeNames.exe[5]
    elif ((taskName == self.fileNature.dos) or 
            (taskName == self.fileNature.pacs) or 
            (taskName == self.fileNature.sige)):
      executable = self.exeNames.exe[6]
    elif (taskName == self.fileNature.wave):
      executable = self.exeNames.exe[7]

    # Call the executable
# So this will probably behave weird
##
##
## Make sure this works
    exestr = self.exeNames.exeMechanism + " " + self.env.obin + "/" + \
        executable + " " + taskCLP + " 2>&1 > ../runtime"
    retval = os.system(exestr)
#    self.runtime.write(retval) 

    # In certain cases a secondary job needs to be run immediately 
    # afterwards
    if (taskName == self.fileNature.sybd):
      systr = self.env.obin + '/makeSYBD -dat fort.5 -out fort.20 ' + \
          '-raw fort.31 -plot fort.41'
      os.system(systr)
      if (self.control.spinPol == 2): # spin polarized case
        systr = self.env.obin + '/makeSYBD -dat fort.5 -out fort.20 ' + \
            '-raw fort.32 -plot fort.42'
        os.system(systr)
    elif(taskName == self.fileNature.optc):
      # Perform the kramers-kronig conversion of eps2 to eps1 and compute
      # the energy loss function (all done in OLCAOkkc). We need to pass the
      # program the number of lines in the input file and which set of file
      # numbers to use. A 1 means use the spin up or default file numbers,
      # and a 2 means use the spin down file numbers.
      f50 = open('fort.50','r')
      linecount = str(len(f50.readlines()))
      
      # We always use the #1 file numbers.
      systr = self.env.obin + '/OLCAOkkc ' + linecount + '1 2>&1'
      os.system(systr)

      # We only use the #2 set of file numbers for the spin polarized
      # calculations.
      if (self.control.spinPol == 2):
        systr = self.env.obin + '/OLCAOkkc ' + linecount + '2 2>&1'
        os.system(systr)

    # Check for the existance fo the fort.2 file that signals completion of
    # the fortran executable without abortive error.
    if (not os.path.isfile('fort.2')):
      self.runtime.write('Fortran success file missing. Exiting Script.')
      sys.exit()
    else:
      os.system('rm -f fort.2')

  def copy(self, tocopy, dest):
    os.system("cp " + tocopy + " " + dest)

  def cpInter(self, tocopy, nature, tail, taskEdge, taskBasis):
    dest = taskEdge + nature + taskBasis + tail
    self.copy(tocopy, dest)

  def move(self, tomove, dest):
    os.system("mv " + tomove + " " + dest)

  def mvInter(self, tomove, nature, tail, taskEdge, taskBasis):
    dest = taskEdge + nature + taskBasis + tail
    self.move(tomove, dest)

  def mvProjHome(self, tomove, nature, tail, taskEdge, taskBasis):
    dest = self.ioNames.proj_home + '/' + taskEdge + nature + taskBasis + \
        tail
    self.move(tomove, dest)

  def manageOutput(self,taskName,taskBasis,taskEdge,taskKP,alt):
    tempBasis = "-temp"

    # Perform these specialized moved and copies for each possible task.
    if (taskName == self.fileNature.setup):
      start = self.fileNature.setup + tempBasis + self.fExtensions.hdf
      dest = self.fileNature.setup + taskBasis + alt + \
          self.fExtensions.hdf
      self.move(start,dest)

    elif (taskName == self.fileNature.main):
      tail = '.' + self.ioNames.iteration + alt + '.7'
      self.cpInter("fort.7", self.fileNature.main, tail,taskEdge,taskBasis)
      tail = alt + self.fExtensions.dat
      self.mvProjHome("fort.7", self.ioNames.iteration, tail,taskEdge,taskBasis)
      tail = '.' + self.ioNames.iteration + alt + '.14'
      self.cpInter("fort.14", self.fileNature.main, tail,taskEdge,taskBasis)
      tail = alt + self.fExtensions.dat
      self.mvProjHome("fort.14", self.ioNames.energy, tail,taskEdge,taskBasis)
      tail = '.SCFV' + alt + '.8'
      self.cpInter("fort.8", self.fileNature.main, tail,taskEdge,taskBasis)
      tail = alt + self.fExtensions.dat
      self.mvProjHome("fort.8", self.ioNames.initPot, tail,taskEdge,taskBasis)

      if (os.path.isfile("fort.1000")):
        tail = '.iterTDOS' + alt + '.1000'
        self.cpInter("fort.1000", self.fileNature.main, tail,taskEdge,taskBasis)
        tail = '.iterTDOS' + alt + self.fExtensions.plot
        self.mvProjHome('fort.1000', self.fExtensions.main, tail,taskEdge,taskBasis)
      
      start = self.fileNature.setup + tempBasis + self.fExtensions.hdf
      dest = self.fileNature.setup + taskBasis + alt + \
          self.fExtensions.hdf
      self.move(start,dest)
      start = self.fileNature.main + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.main + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      if (self.jobInfo.dosRunSCF == 1):
        if (self.control.spinPol == 2):
          ttail = self.auxSuffix.tdos + alt
          ptail = self.auxSuffix.pdos + alt
          ltail = self.auxSuffix.loci + alt
         
          tail = ttail + self.fExtensions.plot + '.60'
          self.cpInter("fort.60", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ttail + self.auxSuffix.up + self.fExtensions.plot
          self.mvProjHome("fort.60", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.fExtensions.raw + '.70'
          self.cpInter("fort.70", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.auxSuffix.up + self.fExtensions.raw
          self.mvProjHome("fort.70", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.fExtensions.plot + '.80'
          self.cpInter("fort.80", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.auxSuffix.up + self.fExtensions.plot
          self.mvProjHome("fort.80", self.fileNature.dos, tail,taskEdge,taskBasis)
          
          tail = ttail + self.fExtensions.plot + '.61'
          self.cpInter("fort.61", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ttail + self.auxSuffix.dn + self.fExtensions.plot
          self.mvProjHome("fort.61", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.fExtensions.raw + '.71'
          self.cpInter("fort.71", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.auxSuffix.dn + self.fExtensions.raw
          self.mvProjHome("fort.71", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.fExtensions.plot + '.81'
          self.cpInter("fort.81", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.auxSuffix.dn + self.fExtensions.plot
          self.mvProjHome("fort.81", self.fileNature.dos, tail,taskEdge,taskBasis)
        else:
          ttail = self.auxSuffix.tdos + alt
          ptail = self.auxSuffix.pdos + alt
          ltail = self.auxSuffix.loci + alt

          tail = ttail + self.fExtensions.plot + '.60'
          self.cpInter("fort.60", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ttail + self.fExtensions.plot
          self.mvProjHome("fort.60", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.fExtensions.raw + '.70'
          self.cpInter("fort.70", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ptail + self.fExtensions.raw
          self.mvProjHome("fort.70", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.fExtensions.plot + '.80'
          self.cpInter("fort.80", self.fileNature.dos, tail,taskEdge,taskBasis)
          tail = ltail + self.fExtensions.plot
          self.mvProjHome("fort.80", self.fileNature.dos, tail,taskEdge,taskBasis)

      if (self.jobInfo.bondRunSCF == 1):
        if (self.control.spinPol == 2):
          tail = alt + self.fExtensions.raw
          self.cpInter("fort.10",self.fileNature.bond, tail + ".10",taskEdge,taskBasis)
          tail = alt + self.auxSuffix.up + \
                self.fExtensions.plot
          self.mvProjHome("fort.10",self.fExtensions.bond, tail,taskEdge,taskBasis)

          tail = alt + self.fExtensions.raw
          self.cpInter("fort.11",self.fileNature.bond, tail + ".11",taskEdge,taskBasis)
          tail = alt + self.auxSuffix.dn + \
                self.fExtensions.plot
          self.mvProjHome("fort.11",self.fExtensions.bond, tail,taskEdge,taskBasis)
        else:
          tail = self.fExtensions.raw + '.10'
          self.cpInter("fort.10", self.fileNature.bond, tail,taskEdge,taskBasis)
          tail =  self.fExtensions.raw
          self.mvProjHome("fort.10", self.fileNature.bond, tail,taskEdge,taskBasis)

    elif (taskName == self.fileNature.intg):
      start = self.fileNature.intg + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.intg + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      dest = taskEdge + self.fileNature.intg + taskBasis + '.SCFV' + \
          alt + '.8'
      self.move('fort.8',dest)

    elif (taskName == self.fileNature.band):
      start = self.fileNature.intg + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.intg + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      start = self.fileNature.band + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.band + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

    elif (taskName == self.fileNature.sybd):
      start = self.fileNature.intg + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.intg + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      if (self.control.spinPol == 2):
        tail = alt + self.fExtensions.raw
        self.mvInter("fort.31", self.fileNature.sybd, tail + '.31',taskEdge,taskBasis) 
        self.mvInter("fort.32", self.fileNature.sybd, tail + '.32',taskEdge,taskBasis) 

        tail = alt + self.auxSuffix.up + \
            self.fExtensions.raw
        self.mvProjHome('fort.41', self.fileNature.sybd, tail,taskEdge,taskBasis)
        
        tail = alt + self.auxSuffix.dn + \
            self.fExtensions.raw
        self.mvProjHome('fort.42', self.fileNature.sybd, tail,taskEdge,taskBasis)

      else:
        tail = alt + self.fExtensions.raw
        self.mvInter("fort.31", self.fileNature.sybd, tail + '.31',taskEdge,taskBasis) 
        self.mvProjHome('fort.41', self.fileNature.sybd, tail,taskEdge,taskBasis)

    elif (taskName == self.fileNature.dos):
      start = self.fileNature.band + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.band + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      if (self.control.spinPol == 2):
        ttail = self.auxSuffix.tdos + alt + \
            self.auxSuffix.up + self.fExtensions.plot
        ptail = self.auxSuffix.pdos + alt + \
            self.auxSuffix.up + self.fExtensions.plot
        ltail = self.auxSuffix.loci + alt + \
            self.auxSuffix.up + self.fExtensions.plot

        tail = self.fExtensions.plot + alt + '.60'
        self.cpInter('fort.60', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.60', self.fileNature.dos, ttail,taskEdge,taskBasis)
        tail = alt + self.fExtensions.raw + '.70'
        self.cpInter('fort.70', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.70', self.fileNature.dos, ptail,taskEdge,taskBasis)
        tail = self.fExtensions.plot + alt + '.80'
        self.cpInter('fort.80', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.80', self.fileNature.dos, ltail,taskEdge,taskBasis)

        ttail = ttail.replace(self.auxSuffix.up,self.auxSuffix.dn)
        ptail = ttail.replace(self.auxSuffix.up,self.auxSuffix.dn)
        ltail = ttail.replace(self.auxSuffix.up,self.auxSuffix.dn)
        tail = self.fExtensions.plot + alt + '.61'
        self.cpInter('fort.61', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.61', self.fileNature.dos, ttail,taskEdge,taskBasis)
        tail = alt + self.fExtensions.raw + '.71'
        self.cpInter('fort.71', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.61', self.fileNature.dos, ptail,taskEdge,taskBasis)
        tail = self.fExtensions.plot + alt + '.81'
        self.cpInter('fort.81', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.81', self.fileNature.dos, ltail,taskEdge,taskBasis)

      else:
        ttail = self.auxSuffix.tdos + alt + \
            self.fExtensions.plot
        ptail = self.auxSuffix.pdos + alt + \
            self.fExtensions.plot
        ltail = self.auxSuffix.loci + alt + \
            self.fExtensions.plot

        tail = self.fExtensions.plot + alt + '.60'
        self.cpInter('fort.60', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.60', self.fileNature.dos, ttail,taskEdge,taskBasis)
        tail = alt + self.fExtensions.raw + '.70'
        self.cpInter('fort.70', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.70', self.fileNature.dos, ptail,taskEdge,taskBasis)
        tail = self.fExtensions.plot + alt + '.80'
        self.cpInter('fort.80', self.fileNature.dos, tail,taskEdge,taskBasis)
        self.mvProjHome('fort.80', self.fileNature.dos, ltail,taskEdge,taskBasis)

    elif (taskName == self.fileNature.bond):
      start = self.fileNature.band + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.band + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      if (self.control.spinPol == 2):
        tail = alt + self.fExtensions.raw
        tail2 = alt + self.auxSuffix.up + \
            self.fExtensions.raw
        self.cpInter("fort.10", self.fileNature.bond, tail + '.10',taskEdge,taskBasis)
        self.mvProjHome('fort.10',self.fileNature.bond,tail2,taskEdge,taskBasis)
        self.cpInter("fort.11", self.fileNature.bond, tail + '.11',taskEdge,taskBasis)
        tail2 = tail2.replace(self.auxSuffix.up,self.auxSuffix.dn)
        self.mvProjHome('fort.11',self.fileNature.bond,tail2,taskEdge,taskBasis)
      else:
        tail = alt + self.fExtensions.raw
        self.cpInter("fort.10", self.fileNature.bond, tail + '.10',taskEdge,taskBasis)
        self.mvProjHome('fort.10',self.fileNature.bond,tail,taskEdge,taskBasis)
    elif (taskName == self.fileNature.optc):
      start = self.fileNature.intg + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.intg + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)

      start = self.fileNature.band + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.band + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)
      if (self.control.spinPol == 2):
        altupplot = alt + self.auxSuffix.up + \
            self.fExtensions.plot
        altdnplot = alt + self.auxSuffix.dn + \
            self.fExtensions.plot

        tail = self.auxSuffix.cond + alt + \
          self.auxSuffix.up+ '.40'
        self.cpInter("fort.40", self.fileNature.optc,tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + alt + \
          self.auxSuffix.up+ '.50'
        self.cpInter("fort.50", self.fileNature.optc,tail,taskEdge,taskBasis)

        tail = self.auxSuffix.cond + altupplot
        self.mvProjHome("fort.40", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + altupplot
        self.mvProjHome("fort.50", self.fileNature.optc, tail,taskEdge,taskBasis)
        self.mvProjHome("fort.100", self.fileNature.optc, altupplot)
        tail = self.auxSuffix.eps1 + altupplot
        self.mvProjHome("fort.110", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.elf + altupplot
        self.mvProjHome("fort.120", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.ref + altupplot
        self.mvProjHome("fort.130", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps1i + altupplot
        self.mvProjHome("fort.140", self.fileNature.optc, tail,taskEdge,taskBasis)
        
        tail = self.auxSuffix.cond + alt + \
          self.auxSuffix.dn + '.41'
        self.cpInter("fort.41", self.fileNature.optc,tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + alt + \
          self.auxSuffix.dn + '.51'
        self.cpInter("fort.51", self.fileNature.optc,tail,taskEdge,taskBasis)

        tail = self.auxSuffix.cond + altdnplot
        self.mvProjHome("fort.41", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + altdnplot
        self.mvProjHome("fort.51", self.fileNature.optc, tail,taskEdge,taskBasis)
        self.mvProjHome("fort.101", self.fileNature.optc, altupplot)
        tail = self.auxSuffix.eps1 + altdnplot
        self.mvProjHome("fort.111", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.elf + altdnplot
        self.mvProjHome("fort.121", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.ref + altdnplot
        self.mvProjHome("fort.131", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps1i + altdnplot
        self.mvProjHome("fort.141", self.fileNature.optc, tail,taskEdge,taskBasis)

      else:
        altplot = alt + self.fExtensions.plot

        tail = self.auxSuffix.cond + alt + '.40'
        self.cpInter("fort.40", self.fileNature.optc,tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + alt + '.50'
        self.cpInter("fort.50", self.fileNature.optc,tail,taskEdge,taskBasis)

        tail = self.auxSuffix.cond + altplot
        self.mvProjHome("fort.40", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps2 + altplot
        self.mvProjHome("fort.50", self.fileNature.optc, tail,taskEdge,taskBasis)
        self.mvProjHome("fort.100", self.fileNature.optc, altupplot)
        tail = self.auxSuffix.eps1 + altplot
        self.mvProjHome("fort.110", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.elf + altplot
        self.mvProjHome("fort.120", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.ref + altplot
        self.mvProjHome("fort.130", self.fileNature.optc, tail,taskEdge,taskBasis)
        tail = self.auxSuffix.eps1i + altplot
        self.mvProjHome("fort.140", self.fileNature.optc, tail,taskEdge,taskBasis)

    elif (taskName == self.fileNature.wave):
      start = self.fileNature.band + tempBasis + self.fExtensions.hdf
      dest = taskEdge + self.fileNature.band + taskBasis + \
        alt + self.fExtensions.hdf
      self.move(start,dest)
      
      tail = '.SCFV' + alt + '.8'
      self.mvInter("fort.8", self.fileNature.intg, tail,taskEdge,taskBasis)

      # The profile files are always created.
      altprof = alt + self.auxSuffix.profile
      tail = altprof + "-a" + self.fExtensions.dat
      self.mvProjHome("fort.30",self.fileNature.wave, tail,taskEdge,taskBasis)
      tail = altprof + "-b" + self.fExtensions.dat
      self.mvProjHome("fort.31",self.fileNature.wave, tail,taskEdge,taskBasis)
      tail = altprof + "-c" + self.fExtensions.dat
      self.mvProjHome("fort.32",self.fileNature.wave, tail,taskEdge,taskBasis)

      if (os.path.isfile("fort.56")): # Only move openDx files if exists
          dest = self.ioNames.proj_home + ',' + self.ioNames.atomPos + \
              self.fExtensions.dx
          self.move("fort.56",dest)
          dest = self.ioNames.proj_home + ',' + self.ioNames.lattic + \
              self.fExtensions.dx
          self.move("fort.56",dest)

          if (self.control.spinPol == 2):
            tail = alt + self.auxSuffix.valeRho + \
                self.auxSuffix.upPdn + self.fExtensions.dx
            self.mvProjHome("fort.58", self.fileNature.wave, tail,taskEdge,taskBasis)
            tail = alt + self.auxSuffix.valeRho + \
                self.auxSuffix.upMdn + self.fExtensions.dx
            self.mvProjHome("fort.59", self.fileNature.wave, tail,taskEdge,taskBasis)
            tail = alt + self.auxSuffix.valeRho + self.auxSuffix.upPdn + \
                self.auxSuffix.minusNeut + self.fExtensions.dx
            self.mvProjHome("fort.60", self.fileNature.wave, tail,taskEdge,taskBasis)
            tail = alt + self.auxSuffix.pot + \
                self.auxSuffix.up + self.auxSuffix.dx 
            self.mvProjHome("fort.61", self.fileNature.wave, tail,taskEdge,taskBasis)
            tail = alt + self.auxSuffix.pot + \
                self.auxSuffix.dn + self.auxSuffix.dx 
            self.mvProjHome("fort.62", self.fileNature.wave, tail,taskEdge,taskBasis)
          else:
            tail = alt + self.auxSuffix.valeRho + self.fExtensions.dx
            self.mvProjHome("fort.58", self.fileNature.wave, tail,taskEdge,taskBasis)
            tail = alt + self.auxSuffix.valeRho + \
                self.auxSuffix.minusNeut + self.fExtensions.dx
            self.mvProjHome("fort.59", self.fileNature.wave, tail),taskEdge,taskBasis
            tail = alt + self.auxSuffix.pot +  self.fExtensions.dx
            self.mvProjHome("fort.60", self.fileNature.wave, tail,taskEdge,taskBasis)

    # Perform these moves for all tasks.
    self.mvInter("fort.5", taskName, self.fExtensions.dat + ".5",taskEdge,taskBasis)

    self.move("fort.15", taskKP+".15")

    dest = self.ioNames.structure + self.fExtensions.dat + ".4"
    self.move("fort.4",dest)

    dest = alt + self.fExtensions.out + ".20"
    self.cpInter("fort.20",taskName,dest,taskEdge,taskBasis)

    dest = alt + self.fExtensions.out
    self.mvProjHome("fort.20", taskName, dest,taskEdge,taskBasis)

  def updatePACS(self, finMainEdge, updateSCFBasis):
    # This subroutine is run after both main calculations for the XAS
    # total calculation are complete. Then, we can determine the
    # approximate location of the edge onset based on the total energy
    # difference between the two states (initial/final).

    os.chdir(self.ioNames.proj_home)

    # Find the total energy of the intitial state.
    tmpfn = "gs_" + self.fileNature.main + updateSCFBasis + \
        self.fExtensions.out
    initmain = open(tmpfn,'r')

    for line in initmain.readlines():
      if ('TOTAL ENERGY' in line):
        vals = (line.strip()).split()
        initEnergy = vals[-1]

    # Find the total energy of the final state.
    tmpfn = finMainEdge + self.fileNature.main + updateSCFBasis + \
        self.fExtensions.out
    finmain = open(tmpfn,'r')

    for line in finmain.readlines():
      if ('TOTAL ENERGY' in line):
        vals = (line.strip()).split()
        finEnergy = vals[-1]
    finmain.close()

    # Use the collected info to determine the total energy difference
    # (TEdiff).
    TEdiff = abs(finEnergy-initEnergy) * self.consts.eVhartree

    # It is assumed that the pacs input file already exists at this point.
    # We are just going to replace the "UNKOWN" values.
    pacsold = open(self.ioNames.olcao+self.fExtensions.dat,'r')
    pacsnew = open('pacs_temp', 'w')

    numlines = len(pacsold.readlines())
    
    found = 0
    for x in xrange(numlines):
      if ('PACS_INPUT_DATA' in polines[x]):
        line = pacsold.readline()
        pacsnew.write(line)
        line = pacsold.readline()
        pacsnew.write(line)
        line = pacsold.readline()
        pacsnew.write(line)
        line = pacsold.readline()
        pacsnew.write(line)
        line = pacsold.readline()
        pacsnew.write(line)

        vals = (line.strip()).split()
        numCoreOrbitals = vals[0]

        for y in range(1,numCoreOrbitals+1):
          line = pacsold.readline()
          vals = (line.strip()).split()
          if ((vals[4] == -1) and (vals[0] == self.control.QN_n) and \
              (vals[1] == self.control.QN_l)):
            
            wrtstr = vals[0] + " " + vals[1] + " " + vals[2] + " " + \
                     vals[3] + " " + TEdiff + " " + vals[5] + " " + \
                     " # QN_n QN_l Init1 Init2 TEDiff\n"
            pacsnew.write(wrtstr)
            found = 1
          else:
            line = pacsold.readline()
            pacsnew.write(line)
      else:
        line = pacsold.readline()
        pacsnew.write(line)

    pacsold.close()
    pacsnew.close()

    # Update the version of the pacs input file in the inputs directory
    os.system("mv pacs_temp " + "olcao.dat");

    # Only print a notice if we changed anything.
    if (found == 1):
      self.runtime.write("olcao.dat Updated.\n")
