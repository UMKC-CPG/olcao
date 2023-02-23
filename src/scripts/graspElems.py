#!/usr/bin/python

# PROGRAM: makeGraspScripts
# AUTHOR:  Patrick Ryan Thomas  (prtnpb@mail.umkc.edu)
#
# The following scrip will take in command line inputs
#   and generate the necessary bash scripts to run
#   Grasp2K calculation.
#
#
# Example Inputs:
#   1) All elements all bases:
#      ./makeGraspScripts total
#
#   2) All elements specific basis:
#      ./makeGraspScripts total MB
#      ./makeGraspScripts total FB
#      ./makeGraspScripts total EB
#
#   3) Specific atom specific basis
#      ./makeGraspScripts F MB H EB
#

import numpy as np
import time
import os
import sys
import copy


def getAtomNum(symb):
  elems={'H': 1 , 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103}
  return elems[symb]



def getAtomList():

  l=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

  return l


def parseInputs():
  tmparr=[]
  basis=['MB','FB','EB']
  elem_list=getAtomList()
  commandout=[]
  global numCore
  numCore=0


  if "-nc" in sys.argv:
    numCore=int(sys.argv[sys.argv.index("-nc")+1])
  else:
    numCore=1

  for i in range(1,len(sys.argv)):
    commandout=[]
    if(sys.argv[i]=='total'):
      try:
        sys.argv[i+1]
      except IndexError:
        for j in range(0,len(elem_list)):
          commandout=[elem_list[j],'MB',elem_list[j],'FB',elem_list[j],'EB']
          tmparr.extend(commandout)
        return tmparr
      else:
        if sys.argv[i+1]=="-nc":
          for j in range(0,len(elem_list)):
            commandout=[elem_list[j],'MB',elem_list[j],'FB',elem_list[j],'EB']
            tmparr.extend(commandout)
          return tmparr
        else:
          for j in range(0,len(elem_list)):
            commandout=[elem_list[j],sys.argv[i+1]]
            tmparr.extend(commandout)
          return tmparr
    elif(sys.argv[i] in elem_list):
      try:
        sys.argv[i+1]
      except IndexError:
        print('Incorrect Submission')
        break
      else:
        if(sys.argv[i+1]=='all'):
          commandout=[sys.argv[i],'MB',sys.argv[i],'FB',sys.argv[i],'EB']
          tmparr.extend(commandout)
        else:
          commandout=[sys.argv[i],sys.argv[i+1]]
          tmparr.extend(commandout)

  return tmparr

def getFileContents(filename):
  f=open(filename,'r')
  content=f.readlines()
  f.close()
  return content


class elemData:

  def __init__(self):
    self.atomSyms=[]
    self.atomNums=[]
    self.coreAtoms=[]
    self.mbValAtoms=[]
    self.fbValAtoms=[]
    self.ebValAtoms=[]
    self.occAtoms=[]
    self.mbNumSubs=[]
    self.fbNumSubs=[]
    self.ebNumSubs=[]
    self.mbERWF=[]
    self.fbERWF=[]
    self.ebERWF=[]
    self.isotopes=[]
    self.atomicMass=[]
    self.spin=[]
    self.magMoment=[]
    self.quadMoment=[]



  def getElemData(self):
    #Open the isotope data file.
    atomInfo=getFileContents("isotopedata.dat")

    n=0
    for i in range(54,len(atomInfo),18):
      line=(atomInfo[i].strip()).split()
      if line[0]=="ATOM_SYM:":
        n+=1

        self.atomSyms.append(line[1])
        self.atomNums.append((atomInfo[i+1].split())[1]) 
        self.coreAtoms.append(atomInfo[i+3].strip())
        self.mbValAtoms.append(atomInfo[i+5].strip())
        self.fbValAtoms.append(atomInfo[i+6].strip())
        self.ebValAtoms.append(atomInfo[i+7].strip())
 
        self.occAtoms.append(atomInfo[i+9].strip())
 
        line2=(atomInfo[i+11].strip()).split()
        self.mbNumSubs.append(line2[0])
        self.fbNumSubs.append(line2[1])
        self.ebNumSubs.append(line2[2])
        
        line2=(atomInfo[i+13].strip()).split()
        self.mbERWF.append(line2[0])
        self.fbERWF.append(line2[1])
        self.ebERWF.append(line2[2])
        
        line2=(atomInfo[i+14].strip()).split()
        self.isotopes.append(line2[0])
        self.atomicMass.append(line2[2])
        if line2[3]=="#":
          self.spin.append('0')
        else:
          self.spin.append(line2[3])
        if line2[4]=="#":
          self.magMoment.append('1')
        else:
          self.magMoment.append(line2[4])
        if line2[5]=="#":
          self.quadMoment.append('1')
        else:
          self.quadMoment.append(line2[5])

      


def makeExecs(commList):
  elem_list=getAtomList()
  atomDat=elemData()
  atomDat.getElemData()
  execLocs=[]
  batchHead=[]
  batchFileNum=1
  elemCount=1


  if os.path.exists("SBATCHhead"):
    f=open("SBATCHhead",'r')
    batchHead=f.readlines()
    f.close()
    jobFile=open("JobSBATCH-"+str(batchFileNum),'w')
    jobFile2=open("JobRSCF-"+str(batchFileNum),'w')
    for i in range(len(batchHead)):
      if i==1:
        jobFile.write(batchHead[i].strip()+str(batchFileNum)+"\n")
        jobFile2.write(batchHead[i].strip()+str(batchFileNum)+"RSCF\n")
      elif i==2:
        jobFile.write(batchHead[i].strip()+str(batchFileNum)+".o%j\n")
        jobFile2.write(batchHead[i].strip()+str(batchFileNum)+"RSCF.o%j\n")
      else:
        jobFile.write(batchHead[i]) 
        jobFile2.write(batchHead[i]) 
    jobFile.write("\n")
    jobFile2.write("\n")


  for i in range(len(commList)):
    if commList[i] in elem_list:
        
      if os.path.isdir(commList[i+1]):
        os.mkdir(commList[i+1]+"/"+commList[i])
        os.chdir(commList[i+1]+"/"+commList[i])

      else:
        os.mkdir(commList[i+1])
        os.mkdir(commList[i+1]+"/"+commList[i])
        os.chdir(commList[i+1]+"/"+commList[i])

      if os.path.exists("../../SBATCHhead"):
        if elemCount%numCore==0 and elemCount !=1:
          jobFile.write("\nwait\n")
          jobFile2.write("\nwait\n")
          batchFileNum+=1
          jobFile.close()
          jobFile2.close()
          jobFile=open("../../JobSBATCH-"+str(batchFileNum),'w')
          jobFile2=open("../../JobRSCF-"+str(batchFileNum),'w')
          for k in range(len(batchHead)):
            if k==1:
              jobFile.write(batchHead[k].strip()+str(batchFileNum)+"\n")
              jobFile2.write(batchHead[k].strip()+str(batchFileNum)+"RSCF\n")
            elif k==2:
              jobFile.write(batchHead[k].strip()+str(batchFileNum)+".o%j\n")
              jobFile2.write(batchHead[k].strip()+str(batchFileNum)+"RSCF.o%j\n")
            else:
              jobFile.write(batchHead[k]) 
              jobFile2.write(batchHead[k]) 
          
        jobFile.write("cd "+commList[i+1]+"/"+commList[i]+"\n")
        jobFile.write("./run"+commList[i]+" > "+commList[i]+"out &\n")
        jobFile.write("cd ../..\n")

        jobFile2.write("cd "+commList[i+1]+"/"+commList[i]+"\n")
        jobFile2.write("./runRSCF > "+commList[i]+"RSCFout &\n")
        jobFile2.write("cd ../..\n")

      atNum=getAtomNum(commList[i])
     
      #ISOTOPE
      g=open("run"+commList[i],'w')
      g.write(batchHead[0])
      g.write("set -x\n")
      g.write("\n$GRASP/bin/iso <<S1\n")
      g.write(str(atNum)+"\n")
      g.write((atomDat.isotopes[atNum-1]).replace(commList[i],''))
      g.write('\nn\n')
      g.write(atomDat.atomicMass[atNum-1]+'\n')
      g.write(atomDat.spin[atNum-1]+'\n')
      g.write(atomDat.magMoment[atNum-1]+'\n')
      g.write(atomDat.quadMoment[atNum-1]+'\n')
      g.write('S1\n\n')

      #CSL
      g.write("$GRASP/bin/csl <<S1\n")
      g.write("y\n")
      g.write(atomDat.coreAtoms[atNum-1]+"\n")
      if commList[i+1]=="MB":
        g.write(atomDat.mbValAtoms[atNum-1]+"\n")
      elif commList[i+1]=="FB":
        g.write(atomDat.fbValAtoms[atNum-1]+"\n")
      elif commList[i+1]=="EB":
        g.write(atomDat.ebValAtoms[atNum-1]+"\n")
      g.write("n\nn\n")
      g.write(atomDat.occAtoms[atNum-1]+"\n")
      g.write("\n\n")
      if commList[i+1]=="MB":
        g.write(atomDat.mbNumSubs[atNum-1]+"\n")
      elif commList[i+1]=="FB":
        g.write(atomDat.fbNumSubs[atNum-1]+"\n")
      elif commList[i+1]=="EB":
        g.write(atomDat.ebNumSubs[atNum-1]+"\n")
      g.write("S1\n\n")

      g.write("mv rcsl.out rcsl.inp\n\n")

      #JSPLIT
      g.write("$GRASP/bin/jsplit <<S1\n")
      g.write("y\n")
      g.write("S1\n\n")

      g.write('rm rcsl.inp\n')
      g.write('mv rcsl.out rcsl.inp\n\n')

      #MCP3
      g.write("$GRASP/bin/mcp3 <<S1\n")
      g.write("y\n")
      g.write("S1\n\n")

      #ERWF
      g.write("$GRASP/bin/erwf <<S1\n")
      g.write("y\n") 
      if commList[i+1]=="MB":
        g.write(atomDat.mbERWF[atNum-1]+"\n")
      elif commList[i+1]=="FB":
        g.write(atomDat.fbERWF[atNum-1]+"\n")
      elif commList[i+1]=="EB":
        g.write(atomDat.ebERWF[atNum-1]+"\n")
      g.write("*\n")
      g.write("S1\n")      

      g.close()
      comm="chmod +x run"+commList[i]
      os.system(comm)

      g=open("runRSCF",'w')

      g.write(batchHead[0])
      g.write("set -x\n")

      g.write("$GRASP/bin/rscf2 <<S1\n")
      g.write("y\ny\n#BLOCKS AS 1's\n")
      g.write("*\n\n100000\nS1\n\n")

     
      g.write("mv rwfn.inp rwfn-orig.inp\n")
      g.write("mv rwfn.out rwfn.inp\n\n")

      g.write("$GRASP/bin/readrwf <<S1\n")
      g.write("1\n")
      g.write("rwfn.inp\n")
      g.write("rwfn.out\n")
      g.write("S1\n")
      
      g.close()
      comm="chmod +x runRSCF"
      os.system(comm)
      os.chdir('../..')
      elemCount+=1

      
  jobFile.write("\nwait\n")
  jobFile2.write("\nwait\n")

  jobFile.close()
  jobFile2.close()




  

elemComs=parseInputs()
makeExecs(elemComs)


