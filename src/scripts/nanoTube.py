#!/usr/bin/env python

'''
Program: nanotube.py
Author: Patrick Ryan Thomas
        prtnpb@mail.umkc.edu

Purpose: The purpose of this program is to take the smallest
         translationally invariant portion of a structure and
         through a series of translations and rotations, make a
         nanotube structure that meets periodic boundary conditions
         in the z-direction

Inputs: The program takes in a list of element names corresponding
        to a list of a list of atomic position. The program also
        requires the initial cell dimensions.

Calculations: The program will then calculate all of the plausible 
              radii and rotation angles by rules of polygon geometry
              and trigonometry.

              Given the y-width of the cell as the length of one
              side of a polygon.  The total sum of the angles of the
              polygon is pi*(n-2) where n is the number of sides.

              Dividing this angle by the number of angles which is
              equal to n, we get each internal angle around the
              perimeter of the polygon. Further, dividing by two
              creates a right triangle. Therefore,

                 ( pi     )   r 
              tan| --(n-2)| = --
                 ( 2n     }   l/2

              For the angle of rotation, since we have a right
              triangle:

              angle_of_rotation=2(pi/2 - (pi/2n)*(n-2))

              Traditional rotation matrices are then used for
              rotations about the principle axis both in setting up
              the initial structures and in the duplication and 
              rotation of cells.

              Arbitrarily number of sides can be calculated for all
              the orientations from 15 to 80 sides. It then becomes
              the job of the investigator to choose the correct
              structure.

Outputs: Skeleton file structure for use in OLCAO.
'''

import os
import sys
import numpy as np
import time
import copy
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D


start=time.time()

'''
  Rotation Descriptions:

  Rotation 1:
    Since all of the positions are in the first octant, i.e.
    (0<x,0<y,0<z). We simply take all positions in the provided
    translationally invariant structure and subtract y-width/2
    from all of the y values.

  Rotation 2:
   For this rotation, we do the following:
     1) Rotate by pi/2 about the y-axis
     2) Rotate by pi about the z-axis
     3) Add initial y-width/2 to all atomic positions

  Rotation 3:
    For this rotation, we do the following:
      1) Rotate about the x-axis by pi/2
      2) Shift by the initial z-width/2

  Rotation 4:
    For this rotation, we do the following:
      1)
'''




class initStructure:
  '''
  CLASS VARIABLES
    initPositions: list of positions [[x_1,y_1,z_1],[x_2,y_2,z_2]...]
    corr_atoms: list of element names associated with initPositions
    width: vector containing the initial cell vectors.
    orientpos1: similar to initPositions, but with Rotation 1 
                (see above) applied
    width1: In this case, the x,y,z widths do not change because we aren't
            rotating the cell so this is the same as the initial.
 
    orientpos2: similar to initPositions, but with Rotation 2 applied.
    width2: x,y, and z widths after Rotation 2 is applied
    orientpos3: similar to initPositions, but with Rotation 3 applied.
    width3: x,y, and z widths after Rotation 3 is applied

  CLASS FUNCTIONS
    Initalization is accomplished by providing the inital cell widths.

    addAtom: This allows for adding in atoms and positions individually.
             Takes string of element type and list of x,y,z position.

    addEntireStructure: Allows for passing of list of element names and
                        list of list of atomic positions. 
                        See initPositions.
                        Takes list of element names and list of list
                        of atomic positions.

    setinitPos1: Performs the shifting described Rotation 1 on atoms in the
                 initial structure
      
    setinitPos2: Peforms the rotations and shifting described in
                 Rotation 2.

    setinitPos3: Performs the rotations and shifting described in
                 Rotation 3.

    setAllThreePos: Peforms all three functions of setinitPos1, 2, and 3.

    doRot1/2/3: 1) The loop determines the radii and the angle of rotation.
                2) Shift all x positions by the radii amount.
                3) Peform number sides -1 rotations about the z-axis.
    
    doAllRots: Will make all three functions

  '''
  
  def __init__(self,cellx,celly,cellz):
    self.initPositions=[]           
    self.corr_atoms=[]             
    self.possibleRadii=[]
    self.width=[cellx,celly,cellz]
    self.orientpos1=[]
    self.width1=[]
    self.orientpos2=[]
    self.width2=[0,0,0]
    self.orientpos3=[]
    self.width3=[0,0,0]
    self.orientpos4=[]
    self.width4=[0.0,0.0,0.0]

  def addAtom(self,atom,pos):
    self.initPositions.append(pos[:])
    self.corr_atoms.append(atom)


  def addEntireStructure(self,atomlist,poslist):
    self.initPositions=list(copy.deepcopy(poslist))
    self.corr_atoms=list(copy.deepcopy(atomlist))


  def setinitPos1(self):
    '''
    This particular function simply shifts the structure cell by
    half of the current y-value such that it is centered on the x-axis.
    '''
    self.orientpos1=list(copy.deepcopy(self.initPositions))
    for i in range(len(self.orientpos1)):
      self.orientpos1[i][1]=self.orientpos1[i][1]-(self.width[1]/2.0)
    self.width1=list(copy.deepcopy(self.width))


  def setinitPos2(self):
    '''
    This function rotates about the y-axis and then the z-axis. Then shifts
    such that the structure is centered on the x-axis.
    '''
    self.orientpos2=list(copy.deepcopy(self.initPositions))
    Rot=np.zeros((3,3),dtype='d')
    Rot[0][2]=1.0
    Rot[1][1]=-1.0
    Rot[2][0]=1.0
    for i in range(len(self.orientpos2)):
      self.orientpos2[i]=(np.dot(Rot,np.array(self.orientpos2[i],dtype='d'))).tolist()
      self.orientpos2[i][1]=self.orientpos2[i][1]+(self.width[1]/2.0)
    self.width2[0]=self.width[2]
    self.width2[1]=self.width[1]
    self.width2[2]=self.width[0]

  def setinitPos3(self):
    '''
    This function rotates about the x-axis and then shifts such that
    the structure is centered on the x-axis.
    '''
    self.orientpos3=list(copy.deepcopy(self.initPositions))
    Rot=np.zeros((3,3),dtype='d')
    Rot[0][0]=1.0
    Rot[1][2]=-1.0
    Rot[2][1]=1.0
    for i in range(len(self.orientpos3)):
      self.orientpos3[i]=(np.dot(Rot,np.array(self.orientpos3[i],dtype='d'))).tolist()

      self.orientpos3[i][1]=self.orientpos3[i][1]+(self.width[2]/2.0)
    self.width3[0]=self.width[0]
    self.width3[1]=self.width[2]
    self.width3[2]=self.width[1]


  def setinitPos4(self):
    self.orientpos4=list(copy.deepcopy(self.initPositions))

    Rot=np.zeros((3,3),dtype='d')
    Rot[0][1]=1.0
    Rot[1][0]=-1.0
    Rot[2][2]=1.0

    for i in range(len(self.orientpos4)):
      self.orientpos4[i]=(np.dot(Rot,np.array(self.orientpos4[i],dtype='d'))).tolist()
      self.orientpos4[i][1]+=self.width[0]/2.0
    self.width4[0]=self.width[1]
    self.width4[1]=self.width[0]
    self.width4[2]=self.width[2]



  def setAllThreePos(self):
    self.setinitPos1()
    self.setinitPos2()
    self.setinitPos3()
    self.setinitPos4()



  def doRot1(self):
    tuberad=0.0
    ang=0.0
    rotang=0.0
    RotMat=np.zeros((3,3),dtype='d')

    for i in range(15,81):
      # Get rotation angle and tube radius
      ang=(np.pi*(i-2.0))/(2.0*(i))
      rotang=2*((np.pi/2)-ang)
      tuberad=(self.width1[1]/2.0)*np.tan(ang)

      # Make a temporary copy of the atomic positions
      temp=np.array(list(copy.deepcopy(self.orientpos1)),dtype='d')

      # Shift the cell x-values out such that the inner cell face is
      #  at the specified radius
      for j in range(len(temp)):
        temp[j][0]=temp[j][0]+tuberad


      # Make the tube coordinates the new coordinates
      tubecoords=temp[:][:]
      elemlist=list(copy.deepcopy(self.corr_atoms)) 

      # Do n-1 duplication / rotations of the structure until the ring is formed.
      for j in range(1,i):
        RotMat[0][0]=np.cos(rotang*j)
        RotMat[0][1]=-1.0*np.sin(rotang*j)
        RotMat[1][0]=np.sin(rotang*j)
        RotMat[1][1]=np.cos(rotang*j)
        RotMat[2][2]=1.0
        elemlist.extend(list(copy.deepcopy(self.corr_atoms)))


        for k in range(len(temp)):
          tubecoords=np.vstack([tubecoords,np.dot(RotMat,np.array(temp[k],dtype='d'))])
            
      #fig=plt.figure()
      #ax = fig.add_subplot(111, projection='3d')
      #plt.ion()
      #plt.show()
      #for j in range(len(tubecoords)):
      #  ax.scatter(tubecoords[j][0],tubecoords[j][1],tubecoords[j][2])
      #  plt.draw()

      # STARTING OUTPUT HERE

      xmax=0.0
      ymax=0.0

      #for j in range(len(tubecoords)):
      #  if abs(tubecoords[j][0])>xmax:
      #    xmax=abs(tubecoords[j][0])
      #  if abs(tubecoords[j][1])>ymax:
      #    ymax=abs(tubecoords[j][1])

      newcellwidth=2.0*tuberad+2.0*self.width1[0]+30.0

      for j in range(len(tubecoords)):
        tubecoords[j][0]+=newcellwidth/2.0
        tubecoords[j][1]+=newcellwidth/2.0
        
      # Uncomment this for another z layer.
      #tubetemp=np.array(list(copy.deepcopy(tubecoords)),dtype='d')
      #for j in range(len(tubetemp)):
      #  tubetemp[j][2]+=self.width1[2] 
      #tubecoords=np.vstack([tubecoords,tubetemp])
      #elemlist.extend(list(copy.deepcopy(elemlist)))
            

      fileout=open('NT_TiO2_'+str(int(np.ceil((2*np.pi)/rotang)))+'_123_OC0.dat','w')
      fileout.write('title\n')
      fileout.write('structure with radius '+str(tuberad)+'\n')
      fileout.write('end\n')
      fileout.write('cell\n')

      # Add a *2 on the z term if doing another layer.
      fileout.write(str(newcellwidth)+'  '+str(newcellwidth)+'  '+str(self.width1[2])+' 90.000 90.000 90.000\n')
      fileout.write('cartesian '+str(len(elemlist))+'\n')
      for j in range(len(tubecoords)):
        fileout.write(elemlist[j]+'\t'+str(tubecoords[j][0])+'\t'+str(tubecoords[j][1])+'\t'+str(tubecoords[j][2])+'\n')
      fileout.write('space 1_a\n')
      fileout.write('supercell 1 1 1\n')
      fileout.write('full\n')
      fileout.close()

            
  def doRot2(self):
    tuberad=0.0
    ang=0.0
    rotang=0.0
    RotMat=np.zeros((3,3),dtype='d')

    for i in range(15,81):
      ang=(np.pi*(i-2.0))/(2.0*(i))
      rotang=2*((np.pi/2)-ang)
      tuberad=(self.width2[1]/2.0)*np.tan(ang)

      temp=np.array(list(copy.deepcopy(self.orientpos2)),dtype='d')

      for j in range(len(temp)):
        temp[j][0]=temp[j][0]+tuberad

      tubecoords=temp[:][:]
      elemlist=list(copy.deepcopy(self.corr_atoms)) 
      for j in range(1,i):
        RotMat[0][0]=np.cos(rotang*j)
        RotMat[0][1]=-1.0*np.sin(rotang*j)
        RotMat[1][0]=np.sin(rotang*j)
        RotMat[1][1]=np.cos(rotang*j)
        RotMat[2][2]=1.0
        elemlist.extend(list(copy.deepcopy(self.corr_atoms)))


        for k in range(len(temp)):
          tubecoords=np.vstack([tubecoords,np.dot(RotMat,np.array(temp[k],dtype='d'))])
            
      #fig=plt.figure()
      #ax = fig.add_subplot(111, projection='3d')
      #plt.ion()
      #plt.show()
      #for j in range(len(tubecoords)):
      #  ax.scatter(tubecoords[j][0],tubecoords[j][1],tubecoords[j][2])
      #  plt.draw()
      xmax=0.0
      ymax=0.0

      #for j in range(len(tubecoords)):
      #  if abs(tubecoords[j][0])>xmax:
      #    xmax=abs(tubecoords[j][0])
      #  if abs(tubecoords[j][1])>ymax:
      #    ymax=abs(tubecoords[j][1])

      newcellwidth=2.0*tuberad+2.0*self.width2[0]+30.0

      for j in range(len(tubecoords)):
        tubecoords[j][0]+=newcellwidth/2.0
        tubecoords[j][1]+=newcellwidth/2.0
        

      #tubetemp=np.array(list(copy.deepcopy(tubecoords)),dtype='d')
      #for j in range(len(tubetemp)):
      #  tubetemp[j][2]+=self.width2[2] 
      #tubecoords=np.vstack([tubecoords,tubetemp])
      #elemlist.extend(list(copy.deepcopy(elemlist)))
            

      fileout=open('NT_TiO2_'+str(int(np.ceil((2*np.pi)/rotang)))+'_321_OC0.dat','w')
      fileout.write('title\n')
      fileout.write('structure with radius '+str(tuberad)+'\n')
      fileout.write('end\n')
      fileout.write('cell\n')
      fileout.write(str(newcellwidth)+'  '+str(newcellwidth)+'  '+str(self.width2[2])+' 90.000 90.000 90.000\n')
      fileout.write('cartesian '+str(len(elemlist))+'\n')
      for j in range(len(tubecoords)):
        fileout.write(elemlist[j]+'\t'+str(tubecoords[j][0])+'\t'+str(tubecoords[j][1])+'\t'+str(tubecoords[j][2])+'\n')
      fileout.write('space 1_a\n')
      fileout.write('supercell 1 1 1\n')
      fileout.write('full\n')
      fileout.close()


  def doRot3(self):
    tuberad=0.0
    ang=0.0
    rotang=0.0
    RotMat=np.zeros((3,3),dtype='d')

    for i in range(15,81):
      ang=(np.pi*(i-2.0))/(2.0*(i))
      rotang=2*((np.pi/2)-ang)
      tuberad=(self.width3[1]/2.0)*np.tan(ang)

      temp=np.array(list(copy.deepcopy(self.orientpos3)),dtype='d')

      for j in range(len(temp)):
        temp[j][0]=temp[j][0]+tuberad

      tubecoords=temp[:][:]
      elemlist=list(copy.deepcopy(self.corr_atoms)) 
      for j in range(1,i):
        RotMat[0][0]=np.cos(rotang*j)
        RotMat[0][1]=-1.0*np.sin(rotang*j)
        RotMat[1][0]=np.sin(rotang*j)
        RotMat[1][1]=np.cos(rotang*j)
        RotMat[2][2]=1.0
        elemlist.extend(list(copy.deepcopy(self.corr_atoms)))


        for k in range(len(temp)):
          tubecoords=np.vstack([tubecoords,np.dot(RotMat,np.array(temp[k],dtype='d'))])
            
      #fig=plt.figure()
      #ax = fig.add_subplot(111, projection='3d')
      #plt.ion()
      #plt.show()
      #for j in range(len(tubecoords)):
      #  ax.scatter(tubecoords[j][0],tubecoords[j][1],tubecoords[j][2])
      #  plt.draw()
      xmax=0.0
      ymax=0.0

      #for j in range(len(tubecoords)):
      #  if abs(tubecoords[j][0])>xmax:
      #    xmax=abs(tubecoords[j][0])
      #  if abs(tubecoords[j][1])>ymax:
      #    ymax=abs(tubecoords[j][1])

      newcellwidth=2.0*tuberad+2.0*self.width3[0]+30.0

      for j in range(len(tubecoords)):
        tubecoords[j][0]+=newcellwidth/2.0
        tubecoords[j][1]+=newcellwidth/2.0
        

      #tubetemp=np.array(list(copy.deepcopy(tubecoords)),dtype='d')
      #for j in range(len(tubetemp)):
      #  tubetemp[j][2]+=self.width3[2] 
      #tubecoords=np.vstack([tubecoords,tubetemp])
      #elemlist.extend(list(copy.deepcopy(elemlist)))
            
      

      fileout=open('NT_TiO2_'+str(int(np.ceil((2*np.pi)/rotang)))+'_132_OC0.dat','w')
      fileout.write('title\n')
      fileout.write('structure with radius '+str(tuberad)+'\n')
      fileout.write('end\n')
      fileout.write('cell\n')
      fileout.write(str(newcellwidth)+'  '+str(newcellwidth)+'  '+str(self.width3[2])+' 90.000 90.000 90.000\n')
      fileout.write('cartesian '+str(len(elemlist))+'\n')
      for j in range(len(tubecoords)):
        fileout.write(elemlist[j]+'\t'+str(tubecoords[j][0])+'\t'+str(tubecoords[j][1])+'\t'+str(tubecoords[j][2])+'\n')
      fileout.write('space 1_a\n')
      fileout.write('supercell 1 1 1\n')
      fileout.write('full\n')
      fileout.close()


  def doRot4(self):
    tuberad=0.0
    ang=0.0
    rotang=0.0
    RotMat=np.zeros((3,3),dtype='d')

    for i in range(15,81):
      ang=(np.pi*(i-2.0))/(2.0*(i))
      rotang=2.0*((np.pi/2.0)-ang)
      tuberad=(self.width4[1]/2.0)*np.tan(ang)

      temp=np.array(list(copy.deepcopy(self.orientpos4)),dtype='d')

      for j in range(len(temp)):
        temp[j][0]=temp[j][0]+tuberad

      tubecoords=temp[:][:]
      elemlist=list(copy.deepcopy(self.corr_atoms)) 
      for j in range(1,i):
        RotMat[0][0]=np.cos(rotang*j)
        RotMat[0][1]=-1.0*np.sin(rotang*j)
        RotMat[1][0]=np.sin(rotang*j)
        RotMat[1][1]=np.cos(rotang*j)
        RotMat[2][2]=1.0
        elemlist.extend(list(copy.deepcopy(self.corr_atoms)))


        for k in range(len(temp)):
          tubecoords=np.vstack([tubecoords,np.dot(RotMat,np.array(temp[k],dtype='d'))])
            
      #fig=plt.figure()
      #ax = fig.add_subplot(111, projection='3d')
      #plt.ion()
      #plt.show()
      #for j in range(len(tubecoords)):
      #  ax.scatter(tubecoords[j][0],tubecoords[j][1],tubecoords[j][2])
      #  plt.draw()
      xmax=0.0
      ymax=0.0

      #for j in range(len(tubecoords)):
      #  if abs(tubecoords[j][0])>xmax:
      #    xmax=abs(tubecoords[j][0])
      #  if abs(tubecoords[j][1])>ymax:
      #    ymax=abs(tubecoords[j][1])

      newcellwidth=2.0*tuberad+2.0*self.width4[0]+30.0

      for j in range(len(tubecoords)):
        tubecoords[j][0]+=newcellwidth/2.0
        tubecoords[j][1]+=newcellwidth/2.0
        
      tubetemp=np.array(list(copy.deepcopy(tubecoords)),dtype='d')
      initelems=list(copy.deepcopy(elemlist))
      for tubeINC in range(1,6):
        tubetemp1=np.array(list(copy.deepcopy(tubetemp)),dtype='d')
        for j in range(len(tubetemp1)):
          tubetemp1[j][2]+=(self.width4[2]*tubeINC)
        tubecoords=np.vstack([tubecoords,tubetemp1])
        elemlist.extend(list(copy.deepcopy(initelems)))


      delx=0.0
      for j in range(36):
        xmax=np.random.randint(int(((newcellwidth/2.0)-(tuberad+self.width4[0]+delx))*100),int(((newcellwidth/2.0)+(tuberad+self.width4[0]+delx))*100))/100.0
        ymax=((-1)**(np.random.randint(0,2)))*np.sqrt(((tuberad+self.width4[0]+delx)**2.0)-((xmax-(newcellwidth/2.0))**2))+((newcellwidth/2.0))
        zmax=np.random.randint(0,int(self.width4[2]*5*100))/100.0     
        quickvec=np.zeros(3,dtype='d')
        quickvec[0]=xmax
        quickvec[1]=ymax
        quickvec[2]=zmax
        tubecoords=np.vstack([tubecoords,quickvec])
        elemlist.append('li1')


      #fig=plt.figure()
      #ax=fig.add_subplot(111,projection='3d')
      #plt.ion()
      #plt.show()
      #for j in range(len(tubecoords)):
      #  ax.scatter(tubecoords[j][0],tubecoords[j][1],tubecoords[j][2])
      #  plt.draw()
            

      fileout=open('NT_TiO2_'+str(int(np.ceil((2*np.pi)/rotang)))+'_213_OC0.dat','w')
      fileout.write('title\n')
      fileout.write('structure with radius '+str(tuberad)+'\n')
      fileout.write('end\n')
      fileout.write('cell\n')
      fileout.write(str(newcellwidth)+'  '+str(newcellwidth)+'  '+str(self.width4[2]*5)+' 90.000 90.000 90.000\n')
      fileout.write('cartesian '+str(len(elemlist))+'\n')
      for j in range(len(tubecoords)):
        fileout.write(elemlist[j]+'\t'+str(tubecoords[j][0])+'\t'+str(tubecoords[j][1])+'\t'+str(tubecoords[j][2])+'\n')
      fileout.write('space 1_a\n')
      fileout.write('supercell 1 1 1\n')
      fileout.write('full\n')
      fileout.close()


  def doAllRots(self):
    self.doRot1()
    self.doRot2()
    self.doRot3()
    self.doRot4()



''' 
This next class will open an olcao.mi and get the resulting
cell vectors, angles, atom names, positions.


''' 

class olcaomiFile():
  def __init__(self):
    self.atomicList=[]
    self.positionList=[]
    self.cellVec=np.zeros(3,dtype='d')
    self.cellAng=np.zeros(3,dtype='d')
    self.fracorcart=''
    self.numAtoms=0

  def getFile(self):
    olcaomi=open('olcao.mi','r')
    allLines=olcaomi.readlines()
    olcaomi.close()


    for i in range(len(allLines)):

      if len(self.atomicList)==self.numAtoms and self.numAtoms!=0:
        break
    
      if allLines[i].strip()=='cell':
        tempLine=(allLines[i+1].strip()).split()
        self.cellVec[0]=tempLine[0]
        self.cellVec[1]=tempLine[1]
        self.cellVec[2]=tempLine[2]
        self.cellAng[0]=tempLine[3]
        self.cellAng[1]=tempLine[4]
        self.cellAng[2]=tempLine[5]

        tempLine=(allLines[i+2].strip()).split()
        self.fracorcart=tempLine[0]
        self.numAtoms=int(tempLine[1])        

        if self.fracorcart=='cart':
          for j in range(i+3,i+self.numAtoms+3):
            tempLine=(allLines[j].strip()).split()
            self.atomicList.append(tempLine[0])     
            self.positionList.append([tempLine[1],tempLine[2],tempLine[3]])        
        else:
          for j in range(i+3,i+self.numAtoms+3):
            tempLine=(allLines[j].strip()).split()
            self.atomicList.append(tempLine[0])     
            self.positionList.append([float(tempLine[1])*self.cellVec[0],float(tempLine[2])*self.cellVec[1],float(tempLine[3])*self.cellVec[2]])        
 

print "Beginning Nanotube Creation\n"

'''     
Basic Idea:
 1) Acquire the structure
 2) Acquire desired parameters
 3) Place initial piece with correct orientation and at correct distance
 4) Make first ring
 5) Replicate ring upwards
'''



#Need to add in a spot to take in chirality

newOLCAOfile=olcaomiFile()
newOLCAOfile.getFile()

newstruct=initStructure(newOLCAOfile.cellVec[0],newOLCAOfile.cellVec[1],newOLCAOfile.cellVec[2])
newstruct.addEntireStructure(newOLCAOfile.atomicList,newOLCAOfile.positionList)

newstruct.setAllThreePos()
newstruct.doAllRots()





print "Time to complete "+str(time.time()-start)+"s."




