#!/usr/bin/env python3

"""
Author:  Patrick Ryan Thomas
         UMKC - CPG
         prtnpb@mail.umkc.edu

Purpose: The purpose of this script is take in relativistic
         orbital data from grasp2K and perform the following:

         1) WEIGTHED LINEAR INTERPOLATE all orbitals (large and small)
              to be expressed on the same radial grid.

         2) PARARMETER SWEEP over the maximum exponential coefficient
              and number of terms for each l quantum number type.

            a) A single LU DECOMPOSITION can be used to solve
                 all orbitals in a single set.

            b) Exponential coefficients are determined via
                 geometric series by maximum, minimum, and
                 number of terms.

         3) Output orbitals in a format usable by OLCAO.

"""


import numpy as np
import scipy
import time
import copy
import os
import sys
import math
import pandas as pd
import make_veusz_graph
import org_radwavefn
import plotly.express as px

from scipy.linalg import lu_solve 

start=time.time()

#==========================
# PRINTING FUNCTIONS:
#==========================

def print_help():
  print("\n=====================")
  print("HELP FOR INTERFIT.PY")
  print("=====================\n")
  print("Operation modes:")
  print( "  -shrink rc sig")
  print( "    Performs a shrinking function operation where rc")
  print("     is the critical radii and sig is the shrinking")
  print("     parameter.\n")
  print("           |          2")
  print("           |    (r-rc) ")
  print("           |    ------ ")
  print("           |         2 ")
  print("           |    2*sig  ")
  print("     s(r)= |1-e         r .leq. rc ")
  print("           |0           r > rc\n")
  print("     A file (shrink.inf) must be present with a list")
  print("      of all orbitals and corresponding rc sig in comma")
  print("      delimited columns.\n\n")
  print( "  -comp 1/2/4")
  print("    Designates the amount of components to generate")
  print("     one, two, or four.")
  print("    1: Averages together to large components into a")
  print("       single orbital.")
  print("       -comp 1")
  print("    2: Fits only the large components.")
  print("       -comp 2")
  print("    3: Fits all four components.")
  print("       -comp 4\n")

def print_line(number_lines):
  for i in range(number_lines):
    print("==================================================================================")

def print_callout(callout_text):
  print_line(1)
  print(callout_text)

def printStats(MaxA,numberTerms,Weight):
  print(str("%.8E"%MaxA) + "\t" + str(numberTerms[0]) + "\t" +
        str(numberTerms[1]) + "\t" + str(numberTerms[2]) + "\t" +
        str(numberTerms[3]) + "\t" + str(numberTerms[4]) + "\t" +
        str(Weight) + "\n")

def printProgress(minAlpha,maxAlpha,currentA,step):
  perc=((float(currentA)-float(minAlpha))/(float(maxAlpha)-float(minAlpha)))*100.0
  if perc%10==0:
    print(str("%.2f"%perc)+"%\t("+str(minAlpha)+"\t- "+str(currentA)+" -\t"+str(maxAlpha)+")")
  
 

#==========================
# ELEMENT INFO FUNCTIONS:
#==========================
def getAtomNum(symb):
  elems={'H': 1 , 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109}
  return elems[symb]

def lktol(k):
  if k<0:
    return (-1*k)-1
  else:
    return k

def sktol(k):
  if k<0:
    return k*-1
  else:
    return k-1

def ltosk(l):
  if l==0:
    return 1,1
  else:
    return -1*l,l+1

def ltok(l):
  if l==0:
    return -1,-1
  else:
    return l,-1*(l+1)

def ltolsym(l):
  if l==0:
    return "s"
  elif l==1:
    return "p"
  elif l==2:
    return "d"
  elif l==3:
    return "f"
  elif l==4:
    return "g"

def get_max_alpha(Z):
  alphas={'1': 10000.0, '2': 10000.0, '3': 50000.0, '4': 50000.0, '5': 50000.0, '6': 50000.0, '7': 50000.0, '8': 50000.0, '9': 50000.0, '10': 50000.0, '11': 100000.0, '12': 100000.0, '13': 100000.0, '14': 200000.0, '15': 100000.0, '16': 100000.0, '17': 100000.0, '18': 100000.0, '19': 500000.0, '20': 500000.0, '21': 500000.0, '22': 500000.0, '23': 500000.0, '24': 500000.0, '25': 500000.0, '26': 500000.0, '27': 500000.0, '28': 500000.0, '29': 500000.0, '30': 1000000.0, '31': 1000000.0, '32': 1000000.0, '33': 1000000.0, '34': 1000000.0, '35': 1000000.0, '36': 1000000.0, '37': 5000000.0, '38': 5000000.0, '39': 5000000.0, '40': 5000000.0, '41': 5000000.0, '42': 5000000.0, '43': 5000000.0, '44': 5000000.0, '45': 5000000.0, '46': 5000000.0, '47': 5000000.0, '48': 5000000.0, '49': 5000000.0, '50': 5000000.0, '51': 5000000.0, '52': 5000000.0, '53': 5000000.0, '54': 5000000.0, '55': 10000000.0, '56': 10000000.0, '57': 10000000.0, '58': 10000000.0, '59': 10000000.0, '60': 10000000.0, '61': 10000000.0, '62': 10000000.0, '63': 10000000.0, '64': 10000000.0, '65': 10000000.0, '66': 10000000.0, '67': 10000000.0, '68': 10000000.0, '69': 10000000.0, '70': 10000000.0, '71': 50000000.0, '72': 50000000.0, '73': 50000000.0, '74': 50000000.0, '75': 50000000.0, '76': 50000000.0, '77': 50000000.0, '78': 50000000.0, '79': 50000000.0, '80': 50000000.0, '81': 50000000.0, '82': 50000000.0, '83': 50000000.0, '84': 50000000.0, '85': 50000000.0, '86': 50000000.0, '87': 50000000.0, '88': 50000000.0, '89': 50000000.0, '90': 50000000.0, '91': 50000000.0, '92': 50000000.0, '93': 50000000.0, '94': 50000000.0, '95': 50000000.0, '96': 50000000.0, '97': 50000000.0, '98': 50000000.0, '99': 50000000.0, '100': 50000000.0, '101': 50000000.0, '102': 50000000.0, '103': 50000000.0,}

  return alphas[str(Z)]

#==========================
# FILE OPERATION FUNCTIONS:
#==========================


def getFileList(ext):
  # This function returns a list of files names in the current directory
  #  ending with the desired file extension.
  inFiles=[]
  inFiles+=[each for each in os.listdir('.') if each.endswith(str(ext))]
  return inFiles


def getFileContents(filename):
  f=open(filename,'r')
  content=f.readlines()
  f.close()
  return content

def parseFileName(filename):
  k=(((filename.split('.'))[0]).split('_'))[4]
  n=(((filename.split('.'))[0]).split('_'))[2]
  atomSym=(((filename.split('.'))[0]).split('_'))[0]
  return atomSym,n,k




#==========================
# INTERPOLATION FUNCTIONS:
#==========================
def intNewCalc(x1,x2,y1,y2,xn):
  """
  This function calculates the linearly interpolated y-value 
  of a given point (xn) on a line between points, (x1,y1)
  and (x2,y2) Utilized point slope form of a line, i.e.,
     (y2 - y1)
  yn=--------- (xn - x1) + y1    
     (x2 - x1)
  """
  return float(((y2-y1)/(x2-x1))*(xn-x1)+y1)  


def getPoints(rval,Rvals):
  """
  This function takes in a radial grid value (rval) that does not
  exits in the original list (Rvals). Starting in the middle
  of the list, Rvals, we check to see if rval is greater
  than the middle value of list. If it is, we check to see
  if it is less than the next value in the list. If it is,
  we've found the lower value needed. If it is not, perform
  the same check on the next item in the list. However, if we
  reach the end of the list, use the previous value.
  """

  # Make a copy of the passed list.
  RO=copy.deepcopy(Rvals)

  #Initialize output variable
  xval=0
  #Initialize marker to be at middle of the list.
  v=int(len(RO)/2)

  # Loop until RO[v]<rval<RO[v+1]    
  while(1):
    """
    If we are at the end of the list, use the previous
    value.
    """
    if v==len(RO)-1:
      xval=v-1
      break

    # Check if value is bigger than midpoint
    if rval>RO[v]:
      #Check if value is less than midpoint + 1. If yes, done.
      if rval<RO[v+1]:
        xval=v
        break
        #If the value is not less than midpoint, increment by one.
      else:
        v+=1
    #If not bigger than the midpoint, decrease position by one.
    else:
      v=v-1              
  return xval


def doInterp(R,LDR,SDR):
  """
  This function will take a set of radial values in
  the global R list and create the corresponding common
  list of all radial positions for interpolation. Then,
  using the values in the global LDR and SDR lists, the function
  will compute a linear interpolation on the points and
  return lists of new radial values and new large and
  small component values that have been interpolated on the
  new radial grid.
  """
  # Create a temp arrays to hold the new values
  r=[]
  l=[]
  s=[]
  tl=[]
  ts=[]
  interpTime=time.time()  

  # For set of radial values, add it to the temp
  #  list.
  for i in range(len(R)):
    r.extend(R[i])
  # Remove duplicate values.
  r=list(set(r))
  
  # Sort the values
  r.sort()
  
  
  # loop through each large and small component
  #  corresponding to a particular radial grid to generate
  #  new values.
  tl=[None]*len(r)
  ts=[None]*len(r)
  n=0
  print_line(1)
  print ("BEGINNING INTERPOLATION:\nThere are " + str(len(R)) +
         " pairs of orbitals to be interpolated and " + str(len(r)) +
         " points in the new grid.")

  # Loop through each pair of radial grid, large divided by r,
  #   and small divided by r.
  for i in range(len(R)):
    # Loop through each newly calculated radial grid point.
    for j in range(len(r)):
      # If the new grid value is in the original grid,
      #   put the corresponding large and small components
      #   in the temporary lists.
      if r[j] in R[i]:
        tl[j]=LDR[i][R[i].index(r[j])]
        ts[j]=SDR[i][R[i].index(r[j])]
      # If the new grid value is not in the original grid,
      #   we perform a linear interpolation using original grid
      #   values on either side of the new value.
      else:
        n=getPoints(r[j],R[i])
        tl[j]=intNewCalc(R[i][n],R[i][n+1],LDR[i][n],LDR[i][n+1],r[j])
        ts[j]=intNewCalc(R[i][n],R[i][n+1],SDR[i][n],SDR[i][n+1],r[j])
      
    # Append the temp arrays to the returnable arrays
    l.append(copy.deepcopy(tl))
    s.append(copy.deepcopy(ts))
      
    # Re-initalize the temp arrays to empty.
    tl=[None]*len(r)
    ts=[None]*len(r)
  print ("Interpolation Complete: ("+str(time.time()-interpTime)+"s).\n")
  return r,l,s

  
#==========================
# FITTING FUNCTIONS:      
#==========================

def getAlphaList(amin,amax,n,mode):
  """
  This function makes a geometric series of alphas
  given n (number of terms) and (amin,amax) which 
  are the minimum and maximum alpha in the series
  respectively.
  """
  a=[]
  
  if mode==1: 
    for i in range(1,n+1):
      r=amin*((float(amax/amin))**((i-1.0)/(n-1.0)))
      a.append(r)
    return a
  elif mode==2:
    diff=float((amax-amin)/(n-1))
    for i in range(0,n):
      r=amin+(diff*i)
      a.append(r)
    return a

def LUDecomp(A):
  """
  This function performs a pivotless LU Decomposition of matrix A
  and returns LU (top triangle U, bottom triangle L) which is the format
  required by scipy.linalg.lu_solve.
  """
  # Define L and U matrices respectively.
  L=np.zeros((len(A),len(A)),dtype="d")
  U=copy.deepcopy(A)
  
  # Loop through all of the columns.
  for j in range(len(A)-1):
    # Loop through values in the lowerr diagonal.
    for i in range(j,len(A)-1):
      
      tc=(U[i+1][j]/U[j][j])
      if U[j][j]==0 or tc==0:
        continue
      else:
        U[i+1][:]=U[i+1][:]-tc*U[j][:]
        L[i+1][j]=tc
      
      # If a traditional LU Decomp is needed uncomment the next two
      #  lines and the subtraction by identity.
      #L[i][i]=1
      #L[i+1][i+1]=1

  LU=(L+U)#-np.identity(len(A))
  
  return LU


def getA(numS,radialGrid,GCoeffs):
  """
  As part of the least norm squares solution to Ax=B, we define
                         2
           -alpha[j]*r[i]
  A[i][j]=e
  """

  AMAT=np.zeros((len(radialGrid),numS),dtype="d")
  for i in range(numS):
    for j in range(len(radialGrid)):
      AMAT[j][i]=np.exp(-1.0*GCoeffs[i]*(radialGrid[j]**2.0))

  return AMAT


def getCoeffs(ATA,A,AT,B,angVals,als,numSP,numD,numF,numG,numOrbs):
  """
  This function solves for the minimal norm
  least squares solution:

      A^T*Ax=A^TB
      
  where B consists of all orbital function values
  arranged in columns.
  
  As part of the process, when we multiply A^T and
  B, we only use enough rows of A^T to equal the number
  of terms used in the fitting. The rest of the column
  will be left as zeros.
  
  If the g-type orbitals are being fitted with ten terms, then we
  take the ten by ten block of LU found from LU Decomposition and
  lu_solve with said matrix and the corresponding block of g-type
  solutions but not the zeros at the bottom of the columns.
  """
  # Make copies to be sure we aren't editing original values.
  fn=copy.deepcopy(B)
  alphas=copy.deepcopy(als)
  ludc=copy.deepcopy(ATA)
  amat=copy.deepcopy(A)
  atmat=copy.deepcopy(AT)
  angVals=map(int,angVals)
  numterms=[0]*len(angVals)
  termlist=[0,0,0,0,0]
  termlist[0]=numSP
  termlist[1]=numSP
  termlist[2]=numD
  termlist[3]=numF
  termlist[4]=numG  

  for i in range(len(angVals)):
    if angVals[i]==0 or angVals[i]==1:
      numterms[i]=numSP
    elif angVals[i]==2:
      numterms[i]=numD
    elif angVals[i]==3:
      numterms[i]=numF
    elif angVals[i]==4:
      numterms[i]=numG


  # Make a new matrix such that we only multiply
  #  the rows corresponding to number of terms
  #  that we want to fit with. Note that a vector
  #  will need to be made that has
  #  the corresponding number of terms for each orbital
  #  or the logic altered.
  Cmat=np.zeros((len(alphas),len(numterms)),dtype="d")

  for j in range(len(numterms)):
    for i in range(numterms[j]):
      Cmat[i][j]=np.dot(atmat[i,:],fn[:,j])


  # Based on the maximum angular momentum of the 
  #  of the orbitals being fitted, we solve
  #  the corresponding systems and save the coefficients
  #  in the coefficients matrix. 
  coeffs=np.zeros((len(alphas),len(numterms)),dtype="d") 

  for i in range(len(numOrbs)-1,-1,-1):
    if numOrbs[i]!=0:
      x=termlist[i]
      y1=sum(numOrbs[i:])-numOrbs[i]
      y2=y1+numOrbs[i]
      PIV=np.matrix(list(xrange(x)))
      coeffs[:x,y1:y2]=scipy.linalg.lu_solve((ludc[:x,:x],PIV),Cmat[:x,y1:y2])

  RMSE=[0]*len(numterms)
  testMat=np.dot(amat,coeffs)
  RM=np.dot(amat,coeffs)-fn

  
  for j in range(len(numterms)):
    for i in range(numSP):
      RMSE[j]+=RM[i,j]**2.0
    RMSE[j]=np.sqrt(RMSE[j])
 
  return coeffs,RMSE

def checkLess(c,l):
  tol=10**-2
  for i in range(len(c)):
    if c[i]-l[i]>0 and abs(c[i]-l[i])>tol:
      return False
  return True

#==========================
# SWEEPING FUNCTIONS
#==========================
"""
This function performs a sweep over the number of terms in
  the expansion for each angular momementum type as well
  as the maximum alpha in the geometric series (or other
  method, see function getAlphaList(amin,amax,n,mode).
  
  Key concepts:
  * Checking for improvement:
      1) The methodolgy has two checks, first the RMSE
      for each orbital must either decrease or decrease
      within the specified tolerance (see checkLess(c,l)).

      2) Because the number of terms in the cartesian
      expansion of the orbitals increases with an increase
      in angular momentum, we want to weight the improvement
      by a reduction in terms.

      weight = numTerm_g*9+numTerm_f*7+numTerm_d*5+numTerms_sp*4

  Key variables to customize sweeping are:
  1) stepSize: Changing the decimal multiplied by
               maxAlpha will change resolution of the
               parameter sweep. Because of the large
               difference in size between reference
               maximum alphas, percentages make
               a more logical choice.

  2) uM,lM: These variables stand for lower multiplier
            and upper multiplier. Again, these work as
            percentages of the reference maximum alpha
            value.

  3) maxNT,minNT: These variables specify the maximum starting
                  number of terms and the lowest to go.
"""
def doSweep(R,B,orbs,lvals,numls):
  sweeptime=time.time()  # Start the sweeping clock
  minAlpha=0.12          # This comes from OLCAO standard (numerical stability)
  maxAlpha=float(get_max_alpha(Z)) #Get the OLCAO stand. maximum alpha
  stepSize=int(maxAlpha*0.01)    #Calculate the stepsize
  cA=[]              # Current list of coefficients
  fA=[]              # Best list of coefficients
  cAlpha=[]          # Current list of exponential alphas
  fAlpha=[]          # Best list of exponential alphas
  cRMSE=0.0          # Current RMSE
  fRMSE=0.0          # Best RMSE
  cWght=0            # Current weight
  fWght=1000000000   # Best weight (start extremely high for first run)
  nT=[0,0,0,0,0]     # Number of terms
  maxNT=20           # Maximum number of terms in every orbital
  minNT=8            # Minimum number of terms in every orbital
  uM=10              # Upper multiplier
  lM=0.3            # Lower multiplier
  uB=int(uM*maxAlpha) # Upper bound for maxa loop
  lB=int(lM*maxAlpha) # Lower bound for maxa loop
  maxl=int(lvals[0]) # Get the maximum l in the set of orbitals
  itNum=0               # Variable to track number of iterations
  bMat=np.transpose(np.matrix(B)) #Make a copy of the functions values list.
  noImprove=0        # Variable to see if the reduction of terms is working.

  impCount=0         # Count number of improvements.

  print ("\t\tC  A  L  C  U  L  A  T  I  O  N  S")
  print_line(1)
  print("BEGINNING SWEEPING:")
  print("Maximum Alpha Range: " + str(int(lM*maxAlpha)) + " to " +
        str(int(uM*maxAlpha)) + " (Step size: " + str(stepSize) + ")")

  # If the maximum is g-type
  if maxl==4:
    print("Max. l Quantum Num.: g-type")
    print_line(1)
    # Loop over maximum alpha values
    for mA in range(lB,uB,stepSize):
      printProgress(lB,uB,mA,stepSize)
      #Loop over s- and p- orbital values
      for sp in range(maxNT,minNT,-1):
        # Get the list of exponential alphas.
        cAlpha=getAlphaList(minAlpha,mA,sp,1)
        # Get the MNLSS solution A
        A=getA(sp,R,cAlpha)
        # Get the transpose of A.
        AT=np.transpose(copy.deepcopy(A))
        # Dot product of transpose of A and A.
        ATA=np.dot(AT,A)
        # Perform an LU decomposition of ATA.
        #  This is key because as long as the number of terms
        #   for s- an p-type orbitals don't change, we can
        #   use the same LU decomposition to solve.
        LU=LUDecomp(ATA)
        # Loop over d,f, and g number of terms:
        for d in range(sp,minNT,-1):
          for f in range(d,minNT,-1):
            for g in range(f,minNT,-1):
              # Calculate the current weight.
              cWght=4*sp+5*d+7*f+9*g
              # If two reductions in the number of terms don't make
              #  an improvement, break the loop and try for the next.
              if noImprove==2:
                noImprove=0
                break
              if cWght>fWght:
                continue
              # Solve for the new coefficients and RMSE
              cA,cRMSE=getCoeffs(LU,A,AT,bMat,lvals,cAlpha,sp,d,f,g,numls)
              # Check to see if the RMSE has been initialized
              #  or if there is an improvement.
              if fRMSE==0.0 or checkLess(cRMSE,fRMSE):
                # Make the current weight, RMSE, coefficents,
                #  and exponential alphas the best. Store
                #  the corresponding number of terms.
                #  re-initialize improvement check.
                fWght=cWght
                fRMSE=copy.deepcopy(cRMSE)
                fA=copy.deepcopy(cA)
                fAlpha=copy.deepcopy(cAlpha)
                nT[0]=sp
                nT[1]=sp
                nT[2]=d
                nT[3]=f
                nT[4]=g
                noImprove=0
                impCount+=1
              else:
                noImprove+=1
              itNum+=1

  elif maxl==3:
    print("Max. l Quantum Num.: f-type")
    print_line(1)
    for mA in range(lB,uB,stepSize):
      printProgress(lB,uB,mA,stepSize)
      for sp in range(maxNT,minNT,-1):
        cAlpha=getAlphaList(minAlpha,mA,sp,1)
        A=getA(sp,R,cAlpha)
        AT=np.transpose(copy.deepcopy(A))
        LU=LUDecomp(np.dot(AT,A))
        for d in range(sp,minNT,-1):
          for f in range(d,minNT,-1):
            if noImprove==2:
              noImprove=0
              break
            cWght=4*sp+5*d+7*f
            if cWght>fWght:
              continue
            cA,cRMSE=getCoeffs(LU,A,AT,bMat,lvals,cAlpha,sp,d,f,0,numls)
            if fRMSE==0.0 or checkLess(cRMSE,fRMSE):
              fWght=cWght
              fRMSE=copy.deepcopy(cRMSE)
              fA=copy.deepcopy(cA)
              fAlpha=copy.deepcopy(cAlpha)
              nT[0]=sp
              nT[1]=sp
              nT[2]=d
              nT[3]=f
              noImprove=0
              impCount+=1
            else:
              noImprove+=1
            itNum+=1

  elif maxl==2:
    print("Max. l Quantum Num.: d-type")
    print_line(1)
    for mA in range(lB,uB,stepSize):
      printProgress(lB,uB,mA,stepSize)
      for sp in range(maxNT,minNT,-1):
        cAlpha=getAlphaList(minAlpha,mA,sp,1)
        A=getA(sp,R,cAlpha)
        AT=np.transpose(copy.deepcopy(A))
        LU=LUDecomp(np.dot(AT,A))
        for d in range(sp,minNT,-1):
          if noImprove==2:
            noImprove=0
            break
          cWght=4*sp+5*d
          if cWght>fWght:
            continue
          cA,cRMSE=getCoeffs(LU,A,AT,bMat,lvals,cAlpha,sp,d,0,0,numls)
          if fRMSE==0.0 or checkLess(cRMSE,fRMSE):
            fWght=cWght
            fRMSE=copy.deepcopy(cRMSE)
            fA=copy.deepcopy(cA)
            fAlpha=copy.deepcopy(cAlpha)
            nT[0]=sp
            nT[1]=sp
            nT[2]=d
            noImprove=0
            impCount+=1
          else:
            noImprove+=1
          itNum+=1

  elif maxl==1 or maxl==0:
    print("Max. l Quantum Num.: s- or p-type")
    print_line(1)
    for mA in range(lB,uB,stepSize):
      printProgress(lB,uB,mA,stepSize)
      for sp in range(maxNT,minNT,-1):
        cAlpha=getAlphaList(minAlpha,mA,sp,1)
        A=getA(sp,R,cAlpha)
        AT=np.transpose(copy.deepcopy(A))
        LU=LUDecomp(np.dot(AT,A))
        cWght=4*sp
        if cWght>fWght:
          continue
        cA,cRMSE=getCoeffs(LU,A,AT,bMat,lvals,cAlpha,sp,0,0,0,numls)

        if fRMSE==0.0 or checkLess(cRMSE,fRMSE):
          fWght=cWght
          fRMSE=copy.deepcopy(cRMSE)
          fA=copy.deepcopy(cA)
          fAlpha=copy.deepcopy(cAlpha)
          nT[0]=sp
          nT[1]=sp
          impCount+=1
        itNum+=1

  print_line(1)
  print("Max. Alpha\ts\tp\td\tf\tg\tWeight")
  printStats(fAlpha[len(fAlpha)-1],nT,fWght)
  print("Total Iterations: "+str(itNum))
  print("Sweeping Complete: ("+str(time.time()-sweeptime)+")")
  print_line(1)
  return fA,fAlpha,fRMSE,nT


#==========================
# OPERATIONAL  FUNCTIONS:      
#==========================


def graphFunctions(orbitals,lvalues,coeffs,alphas):
  q=open('equations.dat','w')
  temptext=""
  for i in range(len(orbitals)):
    temptext=str(orbitals[i])+' '
    for j in range(len(alphas)):
      temptext = temptext + "((" + str('%.8E'%coeffs[j,i]) + ")*(x**(" + \
          str(lvalues[i]) + "))*e**(-(" + str('%.8E'%alphas[j]) + ")*x**2))"
      if j == (len(alphas)-1):
        temptext=temptext + '\n'
      else:
        temptext=temptext + '+'
    q.write(temptext)
    q.write('\n')

  q.close()


def writeFunctions(orbitals,lvalues,coeffs,alphas):
  h=open('orbitals.txt','a')
  h.write('Exponential Alphas\n')
  for i in range(0,len(alphas)):
    h.write(str('%.8E'%alphas[i]))
    if (i+1)%4==0 and i+1>=4:
      h.write('\n')
    else:
      h.write('\t')

  h.write('\n\n')

  h.write("Orbitals\n")
  for i in range(len(orbitals)-1,-1,-1):
    h.write(str(orbitals[i])+'\n')
    for j in range(len(alphas)):
      h.write(str('%.8E'%coeffs[j,i]))
      if ((j+1)%4==0 and j+1>=4) or j+1==len(alphas):
        h.write('\n')
      else:
        h.write('\t')
  h.close()




# The function returns a geometic series of exponental
#  coefficients.
def get_alpha_list(amin,amax,n,mode):

  a=[]
  # If the mode passed is geometric, create a list from
  #  amin to amax with n number of terms.
  if mode=='geometric':
    for i in range(1,n+1):
      r=amin*((float(amax/amin))**((i-1.0)/(n-1.0)))
      a.append(r)
    return a


# Function to create A matrix for AX=B
def get_A(numS,radialGrid,GCoeffs):

  AMAT=np.zeros((len(radialGrid),numS),dtype="d")

  for j in range(numS):
    for i in range(len(radialGrid)):
      AMAT[i,j]=np.exp(-1.0*GCoeffs[j]*(radialGrid[i]**2.0))

  return AMAT


# Function to do a pivotless LU decomposition for scipy's
#  lu_solve routine.
def pivotless_lu_decomp(A):

  dim=np.shape(A)[0]

  L=np.zeros((dim,dim),dtype="d")
  U=copy.deepcopy(A)

  for j in range(dim-1):
    for i in range(j,dim-1):
      tc=(U[i+1,j]/U[j,j])
      if U[j,j]==0 or tc==0:
        continue
      else:
        U[i+1,:]=U[i+1,:]-tc*U[j,:]
        L[i+1,j]=tc

        # Adding this in to make truly upper triangular
        if U[i+1,j]<=1e-13 and U[i+1,j]>=-1e-13:
          U[i+1,j]=0.0

  return L+U


class orbital_fitting:

  def __init__(self):
    self.apply_shrink=False
    self.shrink_critical_radius=0.0
    self.shrink_sigma=0.0
    self.number_components=0

    self.grasp_description=org_radwavefn.atomic_system('isodata', 'rwfn.out')

    self.min_alpha=0.12
    self.max_alpha=get_max_alpha(str(self.grasp_description.atomic_info.atomic_number))

    self.fitting_columns=['min_alpha', 'max_alpha', 'num_terms', 'num_s',
                          'num_p', 'num_d', 'num_f', 'num_g', 'weight',
                          'total_RMSE', 'RMSE', 'time']
    self.fitting_results=[]

    self.RMSE=1E100
    self.weight_factors=[1,3,5,7,9]
    self.weight=1E100
    self.fits=[]
    self.number_terms=[0,0,0,0,0]
    self.fitted_max_alpha=0

    self.orb_fit_tolerance=10**-3
    self.alpha_step_size=int(self.max_alpha*0.001)

    self.loop_max_alpha=np.arange(0.3*self.max_alpha,
                                  20*self.max_alpha + self.alpha_step_size,
                                  self.alpha_step_size)
    self.max_num_terms=25
    self.min_num_terms=7


  def parse_input(self):

    parse_time=time.time()

    if len(sys.argv)==1:
      print("No options were selected. Please try again or -help.\n")
      sys.exit()

    else:
      #Check for help
      if "-help" in sys.argv:
        print_help()
        if sys.argv[1]=="-help":
          print_help()
          sys.exit()
      #Check for -comp 1/2/4
      if '-comp' in sys.argv:
        self.number_components=int(sys.argv[int(sys.argv.index('-comp'))+1])

      #Check for -shrink rc sig
      if '-shrink' in sys.argv:
        self.shrink_critical_radius=float(sys.argv[int(sys.argv.index('-shrink'))+1])
        self.shrink_sigma=float(sys.argv[int(sys.argv.index('-shrink'))+2])
        self.apply_shrink=True

    print("  Parsing Input Complete:\t" +
          str('%.8E'%(time.time()-parse_time)) + "s. / " +
          str('%.8E'%(time.time()-start)) + "s.")


  def organize_orbitals(self):
    organize_time=time.time()

    if self.number_components==1:
      self.grasp_description.orbital_info.single_component()
    elif self.number_components==2:
      self.grasp_description.orbital_info.two_component()
    elif self.number_components==4:
      self.grasp_description.orbital_info.four_component()

    if self.apply_shrink:
      shrinking_time=time.time()      
      self.grasp_description.orbital_info.apply_shrinking_function(self.shrink_critical_radius,self.shrink_sigma)
      print("  Shrinking Orbitals Complete:\t"+str('%.8E'%(time.time()-shrinking_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")
    else:
      for i in range(len(self.grasp_description.orbital_info.fitting_l_list)):
        for j in range(len(self.grasp_description.orbital_info.interpolated_radial_grid)):
          self.grasp_description.orbital_info.fitting_functions[i][j]=self.grasp_description.orbital_info.fitting_functions[i][j]/(self.grasp_description.orbital_info.interpolated_radial_grid[j]**float(self.grasp_description.orbital_info.fitting_l_list[i]))

    print("  Organizing Orbitals Complete:\t"+str('%.8E'%(time.time()-organize_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")


  def setup_calculation(self):
    calc_time=time.time()
    print(' Calculation Setup:')
    self.parse_input()
    self.organize_orbitals()
    print(" Calculation Setup Complete:\t"+str('%.8E'%(time.time()-calc_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")




  def get_prefactor_coefficients(self,LU_matrix,A_matrix,AT_matrix,B_matrix,num_terms,alpha_list):
    prefactor_time=time.time()


    alphas=copy.deepcopy(alpha_list)
    A=copy.deepcopy(A_matrix)
    ATA=copy.deepcopy(LU_matrix)
    AT=copy.deepcopy(AT_matrix)
    functions=copy.deepcopy(B_matrix)
    num_orbitals=copy.deepcopy(self.grasp_description.orbital_info.fitting_lvals)

    RMSE=[]


    #print(self.grasp_description.orbital_info.fitting_lvals)

    function_number_terms=[]

    for l in self.grasp_description.orbital_info.fitting_l_list:

      function_number_terms.append(num_terms[l])

    AT_B=np.zeros((len(alphas),len(function_number_terms)),dtype="d")

    for j in range(len(function_number_terms)):
      for i in range(function_number_terms[j]):
        AT_B[i,j]=np.dot(AT[i,:],functions[:,j])



    coefficient_matrix=np.zeros((len(alphas),len(function_number_terms)),dtype="d") 

    for i in range(len(num_orbitals)-1,-1,-1):
      if num_orbitals[i]!=0:
        x=num_terms[i]
        y1=sum(num_orbitals[i:])-num_orbitals[i]
        y2=y1+num_orbitals[i]
        PIV=np.arange(x)

        coefficient_matrix[:x,y1:y2]=scipy.linalg.lu_solve((ATA[:x,:x],PIV),AT_B[:x,y1:y2])

    RMSE=[0]*len(function_number_terms)
    #testMat=np.dot(amat,coeffs)
    RM=np.dot(A,coefficient_matrix)-functions

    
    '''
    with open('RM_matrix.txt','a') as f:
      matshape=np.shape(RM)

      for i in range(matshape[0]):
        for j in range(matshape[1]):
          f.write(str('%.1E'%RM[i,j])+" ")
        f.write('\n')
      f.write('\n\n')
    '''
    for j in range(len(function_number_terms)):
      for i in range(num_terms[0]):
        RMSE[j]+=RM[i,j]**2.0
      RMSE[j]=np.sqrt(RMSE[j])

    total_RMSE=0.0

    for i in range(len(RMSE)):
      total_RMSE+=self.weight_factors[self.grasp_description.orbital_info.fitting_l_list[i]]*RMSE[i]

    #print('%.8E'%total_RMSE)

    return (coefficient_matrix,total_RMSE,RMSE,time.time()-prefactor_time)




  def fitting_procedure(self):
    sweep_time=time.time()

    best_weight=self.weight
    current_weight=0.0

    best_RMSE=self.RMSE
    current_RMSE=0.0
    total_RMSE=0.0
    
    convergence_reached=False

    current_num_terms=[]
    best_num_terms=[]

    orbital_RMSEs=[]

    itct=0

    print_line(1)
    print(' Fitting Procedure:')
    print_line(1)

    if self.grasp_description.orbital_info.fitting_max_l==4:
      print('Fitted Orbitals: s,p,d,f,g')

    elif self.grasp_description.orbital_info.fitting_max_l==3:
      print('Fitted Orbitals: s,p,d,f')
      print('Number of Orbitals:'
            + " s:" + str(self.grasp_description.orbital_info.fitting_lvals[0])
            + " p:" + str(self.grasp_description.orbital_info.fitting_lvals[1])
            + " d:" + str(self.grasp_description.orbital_info.fitting_lvals[2])
            + " f:" + str(self.grasp_description.orbital_info.fitting_lvals[3]))

      while(not convergence_reached):
        for maximum_alpha in self.loop_max_alpha:

          itct=0
          for sp in range(self.max_num_terms,self.min_num_terms,-1):


            fit_loop_time=time.time()

            if itct>=5:
              continue

            alphas=get_alpha_list(self.min_alpha,maximum_alpha,sp,'geometric')
            A=get_A(sp, self.grasp_description.orbital_info.interpolated_radial_grid,alphas)
            AT=np.transpose(A)
            ATA=np.dot(AT,A)
            LU=pivotless_lu_decomp(ATA)
            B=np.transpose(np.matrix(copy.deepcopy(
                self.grasp_description.orbital_info.fitting_functions)))

            for d in range(sp,self.min_num_terms,-1):
              if itct>=5:
                break

              for f in range(d,self.min_num_terms,-1):
                if itct>=5:
                  break
                

                

                current_weight=7*f+5*d+4*sp

                if current_weight > best_weight:
                  continue

                current_num_terms=[sp,sp,d,f,0]

                coefficients,total_RMSE,current_RMSE,fit_loop_time = \
                        self.get_prefactor_coefficients(LU,A,AT,B,current_num_terms,alphas)


                orbital_RMSEs.append(current_RMSE)

                
                if total_RMSE < best_RMSE or abs(total_RMSE - best_RMSE) <= self.orb_fit_tolerance:
                  if round((maximum_alpha/max(self.loop_max_alpha))*100,2)%20==0:
                    print(str('%.2f'%((maximum_alpha/max(self.loop_max_alpha))*100))+
                          '%',maximum_alpha,max(self.loop_max_alpha),sp,d,f,'%.5E'%total_RMSE)
                  best_weight=current_weight
                  best_RMSE=total_RMSE
                  self.fitted_max_alpha=maximum_alpha
                  self.number_terms=current_num_terms
                  itct=0
                else:
                  itct+=1

                
                #'min_alpha','max_alpha','num_terms','num_s','num_p','num_d','num_f','num_g','weight','RMSE','time'



                self.fits.append([0.12,maximum_alpha,sp,sp,sp,d,f,0,
                                  current_weight,total_RMSE,current_RMSE,
                                  time.time()-fit_loop_time])


        orb_frame_RMSE=pd.DataFrame(data=orbital_RMSEs,
            columns=self.grasp_description.orbital_info.fitting_orbital_names)
        orb_frame_RMSE.to_csv('RMSE.csv',index=False)
        self.fitting_results=pd.DataFrame(data=self.fits,
                                          columns=self.fitting_columns)

        fig=px.box(self.fitting_results, x="max_alpha", y="total_RMSE",
                   color="weight",points="all")
        fig.show()


        print(self.fitting_results)
        convergence_reached=True
    elif self.grasp_description.orbital_info.fitting_max_l==2:
      print('Fitted Orbitals: s,p,d')
      print('Number of Orbitals:'
            + " s:" + str(self.grasp_description.orbital_info.fitting_lvals[0])
            + " p:" + str(self.grasp_description.orbital_info.fitting_lvals[1])
            + " d:" + str(self.grasp_description.orbital_info.fitting_lvals[2])
            + " f:" + str(self.grasp_description.orbital_info.fitting_lvals[3]))

      while(not convergence_reached):
        for maximum_alpha in self.loop_max_alpha:

          itct=0
          for sp in range(self.max_num_terms,self.min_num_terms,-1):


            fit_loop_time=time.time()

            if itct>=5:
              continue

            alphas=get_alpha_list(self.min_alpha,maximum_alpha,sp,'geometric')
            A=get_A(sp, self.grasp_description.orbital_info.interpolated_radial_grid,alphas)
            AT=np.transpose(A)
            ATA=np.dot(AT,A)
            LU=pivotless_lu_decomp(ATA)
            B=np.transpose(np.matrix(copy.deepcopy(
                self.grasp_description.orbital_info.fitting_functions)))

            for d in range(sp,self.min_num_terms,-1):
              if itct>=5:
                break

              current_weight = 5*d + 4*sp

              if current_weight > best_weight:
                continue

              current_num_terms=[sp,sp,d,0,0]

              coefficients,total_RMSE,current_RMSE, \
                      fit_loop_time=self.get_prefactor_coefficients(
                              LU,A,AT,B,current_num_terms,alphas)


              orbital_RMSEs.append(current_RMSE)

              
              if total_RMSE < best_RMSE or abs(total_RMSE - best_RMSE) <= self.orb_fit_tolerance:
                if round((maximum_alpha/max(self.loop_max_alpha))*100,2)%20==0:
                  print(str('%.2f'%((maximum_alpha/max(self.loop_max_alpha))*100))+'%',
                        maximum_alpha,max(self.loop_max_alpha),sp,d,'%.5E'%total_RMSE)
                best_weight=current_weight
                best_RMSE=total_RMSE
                self.fitted_max_alpha=maximum_alpha
                self.number_terms=current_num_terms
                itct=0
              else:
                itct+=1

              
              #'min_alpha','max_alpha','num_terms','num_s','num_p','num_d','num_f','num_g','weight','RMSE','time'



              self.fits.append([0.12,maximum_alpha,sp,sp,sp,d,0,0,
                                current_weight,total_RMSE,current_RMSE,
                                time.time()-fit_loop_time])


        orb_frame_RMSE=pd.DataFrame(data=orbital_RMSEs,
                columns=self.grasp_description.orbital_info.fitting_orbital_names)
        orb_frame_RMSE.to_csv('RMSE.csv',index=False)
        self.fitting_results=pd.DataFrame(data=self.fits,columns=self.fitting_columns)

        print ("coefficients")
        print (coefficients)
        print ("alphas")
        print (alphas)
        #fig=px.box(self.fitting_results, x="max_alpha", y="total_RMSE", color="weight",points="all")
        #fig.show()


        print(self.fitting_results)
        convergence_reached=True
    elif self.grasp_description.orbital_info.fitting_max_l==1 or \
            self.grasp_description.orbital_info.fitting_max_l==0:
      print('Fitted Orbitals: s,p')
      print('Number of Orbitals:'
            + " s:" + str(self.grasp_description.orbital_info.fitting_lvals[0])
            + " p:" + str(self.grasp_description.orbital_info.fitting_lvals[1])
            + " d:" + str(self.grasp_description.orbital_info.fitting_lvals[2])
            + " f:" + str(self.grasp_description.orbital_info.fitting_lvals[3]))

      #print ("self.loop_max_alpha", self.loop_max_alpha)
      while(not convergence_reached):
        for maximum_alpha in self.loop_max_alpha:

          #print("max_num_terms min_num_terms", self.max_num_terms, self.min_num_terms)
          itct=0
          for sp in range(self.max_num_terms,self.min_num_terms,-1):


            fit_loop_time=time.time()

            #print ("itct sp", itct, sp)
            if itct>=5:
              break

            alphas=get_alpha_list(self.min_alpha,maximum_alpha,sp,'geometric')
            A=get_A(sp, self.grasp_description.orbital_info.interpolated_radial_grid,alphas)
            AT=np.transpose(A)
            ATA=np.dot(AT,A)
            LU=pivotless_lu_decomp(ATA)
            B=np.transpose(np.matrix(copy.deepcopy(
                self.grasp_description.orbital_info.fitting_functions)))

            current_weight = 4*sp

            if current_weight > best_weight:
              continue

            current_num_terms=[sp,sp,0,0,0]

            coefficients,total_RMSE,current_RMSE, \
                    fit_loop_time = self.get_prefactor_coefficients(
                            LU,A,AT,B,current_num_terms,alphas)


            orbital_RMSEs.append(current_RMSE)

            #print("totRMSE bestRMSE fitTol", total_RMSE, best_RMSE, self.orb_fit_tolerance)
            
            if total_RMSE < best_RMSE or abs(total_RMSE - best_RMSE) <= self.orb_fit_tolerance:
              #print ("maxAlpha maxLoopAlpha", maximum_alpha, max(self.loop_max_alpha))
              if round((maximum_alpha/max(self.loop_max_alpha))*100,2)%20==0:
                print(str('%.2f'%((maximum_alpha/max(self.loop_max_alpha))*100))+
                      '%',maximum_alpha,max(self.loop_max_alpha),sp,'%.5E'%total_RMSE)
              best_weight=current_weight
              best_RMSE=total_RMSE
              self.fitted_max_alpha=maximum_alpha
              self.number_terms=current_num_terms
              itct=0
            else:
              itct+=1

           
            #'min_alpha','max_alpha','num_terms','num_s','num_p','num_d','num_f','num_g','weight','RMSE','time'



            self.fits.append([0.12,maximum_alpha,sp,sp,sp,0,0,0,
                              current_weight,total_RMSE,current_RMSE,
                              time.time()-fit_loop_time])


        orb_frame_RMSE=pd.DataFrame(data=orbital_RMSEs,
            columns=self.grasp_description.orbital_info.fitting_orbital_names)
        orb_frame_RMSE.to_csv('RMSE.csv',index=False)
        self.fitting_results=pd.DataFrame(data=self.fits,
                                          columns=self.fitting_columns)

        print ("coefficients")
        print (coefficients)
        print ("alphas")
        print (alphas)
        #fig=px.box(self.fitting_results, x="max_alpha", y="total_RMSE", color="weight",points="all")
        #fig.show()


        print(self.fitting_results)
        convergence_reached=True


    print(" Fitting Orbitals Complete:\t"+str('%.8E'%(time.time()-sweep_time))+"s. / "+str('%.8E'%(time.time()-start))+"s.")





#==========================
# PROGRAM OPERATIONS:      
#==========================

if __name__=="__main__":
  print_line(1)
  print_callout("\t\tI N T E R F I T    F I T T I N G    P R O G R A M")
  print_line(1)


  element_system=orbital_fitting()

  element_system.setup_calculation()

  element_system.fitting_procedure()
  
  '''
  # Set up the orbitals for calculation
  radialGrid,functionVals,orbList,lCounts,numLV=setUpCalc()

  # Perform parameter sweep:
  coeffVals,alphaList,RMSE,numTerms=doSweep(radialGrid,functionVals,orbList,lCounts,numLV)

  # Print out the values
  print ("Orbital\tl\tRMSE")
  for i in range(len(orbList)):
    print(str(orbList[i])+"\t"+str(lCounts[i])+"\t"+str("%.8E"%RMSE[i]))

  graphFunctions(orbList,lCounts,coeffVals,alphaList)
  writeFunctions(orbList,lCounts,coeffVals,alphaList)
  '''
 
  print_callout("Total Time: "+str('%.8E'%(time.time()-start))+"s.")
  print_line(1)
  print_line(1)
