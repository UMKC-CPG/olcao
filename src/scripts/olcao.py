#!/bin/env python

import os
import sys
from olcaoObject import olcaoJob

############# Subroutines to manage the runtime #########################
def printHelp():
  print "Help Stuff"

# Main Subroutine
def main():
  # Get command line parameters and remove the first entry in the list
  # because the first entry is the name of this script
  cmdParams = sys.argv
  cmdParams.pop(0)

  dohelp = cmdParams.count('-help')
  if (dohelp > 0):
    printHelp()
    sys.exit()

  # Initialize the OLCAO Data Structure
  jobdir = os.getcwd()
  print cmdParams
  ojob = olcaoJob(jobdir,cmdParams)

  ojob.run()


# If not imported as a module run the main subroutine 
if __name__ == "__main__":
  main()
