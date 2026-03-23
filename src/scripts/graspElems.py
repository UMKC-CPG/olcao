#!/usr/bin/env python3

# PROGRAM: graspElems
# AUTHOR:  Patrick Ryan Thomas  (prtnpb@mail.umkc.edu)
# CONTACT: Paul Rulis (rulisp@umkc.edu)
#
# Runs the Grasp2K pipeline for specified elements and basis sets,
# or generates SBATCH scripts for cluster submission (--batch mode).
#
# Usage:
#   Direct execution (default):
#     graspElems.py C MB
#     graspElems.py total FB
#     graspElems.py F MB H EB
#     graspElems.py C all
#
#   Batch/script generation:
#     graspElems.py --batch C MB
#     graspElems.py --batch -nc 16 total

import argparse as ap
import os
import sys
import subprocess
import re
from datetime import datetime

from element_data import ElementData

_ed = ElementData()
_ed.init_element_data()


# --------------------------------------------------------------------------- #
#                           Script Settings                                    #
# --------------------------------------------------------------------------- #

class ScriptSettings():
    """The instance variables of this object are the user settings that
       control the program. The variable values are pulled from a list
       that is created within a resource control file and that are then
       reconciled with command line parameters."""


    def __init__(self):
        """Define default values for the parameters by pulling them
        from the resource control file in the default location:
        $OLCAO_RC/graspElemsrc.py or from the current working directory
        if a local copy of graspElemsrc.py is present."""

        # Read default variables from the resource control file.
        rc_dir = os.getenv('OLCAO_RC')
        if not rc_dir:
            sys.exit("Error: $OLCAO_RC is not set. See instructions.")
        sys.path.insert(1, rc_dir)
        from graspElemsrc import parameters_and_defaults
        default_rc = parameters_and_defaults()

        # Assign values to the settings from the rc defaults file.
        self.assign_rc_defaults(default_rc)

        # Parse the command line.
        args = self.parse_command_line()

        # Reconcile the command line arguments with the rc file.
        self.reconcile(args)

        # Record the command line parameters that were used.
        self.recordCLP()


    def assign_rc_defaults(self, default_rc):

        # Execution mode.
        self.batch = default_rc["batch"]

        # SBATCH grouping.
        self.num_cores = default_rc["num_cores"]


    def parse_command_line(self):

        # Create the parser tool.
        prog_name = "graspElems"

        description_text = """
Run the Grasp2K atomic structure pipeline (rnucleus, csl, rangular,
rwfnestimate, rmcdhf, readrwf) for specified elements and basis sets.

By default, programs are executed directly via subprocess. Use --batch
to generate runC/runRSCF bash scripts and SBATCH job files instead.

Element specifications:
  ELEM BASIS      Run one element with one basis (e.g. C MB)
  ELEM all        Run one element with all bases (MB, FB, EB)
  total           Run all elements with all bases
  total BASIS     Run all elements with one basis
"""

        epilog_text = """
Originally written by Patrick Ryan Thomas (prtnpb@mail.umkc.edu)
Please contact Paul Rulis (rulisp@umkc.edu) regarding questions.
Defaults are given in ./graspElemsrc.py or $OLCAO_RC/graspElemsrc.py.
"""

        parser = ap.ArgumentParser(prog = prog_name,
                formatter_class=ap.RawDescriptionHelpFormatter,
                description = description_text,
                epilog = epilog_text)

        # Add arguments to the parser.
        self.add_parser_arguments(parser)

        # Parse the arguments and return the results.
        return parser.parse_args()


    def add_parser_arguments(self, parser):

        # Define the batch argument.
        parser.add_argument('-b', '--batch', dest='batch',
                            action='store_true', default=self.batch,
                            help='Generate SBATCH scripts instead of running '
                                 f'directly. Default: {self.batch}')

        # Define the num_cores argument.
        parser.add_argument('-nc', '--num-cores', dest='num_cores',
                            type=int, default=self.num_cores,
                            help='Number of cores for SBATCH job grouping '
                                 f'(batch mode only). Default: {self.num_cores}')

        # Define the positional arguments (element/basis specifications).
        parser.add_argument('specs', nargs='+', metavar='SPEC',
                            help='Element/basis specifications: '
                                 '"ELEM BASIS", "ELEM all", "total", '
                                 'or "total BASIS"')


    def reconcile(self, args):
        self.batch = args.batch
        self.num_cores = args.num_cores
        self.comm_list = build_comm_list(args.specs)


    def recordCLP(self):
        with open("command", "a") as cmd:
            now = datetime.now()
            formatted_dt = now.strftime("%b. %d, %Y: %H:%M:%S")
            cmd.write(f"Date: {formatted_dt}\n")
            cmd.write(f"Cmnd:")
            for argument in sys.argv:
                cmd.write(f" {argument}")
            cmd.write("\n\n")


# --------------------------------------------------------------------------- #
#                         Atom / Element Lists                                 #
# --------------------------------------------------------------------------- #

def getAtomList():
    l=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si',
       'P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',
       'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb',
       'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
       'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
       'Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg',
       'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np',
       'Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
    return l


def build_comm_list(specs):
    """Convert the positional specs list into element/basis pairs.

    Accepts the same specification conventions as the original parseInputs:
      ['total']            -> all elements, all bases
      ['total', 'MB']      -> all elements, one basis
      ['C', 'MB']          -> one element, one basis
      ['C', 'all']         -> one element, all bases
      ['F', 'MB', 'H', 'EB'] -> multiple element/basis pairs
    """
    elem_list = getAtomList()
    all_bases = ['MB', 'FB', 'EB']
    tmparr = []

    i = 0
    while i < len(specs):
        token = specs[i]

        if token == 'total':
            # Check if a basis follows.
            if i + 1 < len(specs) and specs[i+1] in all_bases:
                basis = specs[i+1]
                for elem in elem_list:
                    tmparr.extend([elem, basis])
                i += 2
            else:
                for elem in elem_list:
                    for basis in all_bases:
                        tmparr.extend([elem, basis])
                i += 1

        elif token in elem_list:
            if i + 1 >= len(specs):
                sys.exit(f"Error: element '{token}' given without a basis.")
            basis = specs[i+1]
            if basis == 'all':
                for b in all_bases:
                    tmparr.extend([token, b])
            else:
                tmparr.extend([token, basis])
            i += 2

        else:
            # Skip tokens that are neither 'total' nor an element
            # (e.g. a basis that was already consumed).
            i += 1

    return tmparr


# --------------------------------------------------------------------------- #
#                         Isotope / Element Data                               #
# --------------------------------------------------------------------------- #

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
        filename = os.getenv('OLCAO_DIR') + "/share/isotopes.dat"
        atomInfo=getFileContents(filename)

        n=0
        for i in range(59,len(atomInfo),18):
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


# --------------------------------------------------------------------------- #
#                       Grasp2K Execution Helpers                              #
# --------------------------------------------------------------------------- #

def count_blocks(rcsf_path):
    """Count unique J/parity blocks in an rcsf.inp file."""
    with open(rcsf_path) as f:
        content = f.read()
    blocks = set(re.findall(r'\b(\d+[+-])\s*$', content, re.MULTILINE))
    return len(blocks)


def run_grasp(program, stdin_text, cwd):
    """Run a Grasp2K program with the given stdin text."""
    grasp_bin = os.environ.get('GRASP_BIN',
                               os.path.join(os.environ['GRASP_DIR'], 'bin'))
    exe = os.path.join(grasp_bin, program)
    print(f"  Running {program}...", flush=True)
    result = subprocess.run(
        [exe],
        input=stdin_text, text=True,
        capture_output=True, cwd=cwd
    )
    if result.returncode != 0:
        print(f"ERROR running {program} in {cwd}:\n{result.stderr}",
              file=sys.stderr)
        sys.exit(1)
    return result


def get_val_atoms(atomDat, basis, atNum):
    """Return the valence atom string for the given basis."""
    if basis == "MB":
        return atomDat.mbValAtoms[atNum-1]
    elif basis == "FB":
        return atomDat.fbValAtoms[atNum-1]
    elif basis == "EB":
        return atomDat.ebValAtoms[atNum-1]


def get_num_subs(atomDat, basis, atNum):
    """Return the numSubs string for the given basis."""
    if basis == "MB":
        return atomDat.mbNumSubs[atNum-1]
    elif basis == "FB":
        return atomDat.fbNumSubs[atNum-1]
    elif basis == "EB":
        return atomDat.ebNumSubs[atNum-1]


# --------------------------------------------------------------------------- #
#                        Direct Execution Mode                                 #
# --------------------------------------------------------------------------- #

def run_element(element, basis, atomDat):
    """Run the full Grasp2K pipeline for one element/basis pair."""
    atNum = _ed.get_element_z(element)
    orig_dir = os.getcwd()

    # Create directory structure.
    os.makedirs(os.path.join(basis, element), exist_ok=True)
    work_dir = os.path.join(orig_dir, basis, element)

    print(f"Processing {element} {basis} in {work_dir}")

    # 1) rnucleus
    isotope_str = atomDat.isotopes[atNum-1].replace(element, '')
    rnucleus_stdin = (
        f"{atNum}\n"
        f"{isotope_str}\n"
        f"n\n"
        f"{atomDat.atomicMass[atNum-1]}\n"
        f"{atomDat.spin[atNum-1]}\n"
        f"{atomDat.magMoment[atNum-1]}\n"
        f"{atomDat.quadMoment[atNum-1]}\n"
        f"n\n"
    )
    print (rnucleus_stdin)
    run_grasp('rnucleus', rnucleus_stdin, work_dir)

    # 2) csl
    val_atoms = get_val_atoms(atomDat, basis, atNum)
    num_subs = get_num_subs(atomDat, basis, atNum)
    csl_stdin = (
        f"y\n"
        f"{atomDat.coreAtoms[atNum-1]}\n"
        f"{val_atoms}\n"
        f"n\nn\n"
        f"{atomDat.occAtoms[atNum-1]}\n"
        f"\n\n"
        f"{num_subs}\n"
    )
    run_grasp('csl', csl_stdin, work_dir)

    # 3) mv rcsl.out -> rcsf.inp
    os.rename(os.path.join(work_dir, 'rcsl.out'),
              os.path.join(work_dir, 'rcsf.inp'))

    # 4) rangular
    run_grasp('rangular', "y\n", work_dir)

    # 5) rwfnestimate
    run_grasp('rwfnestimate', "y\n2\n*\n", work_dir)

    # 6) Count J/parity blocks from rcsf.inp
    n_blocks = count_blocks(os.path.join(work_dir, 'rcsf.inp'))
    print(f"  Found {n_blocks} J/parity block(s)")

    # 7) rmcdhf
    rmcdhf_stdin = "y\ny\n" + "1\n" * n_blocks + "*\n\n100000\n"
    run_grasp('rmcdhf', rmcdhf_stdin, work_dir)

    # 8) Rename wave function files.
    os.rename(os.path.join(work_dir, 'rwfn.inp'),
              os.path.join(work_dir, 'rwfn-orig.inp'))
    os.rename(os.path.join(work_dir, 'rwfn.out'),
              os.path.join(work_dir, 'rwfn.inp'))

    # 9) readrwf
    run_grasp('readrwf', "1\nrwfn.inp\nrwfn.out\n", work_dir)

    print(f"  Done: {element} {basis}")


def run_pipeline(commList):
    """Run the Grasp2K pipeline directly for each element/basis pair."""
    elem_list = getAtomList()
    atomDat = elemData()
    atomDat.getElemData()

    for i in range(len(commList)):
        if commList[i] in elem_list:
            run_element(commList[i], commList[i+1], atomDat)


# --------------------------------------------------------------------------- #
#                        Batch / Script Generation Mode                        #
# --------------------------------------------------------------------------- #

def make_scripts(commList, num_cores):
    """Generate runC/runRSCF bash scripts and SBATCH job files."""
    elem_list=getAtomList()
    atomDat=elemData()
    atomDat.getElemData()
    execLocs=[]
    batchHead=[]
    batchFileNum=1
    elemCount=1

    if not os.path.exists("SBATCHhead"):
        f=open("SBATCHhead",'w')
        temp = """#!/bin/bash
#SBATCH -J myJob
#SBATCH -o myJob
#SBATCH -n 16
#SBATCH -p normal
#SBATCH -t 30:00:00
#SBATCH --mail-user=rulisp@umkc.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

source "$GRASP_DIR/make-environment_ifort"

cd "$WORK/ELEM_FILES/graspElems"
"""
        f.write(temp)
        f.close()

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
              if elemCount%num_cores==0 and elemCount !=1:
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

          atNum=_ed.get_element_z(commList[i])

          #ISOTOPE
          g=open("run"+commList[i],'w')
          g.write(batchHead[0])
          g.write("set -x\n")
          g.write("\n$GRASP_BIN/rnucleus <<S1\n")
          g.write(str(atNum)+"\n")
          g.write((atomDat.isotopes[atNum-1]).replace(commList[i],''))
          g.write('\nn\n')
          g.write(atomDat.atomicMass[atNum-1]+'\n')
          g.write(atomDat.spin[atNum-1]+'\n')
          g.write(atomDat.magMoment[atNum-1]+'\n')
          g.write(atomDat.quadMoment[atNum-1]+'\n')
          g.write('n\n')
          g.write('S1\n\n')

          #CSL
          g.write("$GRASP_BIN/csl <<S1\n")
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

          g.write("mv rcsl.out rcsf.inp\n\n")

          # RANGULAR
          g.write("$GRASP_BIN/rangular <<S1\n")
          g.write("y\n")
          g.write("S1\n\n")

          # RWFNESTIMATE
          g.write("$GRASP_BIN/rwfnestimate <<S1\n")
          g.write("y\n")
          g.write("2\n")
          g.write("*\n")
          g.write("S1\n\n")

          g.close()
          comm="chmod +x run"+commList[i]
          os.system(comm)

          g=open("runRSCF",'w')

          g.write(batchHead[0])
          g.write("set -x\n")

          g.write("$GRASP_BIN/rmcdhf <<S1\n")
          g.write("y\ny\n1\n*\n\n100000\nS1\n\n")

          g.write("mv rwfn.inp rwfn-orig.inp\n")
          g.write("mv rwfn.out rwfn.inp\n\n")

          g.write("$GRASP_BIN/readrwf <<S1\n")
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


# --------------------------------------------------------------------------- #
#                              Main Program                                    #
# --------------------------------------------------------------------------- #

def main():

    # Get script settings from a combination of the resource control file
    #   and parameters given by the user on the command line.
    settings = ScriptSettings()

    # Execute based on mode.
    if settings.batch:
        make_scripts(settings.comm_list, settings.num_cores)
    else:
        run_pipeline(settings.comm_list)


if __name__ == '__main__':
    # Everything before this point was a subroutine definition or a request
    #   to import information from external modules. Only now do we actually
    #   start running the program. The purpose of this is to allow another
    #   python program to import *this* script and call its functions
    #   internally.
    main()
