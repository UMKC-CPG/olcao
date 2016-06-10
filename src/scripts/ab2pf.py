#!/usr/bin/env python
################################################################
#
# ab2pf:
#   convert abinit output data into potfit reference configurations
#
#################################################################
# This program was written by Paul Rulis and was based on the 2013 version of
#   vasp2force.


import argparse
import gzip
import os
import sys
import subprocess

##########################
# DEFINE ALL SUBROUTINES #
##########################

# walk directory (recursively) and return all abinit output files
def get_abinitout_files(directory, recursive):
    sys.stderr.write('Searching directory %s for all abinit output files ...\n'
                    % directory)
    abinitOuts = []
    if not recursive:
        for item in os.listdir(directory):
            # Do a test to determine if the file is ASCII text or not.
            if "ASCII" in subprocess.check_output(["file",item]).split()[1:]:
                f = open(item,'r');
                line = f.readline();
                line = f.readline();
                f.close();
                if (line.count('Version') == 1 and line.count('ABINIT') == 1):
                    abinitOuts.append(os.path.join(directory,item))
    else:
        for root, SubFolders, files in os.walk(directory):
            for item in files:
                f = open(item,'r');
                line = f.readline();
                line = f.readline();
                f.close();
                if (line.count('Version') == 1 and line.count('ABINIT') == 1):
                    abinitOuts.append(os.path.join(directory,item))
    if len(abinitOuts) == 0:
        sys.stderr.write('Could not find any abinit output '
                         'files in this directory.\n')
    else:
        sys.stderr.write('Found the following files:\n')
        sys.stderr.write('  {}\n'.format(('\n  ').join(abinitOuts)))
        return abinitOuts
    return abinitOuts

# scans an abinit output file and returns a list with
# - number of configurations
# - atom types
def scan_abinitout_file(file_handle):
    # first try TOTAL-FORCE
    numAtoms = 0
    configs = 0
    atom_types = []
    ipt = []
    while 1:
        line = file_handle.readline();
        if not line: break
        if ('natom' in line) and ('mqgrid' in line):
            numAtoms = int(line.split()[5])
        if 'Cartesian forces' in line:
            configs += 1
        if 'atom type' in line:
            zNumber = int(line.split('/')[-1].split('.')[0].
               strip('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'))
            atom_types.append(elements[zNumber-1])
        if ' typat' in line:
            # Each integer in this line connects the type (element) of each
            #   atom to the order it appears in the rest of the file data.
            #   (E.g. coordinate list, forces list, etc.) At least one line is
            #   always read in at this point.
            ipt = [int(s) for s in line.split()[1:]]
            # If there are more than 20 atoms total then we need to read more
            #   lines. The exact number of lines to read depends on the number
            #   of atoms. 
            if (numAtoms%20 > 0) and (numAtoms>20):
                bonusLine = 1
            else:
                bonusLine = 0
            # Read a line for every group of 20 atoms or less (beyond those
            #   contained in the first line).
            for j in range(0,numAtoms/20 + bonusLine - 1):
                line = file_handle.readline();
                ipt += [int(s) for s in line.split()[0:]]
    if atom_types:
        return [configs, atom_types, ipt]
    else:
        sys.stderr.write('Could not determine atom types in file %s.\n'
                        % filename)
        sys.exit()


# return unique list without changing the order
# from http://stackoverflow.com/questions/480214
def uniq(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

##############################
# END DEFINE ALL SUBROUTINES #
##############################


###########################
# BEGIN PROGRAM EXECUTION #
###########################

# Establish program conversion factors.
bohr=0.529177208 # 1 bohr = this many Angstroms
eVPerA3 = 160.2176487 # 1 eV/A^3 = this many GPa
HaPerBohr3 = 29421.010901602753 # 1 Ha/Bohr^3 = this many GPa
pressureRatio = HaPerBohr3 / eVPerA3 # Multiplier to convert Ha/Bohr^3 to eV/A^3
HaPerBohr = 51.42208619083232 # 1 Ha/Bohr = this many eV/A
Ha = 27.21138386 # 1 Ha = this many eV.


# Establish a periodic table of the elements list. Note that the atomic Z
#   number minus one corresponds to the index of the element in the list.
elements = ['H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', 'Na',
            'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' , 'Ca', 'Sc', 'Ti',
            'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ga', 'As',
            'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
            'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
            'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir',
            'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra',
            'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
            'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
            'Rg', 'Cn']

# Establish the possible command line parameters.
parser = argparse.ArgumentParser(
#    formatter_class=argparse.RawDescriptionHelpFormatter,
    description = '''Converts abinit output data into potfit '''
                  '''reference configurations.''')

parser.add_argument('-c', type=str, required=False,
                help='list of chemical species to use, e.g. -c Mg=0,Zn=1')
parser.add_argument('-e', type=str, required=False,
                help='file with single atom energies (NYI)')
parser.add_argument('-f', '--final', action='store_true',
                help='use only the final configuration from abinit output')
parser.add_argument('-l', '--list', action='store_true',
                help='list all chemical species and exit')
parser.add_argument('-r', '--recursive', action='store_true',
                help='scan recursively for abinit output files')
parser.add_argument('-s', '--configs', type=str,
                help='comma separated list of configurations to use from each '
                'file \n (use the same configurations for all abinit output '
                'files)')
parser.add_argument('-w', '--weight', type=float, default=1.0,
                help='set configuration weight for all configurations')
parser.add_argument('files', type=str, nargs='*',
                help='list of abinit output files (plain or gzipped)')

args = parser.parse_args()

# Check for sane arguments
if args.weight < 0:
    sys.stderr.write("The weight needs to be positive!\n")
    sys.exit()

# Determine all abinit output files
abinitOuts = []
# If we didn't ask for specific files on the command line, then go search for
#   them in the current directory (recursively if requested on the command
#   line).
if not args.files:
    abinitOuts = get_abinitout_files('.',args.recursive)
# Use the list of files and directories given on the command line.
for item in args.files:
    if os.path.isdir(item): # Recursively search directories
        abinitOuts += get_abinitout_files(item, args.recursive)
    else: # Directly append files.
        abinitOuts.append(os.path.abspath(item))

# remove duplicate entries from the abinitOuts table
abinitOuts = uniq(abinitOuts)

# read metadata from all abinit output files
data= []
max_types = 1
for item in abinitOuts:
    if item.endswith('.gz'):
        f = gzip.open(item, 'rb')
    else:
        f = open(item, 'r')
    data.append(scan_abinitout_file(f))
    f.close()
    # Look at the the second element *of* the last list element. This is the
    #   list of elements names in this configuration. Count the number of
    #   names and this becomes the possible max number of types.
    max_types = max(max_types, len(data[-1][1]))

# print list of all files (and meta-data (# of types (elements) and # of
#   configurations) about the files) and exit
if args.list:
    print ""
    for i in range(len(abinitOuts)):
        print ("Found {} atom types in {} configurations in file {}:"
                        .format(len(data[i][1]), data[i][0], abinitOuts[i]))
        for j in range(len(data[i][1])):
            sys.stdout.write(" {}={}".format(data[i][1][j],j))
        print ""
    sys.exit()

# TODO: account for SAE file
if args.e:
    sys.stderr.write("\nSingle atom energies are (not yet) supported\n")
    sys.exit()

# check the config string for appropriate and valid values
configs = []
if args.configs:
    for item in args.configs.split(','):
        if not item:
            sys.stderr.write("ERROR: Could not read the -s string.\n")
            sys.exit()
        if len(item.split('-')) == 1:
            try:
                configs.append(int(item))
            except:
                sys.stderr.write("ERROR: Could not read the -s string.\n")
                sys.exit()
        else:
            if len(item.split('-')) != 2:
                sys.stderr.write("ERROR: Could not read the -s string.\n")
                sys.exit()
            else:
                items = item.split('-')
                try:
                    val_min = int(items[0])
                    val_max = int(items[1])
                except:
                    sys.stderr.write("ERROR: Could not read the -s string.\n")
                    sys.exit()
                if val_min > val_max:
                    val_min = val_max
                    val_max = int(items[0])
                for i in range(val_min,val_max+1):
                    configs.append(i)

# abort with an error message if no configurations are specified
if (len(configs) == 0 and not args.final):
    sys.stderr.write("ERROR: \tNo configurations specified.\n")
    sys.stderr.write("\tUse the -s and/or -f command line switches.\n")
    sys.exit()

# check for valid configurations
for i in range(len(data)):
    item = data[i]
    if configs and int(configs[-1]) > int(item[0]):
        sys.stderr.write('\nThere are only %s configurations in %s.\n'
                        % (item[0],abinitOuts[i]))
        sys.stderr.write('Please adjust your -s string!\n')
        sys.exit()

# check the types string
types = dict()
numbers = dict()
if args.c:
    if len(args.c.split(',')) > max_types:
        sys.stderr.write("\nERROR: There are too many items "
                        "in your -c string!\n")
        sys.exit()
    if len(args.c.split(',')) < max_types:
        sys.stderr.write("\nERROR: There are not enough items "
                        "in your -c string!\n")
        sys.exit()
    for item in args.c.split(','):
        if len(item.split('=')) != 2:
            sys.stderr.write("\nERROR: Could not read the -c string.\n")
            sys.stderr.write("Maybe a missing or extra '=' sign?\n")
            sys.exit()
        else:
            try:
                name = str(item.split('=')[0])
                number = int(item.split('=')[1])
            except:
                sys.stderr.write("\nERROR: Could not read the -c string\n")
                sys.exit()
            if number >= max_types:
                sys.stderr.write("\nERROR: The atom type for %s is invalid!\n"
                                % name)
                sys.exit()
            if name in types:
                sys.stderr.write("\nERROR: Duplicate atom type found in -c "
                                "string\n")
                sys.exit()
            if number in numbers:
                sys.stderr.write("\nERROR: Duplicate atom number found in -c "
                                "string\n")
                sys.exit()
            types[name] = number
            numbers[number] = name

# Now we finally get down to the business of making the configuration file.
for i in range(len(abinitOuts)):
    if abinitOuts[i].endswith('.gz'):
        f = gzip.open(abinitOuts[i], 'rb')
    else:
        f = open(abinitOuts[i], 'r')
    natoms = len(data[i][2])
    count = 0
    line = f.readline()
    box_x = []
    box_y = []
    box_z = []
    while line != '':
        line = f.readline()
        if 'STEP' in line: # Note that we found the next ion move iteration.
            energy = 0
            stress = [None]*6
        if 'ETOT' in line: # Get the total energy. (Note that we read every
                # one, but only the last one is actually kept because this
                # repeatedly overwrites the previous values.)
            energy = float(line.split()[2]) / natoms * Ha
        if 'space primitive' in line:
            # Read the lattice parameters and convert Bohr to A.
            line = f.readline()
            box_x = [float(s)*bohr for s in
                            line.replace('-',' -').split()[1:4]]
            line = f.readline()
            box_y = [float(s)*bohr for s in
                            line.replace('-',' -').split()[1:4]]
            line = f.readline()
            box_z = [float(s)*bohr for s in
                            line.replace('-',' -').split()[1:4]]
        if 'stress tensor (hartree/bohr^3)' in line:
            # Convert Ha/Bohr^3 to eV/A^3
            line = f.readline()
            stress[0] = float(line.split()[2])*pressureRatio # xx
            stress[4] = float(line.split()[5])*pressureRatio # yz
            line = f.readline()
            stress[1] = float(line.split()[2])*pressureRatio # yy
            stress[5] = float(line.split()[5])*pressureRatio # xz
            line = f.readline()
            stress[2] = float(line.split()[2])*pressureRatio # zz
            stress[3] = float(line.split()[5])*pressureRatio # xy
        if 'Cartesian coordinates' in line:
            cart = []
            force = []
            count += 1
            # Run a check to be sure that we want this configuration.
            if count in configs or (args.final and count == data[i][0]):
                print "#N %s 1" % natoms
                sys.stdout.write("#C")
                # Record the element name for this atom from either the args.c
                #   list or from the automatically collected data.
                if args.c:
                    sys.stdout.write(" %s" % numbers[0])
                    for j in range(1,max_types):
                        sys.stdout.write(' %s' % numbers[j])
                else:
                    sys.stdout.write(" %s" % data[i][1][0])
                    for j in range(1,max_types):
                        sys.stdout.write(' %s' % data[i][1][j])
                sys.stdout.write("\n")
                print ("## force file generated from file %s" % abinitOuts[i])
                print ("#X {:13.8f} {:13.8f} {:13.8f}"
                                .format(box_x[0],box_x[1],box_x[2]))
                print ("#Y {:13.8f} {:13.8f} {:13.8f}"
                                .format(box_y[0],box_y[1],box_y[2]))
                print ("#Z {:13.8f} {:13.8f} {:13.8f}"
                                .format(box_z[0],box_z[1],box_z[2]))
                print ("#W {:f}".format(args.weight))
                print ("#E {:.10f}".format(energy))
                if stress:
                    sys.stdout.write("#S")
                    for num in range(6):
                        sys.stdout.write(' {:8.7g}'.format(stress[num]))
                    sys.stdout.write('\n')
                print "#F"
                # For the current i configuration, read the coordinates of each
                #   atom, then read the force for each atom, then print the
                #   data for each atom. (One natoms loop after another.)
                for num in range(natoms):
                    line = [float(s)*bohr for s in f.readline().split()]
                    cart.append(line)
                f.readline() # Read the header line for forces.
                for num in range(natoms):
                    line = [float(s)*HaPerBohr for s in f.readline().split()]
                    force.append(line)
                for num in range(natoms):
                    if args.c:
                        sys.stdout.write('{} '.
                                 format(types[data[i][1][data[i][2][num]-1]]))
                    else:
                        sys.stdout.write('{} '.
                                 format(data[i][2][num]-1))
                    sys.stdout.write("{:12.7g} {:12.7g} {:12.7g}".
                              format(cart[num][0],cart[num][1],cart[num][2]))
                    sys.stdout.write("{:>12.7g} {:>12.7g} {:>12.7g}\n".
                              format(force[num][0],force[num][1],force[num][2]))
    f.close()
