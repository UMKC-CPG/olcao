#!/bin/env python

import os

def writeSkl(el,directory,n):
	cwd = os.getcwd()
	datadir = '../Downloads/lammps-28Jun14/Jobs/' + el + '-amorphous'
	os.chdir(datadir)
	
	inp = open('dump.melt')
	dump = inp.readlines()
	inp.close()

	rel = []
	for m in range(len(dump) - n, len(dump)):
		rel.append(dump[m])
	
	listout = []

	lattice = dump[len(dump) - n - 4].split()[1] + ' ' + dump[len(dump) - n - 3].split()[1] + ' ' + dump[len(dump) - n - 2].split()[1]
	lattice += ' 90.0 90.0 90.0\n'
	
	for atom in rel:
		coords = atom.split()
		line = el + ' ' + coords[2] + ' ' + coords[3] + ' ' + coords[4] + '\n'
		listout.append(line)

	os.chdir(cwd)
	os.chdir('../Olcaoskl')
	try:
		os.mkdir(directory)
	except:
		OSError
	os.chdir(directory)

	Skl = open('olcao.skl', 'w')
	Skl.write('title\n')
	title = directory + '\n'
	Skl.write(title)
	Skl.write('end\ncell\n')
	Skl.write(lattice)
	frac = 'fractional ' + str(n) + '\n'
	Skl.write(frac)
	for item in listout:
		Skl.write(item)
	Skl.write('space 1_a\nsupercell 1 1 1\nfull')
	Skl.close()
	
	filename = os.getcwd() + '/olcao.skl'

	os.chdir(cwd)
	return filename
