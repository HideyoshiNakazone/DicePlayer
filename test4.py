#!/opt/local/bin/python3

import sys, math
from math import pi
import numpy as np
from numpy import linalg
from numpy.random import rand, random
from copy import deepcopy

from DPpack.MolHandling import *
from DPpack.PTable import *


natoms = input("# of atoms: ")
natoms = int(natoms)
file = input("xyz file: ")

molecules=[]

with open(file) as fh:
	xyzfile = fh.readlines()

file += ".new"
file2 = file + "2"
fh = open(file, "w")
fh2 = open(file2, "w")

total_atoms = int(xyzfile.pop(0).split()[0])
fh.write("{}\n".format(total_atoms))
comment = xyzfile.pop(0)
fh.write("{}".format(comment))

nmols = round(total_atoms/natoms)
for i in range(nmols):
	molecules.append([])
	for j in range(natoms):
		molecules[i].append({})
		line = xyzfile.pop(0).split()
		molecules[i][j]['na'] = int(line[0])
		molecules[i][j]['rx'] = float(line[1])
		molecules[i][j]['ry'] = float(line[2])
		molecules[i][j]['rz'] = float(line[3])
		molecules[i][j]['mass'] = atommass[molecules[0][j]['na']]

for atom in molecules[0]:
	fh.write("{:>4s} {:>11.5f} {:>11.5f} {:>11.5f}\n".format(atomsymb[atom['na']], 
													  atom['rx'], atom['ry'], atom['rz']))
	fh2.write("{:>4s} {:>11.5f} {:>11.5f} {:>11.5f}\n".format(atomsymb[atom['na']], 
													  atom['rx'], atom['ry'], atom['rz']))

for i in range(1, nmols):
	fh2.write("{}\n".format(natoms))
	fh2.write("{}".format(comment))
	rmsd, projected_mol = rmsd_fit(molecules[i], molecules[0])
	print("{:>9.5f}".format(rmsd))
	for atom in projected_mol:
		fh.write("{:>4s} {:>11.5f} {:>11.5f} {:>11.5f}\n".format(atomsymb[atom['na']], 
													  atom['rx'], atom['ry'], atom['rz']))
		fh2.write("{:>4s} {:>11.5f} {:>11.5f} {:>11.5f}\n".format(atomsymb[atom['na']], 
													  atom['rx'], atom['ry'], atom['rz']))
	
fh.close()
fh2.close()

