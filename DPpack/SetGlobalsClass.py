import os, sys
import math
import shutil
import textwrap
import numpy as np

from DPpack.PTable import *
from DPpack.Misc import *

class Molecule:

	def __init__(self):

		self.atoms = []					# Lista de instancias de Atom
		self.positions = None			# Array Numpy
		self.energy = None				# Array Numpy
		self.gradient = None			# Array Numpy
		self.hessian = None				# Array Numpy

	def add_atom(self, a):

		self.atoms.append(a)			# Inserção de um novo atomo

	def center_of_mass(self):
	
		com = np.zeros(3)
		total_mass = 0.0
		for atom in self.atoms:
			total_mass += atom.mass
			com += atom.mass * np.array([atom.rx, atom.ry, atom.rz])
		
		com = com / total_mass
		
		return com


class Atom:

	def __init__(self, lbl,na,rx,ry,rz,chg,eps,sig):
		
		self.lbl = lbl                  # Integer
		self.na = na                    # Integer
		self.rx = rx                    # Double
		self.ry = ry                    # Double
		self.rz = rz                    # Double
		self.chg = chg                  # Double
		self.eps = eps                  # Double
		self.sig = sig                  # Double
		self.mass = atommass[self.na]   # Double
