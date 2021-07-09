from DPpack.MolHandling import total_mass
import os, sys
import math
import shutil
import textwrap
import sys, math
from copy import deepcopy

import numpy as np
from numpy import linalg

from DPpack.Misc import *
from DPpack.PTable import *
from DPpack.SetGlobals import *

# Usaremos uma nova classe que ira conter toda interação entre moleculas

class System:

	def __init__(self):

		self.molecule = []

	def add_molecule(self, m):

		self.molecule.append(m)

	# Função que calcula a distância entre dois centros de massa
	# e por se tratar de uma função de dois atomos não deve ser
	# inserida dentro de Molecule
	def center_of_mass_distance(self, a, b):
	
		com1 = self.molecule[a].center_of_mass()
		com2 = self.molecule[b].center_of_mass()
		dx = com1[0] - com2[0]
		dy = com1[1] - com2[1]
		dz = com1[2] - com2[2]
		distance = math.sqrt(dx**2 + dy**2 + dz**2)
		
		return distance

	def minimum_distance(self, index1, index2):
	
		distances = []
		for atom1 in self.molecule[index1]:
			if atom1.na != ghost_number:
				for atom2 in self.molecule[index2]:
					if atom2.na != ghost_number:
						dx = atom1.rx - atom2.rx
						dy = atom1.ry - atom2.ry
						dz = atom1.rz - atom2.rz
						distances.append(math.sqrt(dx**2 + dy**2 + dz**2))
		
		return min(distances)

# Classe que conterá toda informação e funções relacionadas a uma unica molecula

class Molecule:

	def __init__(self):

		self.atom = []					# Lista de instancias de Atom
		self.position = None			# Array Numpy
		self.energy = None				# Array Numpy
		self.gradient = None			# Array Numpy
		self.hessian = None				# Array Numpy
		self.total_mass = 0

	def add_atom(self, a):

		self.atom.append(a)			# Inserção de um novo atomo
		self.total_mass += a.mass

	def center_of_mass(self):
	
		com = np.zeros(3)
		total_mass = 0.0

		for atom in self.atom:
			
			total_mass += atom.mass
			com += atom.mass * np.array([atom.rx, atom.ry, atom.rz])
		
		com = com / total_mass
		
		return com

	def center_of_mass_to_origin(self):
	
		com = self.center_of_mass()

		for atom in self.atom:

			atom.rx -= com[0]
			atom.ry -= com[1]
			atom.rz -= com[2]
	
	def charges_and_dipole(self):
	
		eA_to_Debye = 1/0.20819434
		charge = 0
		dipole = np.zeros(3)
		for atom in self.atom:
			position = np.array([ atom.rx, atom.ry, atom.rz ])
			dipole += atom.chg * position
			charge += atom.chg
		
		dipole *= eA_to_Debye
		total_dipole = math.sqrt(dipole[0]**2 + dipole[1]**2 + dipole[2]**2)
		
		return [charge, dipole[0], dipole[1], dipole[2], total_dipole]

	def distances_between_atoms(self):
		
		distances = []
		dim = len(self.atom)
		for atom1 in self.atom:
			if atom1.na != ghost_number:
				for atom2 in self.atom:
					if atom2.na != ghost_number:
						dx = atom1.rx - atom2.rx
						dy = atom1.ry - atom2.ry
						dz = atom1.rz - atom2.rz
						distances.append(math.sqrt(dx**2 + dy**2 + dz**2))
		
		return np.array(distances).reshape(dim, dim)

	def eixos(self):
	
		eixos = np.zeros(3)
		if len(self.atom) == 2:

			position1 = np.array([ self.atom[0].rx, self.atom[0].ry, self.atom[0].rz ])
			position2 = np.array([ self.atom[1].rx, self.atom[1].ry, self.atom[1].rz ])
			eixos = position2 - position1
			eixos /= linalg.norm(eixos)
		
		elif len(self.atom) > 2:
		
			position1 = np.array([ self.atom[0].rx, self.atom[0].ry, self.atom[0].rz ])
			position2 = np.array([ self.atom[1].rx, self.atom[1].ry, self.atom[1].rz ])
			position3 = np.array([ self.atom[2].rx, self.atom[2].ry, self.atom[2].rz ])
			v1 = position2 - position1
			v2 = position3 - position1
			v3 = np.cross(v1, v2)
			v2 = np.cross(v1, v3)
			v1 /= linalg.norm(v1)
			v2 /= linalg.norm(v2)
			v3 /= linalg.norm(v3)
			eixos = np.array([[v1[0], v1[1], v1[2]],
							[v2[0], v2[1], v2[2]],
							[v3[0], v3[1], v3[2]]])
		
		return eixos

	def inertia_tensor(self):
	
		com = self.center_of_mass()
		Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0.0

		for atom in self.atom:

			####  Obtain the displacement from the center of mass
			dx = atom.rx - com[0]
			dy = atom.ry - com[1]
			dz = atom.rz - com[2]
			####  Update the diagonal components of the tensor
			Ixx += atom.mass * (dy**2 + dz**2)
			Iyy += atom.mass * (dz**2 + dx**2)
			Izz += atom.mass * (dx**2 + dy**2)
			####  Update the off-diagonal components of the tensor
			Ixy += atom.mass * dx * dy * -1
			Ixz += atom.mass * dx * dz * -1
			Iyz += atom.mass * dy * dz * -1
			
		return np.array([	[Ixx, Ixy, Ixz],
							[Ixy, Iyy, Iyz],
							[Ixz, Iyz, Izz]	])

	def principal_axes(self):
		
		try:
			evals, evecs = linalg.eigh(self.inertia_tensor())
		except:
			sys.exit("Error: diagonalization of inertia tensor did not converge")
		
		return evals, evecs

	def read_position(self):
	
		position_list = []
		for atom in self.atom:
			position_list.extend([ atom.rx, atom.ry, atom.rz ])
		position = np.array(position_list)
		position *= ang2bohr
		
		return position

	def update_hessian(self, step, cur_gradient):	## According to the BFGS
		
		dif_gradient = cur_gradient - self.gradient
		
		mat1 = 1/np.dot(dif_gradient, step) * np.matmul(dif_gradient.T, dif_gradient)
		mat2 = 1/np.dot(step, np.matmul(self.hessian, step.T).T)
		mat2 *= np.matmul( np.matmul(self.hessian, step.T), np.matmul(step, hessian) )
		
		self.hessian += mat1 - mat2

	def sizes_of_molecule(self):
	
		x_list = []
		y_list = []
		z_list = []

		for atom in self.atom:
			if atom.na != ghost_number:
				x_list.append(atom.rx)
				y_list.append(atom.ry)
				z_list.append(atom.rz)
		
		x_max = max(x_list)
		x_min = min(x_list)
		y_max = max(y_list)
		y_min = min(y_list)
		z_max = max(z_list)
		z_min = min(z_list)
		
		sizes = [x_max - x_min, y_max - y_min, z_max - z_min]
		
		return sizes

	def standard_orientation(self):
		
		self.center_of_mass_to_origin()
		tensor = self.inertia_tensor()
		evals, evecs = self.principal_axes()

		if round(linalg.det(evecs)) == -1:
			evecs[0,2] *= -1
			evecs[1,2] *= -1
			evecs[2,2] *= -1
		if round(linalg.det(evecs)) != 1:
			sys.exit("Error: could not make a rotation matrix while adopting the standard orientation")
		
		rot_matrix = evecs.T
		for atom in self.atom:
			position = np.array([ atom.rx, atom.ry, atom.rz ])
			new_position = np.matmul(rot_matrix, position.T).T
			atom.rx = new_position[0]
			atom.ry = new_position[1]
			atom.rz = new_position[2]

	def translate(self, vector):
	
		new_molecule = deepcopy(self)
		
		for atom in new_molecule.atom:

			atom.rx += vector[0]
			atom.ry += vector[1]
			atom.rz += vector[2]
		
		return new_molecule

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
