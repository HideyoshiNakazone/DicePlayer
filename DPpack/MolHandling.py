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

env = ["OMP_STACKSIZE"]

bohr2ang = 0.52917721092
ang2bohr = 1/bohr2ang

# Usaremos uma nova classe que ira conter toda interação entre moleculas

class System:

	def __init__(self):

		self.molecule = []
		self.nmols = []

	def add_type(self, nmols, m):
		
		self.molecule.append(m)
		self.nmols.append(nmols)

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

	def rmsd_fit(self, index_p, index_r):

		projecting_mol = self.molecule[index_p]
		reference_mol = self.molecule[index_r]
		
		if len(projecting_mol.atom) != len(reference_mol.atom):
			sys.exit("Error in RMSD fit procedure: molecules have different number of atoms")
		dim = len(projecting_mol.atom)
		
		new_projecting_mol = deepcopy(projecting_mol)
		new_reference_mol = deepcopy(reference_mol)
		
		new_projecting_mol.center_of_mass_to_origin()
		new_reference_mol.center_of_mass_to_origin()
		
		x = []
		y = []
		
		for atom in new_projecting_mol:
			x.extend([ atom.rx, atom.ry, atom.rz ])
		
		for atom in new_reference_mol:
			y.extend([ atom.rx, atom.ry, atom.rz ])
		
		x = np.array(x).reshape(dim, 3)
		y = np.array(y).reshape(dim, 3)
		
		r = np.matmul(y.T, x)
		rr = np.matmul(r.T, r)
		
		try:
			evals, evecs = linalg.eigh(rr)
		except:
			sys.exit("Error: diagonalization of RR matrix did not converge")
		
		a1 = evecs[:,2].T
		a2 = evecs[:,1].T
		a3 = np.cross(a1, a2)
		
		A = np.array([ a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], a3[0], a3[1], a3[2] ])
		A = A.reshape(3,3)
		
		b1 = np.matmul(r, a1.T).T		# or np.dot(r, a1)
		b1 /= linalg.norm(b1)
		b2 = np.matmul(r, a2.T).T		# or np.dot(r, a2)
		b2 /= linalg.norm(b2)
		b3 = np.cross(b1, b2)
		
		B = np.array([ b1[0], b1[1], b1[2], b2[0], b2[1], b2[2], b3[0], b3[1], b3[2] ])
		B = B.reshape(3,3).T
		
		rot_matrix = np.matmul(B, A)
		x = np.matmul(rot_matrix, x.T).T
		
		rmsd = 0
		for i in range(dim):
			rmsd += (x[i,0] - y[i,0])**2 + (x[i,1] - y[i,1])**2 + (x[i,2] - y[i,2])**2
		rmsd = math.sqrt(rmsd/dim)
		
		for i in range(dim):
			new_projecting_mol.atom[i].rx = x[i,0]
			new_projecting_mol.atom[i].ry = x[i,1]
			new_projecting_mol.atom[i].rz = x[i,2]
		
		tr_vector = reference_mol.center_of_mass()
		projected_mol = new_projecting_mol.translate(tr_vector)
		
		return rmsd, projected_mol

	def update_molecule(self, position, fh):
	
		position_in_ang = (position * bohr2ang).tolist()
		self.add_molecule(deepcopy(self.molecule[0]))

		for atom in self.molecule[-1].atom:

			atom.rx = position_in_ang.pop(0)
			atom.ry = position_in_ang.pop(0)
			atom.rz = position_in_ang.pop(0)
		
		rmsd, self.molecule[0] = self.rmsd_fit(-1, 0)
		self.molecule.pop(-1)
		
		fh.write("\nProjected new conformation of reference molecule with RMSD fit\n")
		fh.write("RMSD = {:>8.5f} Angstrom\n".format(rmsd))

	def nearest_image(self, index_r, index_m, lx, ly, lz, criterium=None):
		
		if criterium in None:
			criterium = "com"

		if criterium != "com" and criterium != "min":
			sys.exit("Error in value passed to function nearest_image")
		
		min_dist = 1e20

		for i in range(-1, 2):
			for j in range(-1, 2):
				for k in range(-1, 2):
					
					tr_vector = [i * lx, j * ly, k * lz]
					self.add_molecule(self.molecule[index_m].translate(tr_vector))

					if criterium == "com":
						dist = self.center_of_mass_distance(index_r, -1)
					else:
						dist = self.minimum_distance(index_r, -1)
						
					if dist < min_dist:
						min_dist = dist
						nearestmol = deepcopy(self.molecule[-1])
		
		self.molecule.pop(-1)

		return min_dist, nearestmol

	def print_geom(self, cycle, fh):
	
		fh.write("{}\n".format(len(self.molecule[0].atom)))
		fh.write("Cycle # {}\n".format(cycle))
		for atom in self.molecule[0].atom:
			symbol = atomsymb[atom.na]
			fh.write("{:<2s}    {:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(symbol, 
														atom.rx, atom.ry, atom.rz))


	
# Classe que conterá toda informação e funções relacionadas a uma unica molecula

class Molecule:

	def __init__(self, molname):

		self.molname = molname

		self.atom = []					# Lista de instancias de Atom
		self.position = None			# Array Numpy
		self.energy = None				# Array Numpy
		self.gradient = None			# Array Numpy
		self.hessian = None				# Array Numpy

		self.total_mass = 0
		self.com = None
		self.ghost_atoms = [] # Stores the index of the ghost atoms in the atoms array
		self.lp_atoms = []

	def add_atom(self, a):

		self.atom.append(a)			# Inserção de um novo atomo
		self.total_mass += a.mass

		if (a.na == ghost_number):

			self.ghost_atoms.append(self.atom.index(a))

		self.center_of_mass()

	def center_of_mass(self):
	
		self.com = np.zeros(3)

		for atom in self.atom:
			
			self.com += atom.mass * np.array([atom.rx, atom.ry, atom.rz])
		
		self.com = self.com / self.total_mass

	def center_of_mass_to_origin(self):
	
		self.center_of_mass()

		for atom in self.atom:

			atom.rx -= self.com[0]
			atom.ry -= self.com[1]
			atom.rz -= self.com[2]
	
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

	def inertia_tensor(self):
		
		self.center_of_mass()
		Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0.0

		for atom in self.atom:

			####  Obtain the displacement from the center of mass
			dx = atom.rx - self.com[0]
			dy = atom.ry - self.com[1]
			dz = atom.rz - self.com[2]
			####  Update the diagonal components of the tensor
			Ixx += atom.mass * (dy**2 + dz**2)
			Iyy += atom.mass * (dz**2 + dx**2)
			Izz += atom.mass * (dx**2 + dy**2)
			####  Update the off-diagonal components of the tensor
			Ixy += atom.mass * dx * dy * -1
			Ixz += atom.mass * dx * dz * -1
			Iyz += atom.mass * dy * dz * -1
			
		return np.array([[Ixx, Ixy, Ixz],
					 	[Ixy, Iyy, Iyz],
					 	[Ixz, Iyz, Izz]])

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
		mat2 *= np.matmul( np.matmul(self.hessian, step.T), np.matmul(step, self.hessian) )
		
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

	def print_mol_info(self, fh):

		fh.write("    Center of mass = ( {:>10.4f} , {:>10.4f} , {:>10.4f} )\n".format(self.com[0], 
																			self.com[1], self.com[2]))
		inertia = self.inertia_tensor()
		evals, evecs = self.principal_axes()
		
		fh.write("    Moments of inertia =  {:>9E}  {:>9E}  {:>9E}\n".format(evals[0], 
																		evals[1], evals[2]))
		
		fh.write("    Major principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
														evecs[0,0], evecs[1,0], evecs[2,0]))
		fh.write("    Inter principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
														evecs[0,1], evecs[1,1], evecs[2,1]))
		fh.write("    Minor principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
														evecs[0,2], evecs[1,2], evecs[2,2]))
		
		sizes = self.sizes_of_molecule()
		fh.write("    Characteristic lengths = ( {:>6.2f} , {:>6.2f} , {:>6.2f} )\n".format(
																sizes[0], sizes[1], sizes[2]))
		fh.write("    Total mass = {:>8.2f} au\n".format(self.total_mass))
		
		chg_dip = self.charges_and_dipole()
		fh.write("    Total charge = {:>8.4f} e\n".format(chg_dip[0]))
		fh.write("    Dipole moment = ( {:>9.4f} , {:>9.4f} , {:>9.4f} )     Total = {:>9.4f} Debye\n\n".format(
											chg_dip[1], chg_dip[2], chg_dip[3], chg_dip[4]))

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