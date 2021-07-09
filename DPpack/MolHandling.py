import sys, math
import textwrap
from copy import deepcopy

import numpy as np
from numpy import linalg

from DPpack.PTable import *
from DPpack.SetGlobals import *


#######################################  functions  ######################################

# def center_of_mass(molecule):
	
# 	com = np.zeros(3)
# 	total_mass = 0.0
# 	for atom in molecule:
# 		total_mass += atom['mass']
# 		position = np.array([atom['rx'], atom['ry'], atom['rz']])
# 		com += atom['mass'] * position
	
# 	com = com / total_mass
	
# 	return com



# def center_of_mass_distance(molecule1, molecule2):
	
# 	com1 = center_of_mass(molecule1)
# 	com2 = center_of_mass(molecule2)
# 	dx = com1[0] - com2[0]
# 	dy = com1[1] - com2[1]
# 	dz = com1[2] - com2[2]
# 	distance = math.sqrt(dx**2 + dy**2 + dz**2)
	
# 	return distance



# def center_of_mass_to_origin(molecule):
	
# 	com = center_of_mass(molecule)
# 	for atom in molecule:
# 		atom['rx'] -= com[0]
# 		atom['ry'] -= com[1]
# 		atom['rz'] -= com[2]
	
# 	return



# def charges_and_dipole(molecule):
	
# 	eA_to_Debye = 1/0.20819434
# 	charge = 0
# 	dipole = np.zeros(3)
# 	for atom in molecule:
# 		position = np.array([ atom['rx'], atom['ry'], atom['rz'] ])
# 		dipole += atom['chg'] * position
# 		charge += atom['chg']
	
# 	dipole *= eA_to_Debye
# 	total_dipole = math.sqrt(dipole[0]**2 + dipole[1]**2 + dipole[2]**2)
	
# 	return [charge, dipole[0], dipole[1], dipole[2], total_dipole]



# def distances_between_atoms(molecule):
	
# 	distances = []
# 	dim = len(molecule)
# 	for atom1 in molecule:
# 		if atom1['na'] != ghost_number:
# 			for atom2 in molecule:
# 				if atom2['na'] != ghost_number:
# 					dx = atom1['rx'] - atom2['rx']
# 					dy = atom1['ry'] - atom2['ry']
# 					dz = atom1['rz'] - atom2['rz']
# 					distances.append(math.sqrt(dx**2 + dy**2 + dz**2))
	
# 	return np.array(distances).reshape(dim, dim)



# def eixos(molecule):
	
# 	eixos = np.zeros(3)
# 	if len(molecule) == 2:
# 		position1 = np.array([ molecule[0]['rx'], molecule[0]['ry'], molecule[0]['rz'] ])
# 		position2 = np.array([ molecule[1]['rx'], molecule[1]['ry'], molecule[1]['rz'] ])
# 		eixos = position2 - position1
# 		eixos /= linalg.norm(eixos)
# 	elif len(molecule) > 2:
# 		position1 = np.array([ molecule[0]['rx'], molecule[0]['ry'], molecule[0]['rz'] ])
# 		position2 = np.array([ molecule[1]['rx'], molecule[1]['ry'], molecule[1]['rz'] ])
# 		position3 = np.array([ molecule[2]['rx'], molecule[2]['ry'], molecule[2]['rz'] ])
# 		v1 = position2 - position1
# 		v2 = position3 - position1
# 		v3 = np.cross(v1, v2)
# 		v2 = np.cross(v1, v3)
# 		v1 /= linalg.norm(v1)
# 		v2 /= linalg.norm(v2)
# 		v3 /= linalg.norm(v3)
# 		eixos = np.array([[v1[0], v1[1], v1[2]],
# 		                  [v2[0], v2[1], v2[2]],
# 		                  [v3[0], v3[1], v3[2]]])
	
# 	return eixos



# def inertia_tensor(molecule):
	
# 	com = center_of_mass(molecule)
# 	Ixx = Ixy = Ixz = Iyy = Iyz = Izz = 0.0
# 	for atom in molecule:
# 		####  Obtain the displacement from the center of mass
# 		dx = atom['rx'] - com[0]
# 		dy = atom['ry'] - com[1]
# 		dz = atom['rz'] - com[2]
# 		####  Update the diagonal components of the tensor
# 		Ixx += atom['mass'] * (dy**2 + dz**2)
# 		Iyy += atom['mass'] * (dz**2 + dx**2)
# 		Izz += atom['mass'] * (dx**2 + dy**2)
# 		####  Update the off-diagonal components of the tensor
# 		Ixy += atom['mass'] * dx * dy * -1
# 		Ixz += atom['mass'] * dx * dz * -1
# 		Iyz += atom['mass'] * dy * dz * -1
		
# 	return np.array([[Ixx, Ixy, Ixz],
# 					 [Ixy, Iyy, Iyz],
# 					 [Ixz, Iyz, Izz]])



# def minimum_distance(molecule1, molecule2):
	
# 	distances = []
# 	for atom1 in molecule1:
# 		if atom1['na'] != ghost_number:
# 			for atom2 in molecule2:
# 				if atom2['na'] != ghost_number:
# 					dx = atom1['rx'] - atom2['rx']
# 					dy = atom1['ry'] - atom2['ry']
# 					dz = atom1['rz'] - atom2['rz']
# 					distances.append(math.sqrt(dx**2 + dy**2 + dz**2))
	
# 	return min(distances)



def nearest_image(refmol, molecule, lx, ly, lz, criterium="com"):
	
	if criterium != "com" and criterium != "min":
		sys.exit("Error in value passed to function nearest_image")
	min_dist = 1e20
	for i in range(-1, 2):
		for j in range(-1, 2):
			for k in range(-1, 2):
				
				tr_vector = [i * lx, j * ly, k * lz]
				new_molecule = translate(molecule, tr_vector)
				if criterium == "com":
					dist = center_of_mass_distance(refmol, new_molecule)
				else:
					dist = minimum_distance(refmol, new_molecule)
					
				if dist < min_dist:
					min_dist = dist
					nearestmol = deepcopy(new_molecule)
	
	return min_dist, nearestmol



def calculate_step(gradient, hessian, fh):
	
	invhessian = linalg.inv(hessian)
	pre_step = -1 * np.matmul(invhessian, gradient.T).T
	maxstep = np.amax(np.absolute(pre_step))
	factor = min(1, player['maxstep']/maxstep)
	step = factor * pre_step
	
	fh.write("\nCalculated step:\n")
	pre_step_list = pre_step.tolist()
	
	fh.write("-----------------------------------------------------------------------\n"
			 "Center     Atomic                          Step (Bohr)\n"
			 "Number     Number              X                Y                Z\n"
			 "-----------------------------------------------------------------------\n")
	for i in range(len(molecules[0])):
		fh.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
			   												i + 1, molecules[0][i]['na'], 
			   			pre_step_list.pop(0), pre_step_list.pop(0), pre_step_list.pop(0)))
	
	fh.write("-----------------------------------------------------------------------\n")
	
	fh.write("Maximum step is {:>11.6}\n".format(maxstep))
	fh.write("Scaling factor = {:>6.4f}\n".format(factor))
	fh.write("\nFinal step (Bohr):\n")
	step_list = step.tolist()
	
	fh.write("-----------------------------------------------------------------------\n"
			 "Center     Atomic                          Step (Bohr)\n"
			 "Number     Number              X                Y                Z\n"
			 "-----------------------------------------------------------------------\n")
	for i in range(len(molecules[0])):
		fh.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
			   												i + 1, molecules[0][i]['na'], 
			   						step_list.pop(0), step_list.pop(0), step_list.pop(0)))
	
	fh.write("-----------------------------------------------------------------------\n")
	
	step_max = np.amax(np.absolute(step))
	step_rms = np.sqrt(np.mean(np.square(step)))
	
	fh.write("  Max Step = {:>14.9f}      RMS Step = {:>14.9f}\n\n".format(
																	  step_max, step_rms))
	
	return step



# def read_position(molecule):
	
# 	position_list = []
# 	for atom in molecule:
# 		position_list.extend([ atom['rx'], atom['ry'], atom['rz'] ])
# 	position = np.array(position_list)
# 	position *= ang2bohr
	
# 	return position



def update_molecule(position, fh):
	
	position_in_ang = (position * bohr2ang).tolist()
	new_molecule = deepcopy(molecules[0])
	for atom in new_molecule:
		atom['rx'] = position_in_ang.pop(0)
		atom['ry'] = position_in_ang.pop(0)
		atom['rz'] = position_in_ang.pop(0)
	
	rmsd, molecules[0] = rmsd_fit(new_molecule, molecules[0])
	
	fh.write("\nProjected new conformation of reference molecule with RMSD fit\n")
	fh.write("RMSD = {:>8.5f} Angstrom\n".format(rmsd))
		
	return



# def update_hessian(step, cur_gradient, old_gradient, hessian):	## According to the BFGS
	
# 	dif_gradient = cur_gradient - old_gradient
	
# 	mat1 = 1/np.dot(dif_gradient, step) * np.matmul(dif_gradient.T, dif_gradient)
# 	mat2 = 1/np.dot(step, np.matmul(hessian, step.T).T)
# 	mat2 *= np.matmul( np.matmul(hessian, step.T), np.matmul(step, hessian) )
	
# 	hessian += mat1 - mat2
	
# 	return hessian



def populate_asec_vdw(cycle, fh):
		
	asec_charges = []  	# (rx, ry, rz, chg)
	vdw_meanfield = [] 	# (rx, ry, rz, eps, sig)
	
	if dice['nstep'][-1] % dice['isave'] == 0:
		nconfigs = round(dice['nstep'][-1] / dice['isave'])
	else:
		nconfigs = int(dice['nstep'][-1] / dice['isave'])
	
	norm_factor = nconfigs * player['nprocs']
	
	nsitesref = len(molecules[0]) + len(ghost_atoms) + len(lp_atoms)
	
	nsites_total = dice['nmol'][0] * nsitesref
	for i in range(1, len(dice['nmol'])):
		nsites_total += dice['nmol'][i] * len(molecules[i])
	
	thickness = []
	picked_mols = []
		
	for proc in range(1, player['nprocs'] + 1):  ## Run over folders
		
		path = "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
		file = path + os.sep + dice['outname'] + ".xyz" 
		if not os.path.isfile(file):
			sys.exit("Error: cannot find file {}".format(file))
		try:
			with open(file) as xyzfh:
				xyzfile = xyzfh.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(file))
		
		for config in range(nconfigs):  ## Run over configs in a folder
			
			if int( xyzfile.pop(0).split()[0] ) != nsites_total:
				sys.exit("Error: wrong number of sites in file {}".format(file))
			
			box = xyzfile.pop(0).split()[-3:]
			box = [ float(box[0]), float(box[1]), float(box[2]) ]
			sizes = sizes_of_molecule(molecules[0])
			thickness.append( min([ (box[0] - sizes[0])/2, (box[1] - sizes[1])/2, 
														   (box[2] - sizes[2])/2 ]) )
			
			xyzfile = xyzfile[nsitesref:]  ## Skip the first (reference) molecule
			mol_count = 0
			for type in range(len(dice['nmol'])):  ## Run over types of molecules
				
				if type == 0:
					nmols = dice['nmol'][0] - 1
				else:
					nmols = dice['nmol'][type]
				
				for mol in range(nmols):  ## Run over molecules of each type
				
					new_molecule = []
					for site in range(len(molecules[type])):  ## Run over sites of each molecule
					
						new_molecule.append({})
						line = xyzfile.pop(0).split()
						
						if line[0].title() != atomsymb[molecules[type][site]['na']].strip():
							sys.exit("Error reading file {}".format(file))
						
						new_molecule[site]['na'] = molecules[type][site]['na']
						new_molecule[site]['rx'] = float(line[1])
						new_molecule[site]['ry'] = float(line[2])
						new_molecule[site]['rz'] = float(line[3])
						new_molecule[site]['chg'] = molecules[type][site]['chg']
						new_molecule[site]['eps'] = molecules[type][site]['eps']
						new_molecule[site]['sig'] = molecules[type][site]['sig']
						
					dist = minimum_distance(molecules[0], new_molecule)
					if dist < thickness[-1]:
						mol_count += 1
						for atom in new_molecule:
							asec_charges.append({})
							vdw_meanfield.append({})
							
							asec_charges[-1]['rx'] = atom['rx']
							asec_charges[-1]['ry'] = atom['ry']
							asec_charges[-1]['rz'] = atom['rz']
							asec_charges[-1]['chg'] = atom['chg'] / norm_factor
							
							if player['vdwforces'] == "yes":
								vdw_meanfield[-1]['rx'] = atom['rx']
								vdw_meanfield[-1]['ry'] = atom['ry']
								vdw_meanfield[-1]['rz'] = atom['rz']
								vdw_meanfield[-1]['eps'] = atom['eps']
								vdw_meanfield[-1]['sig'] = atom['sig']
					
					####  Read lines with ghosts or lps in molecules of type 0 (reference)
					####  and, if dist < thickness, appends to asec
					if type == 0:
						for ghost in ghost_atoms:
							line = xyzfile.pop(0).split()
							if line[0] != dice_ghost_label:
								sys.exit("Error reading file {}".format(file))
							if dist < thickness[-1]:
								asec_charges.append({})
								asec_charges[-1]['rx'] = float(line[1])
								asec_charges[-1]['ry'] = float(line[2])
								asec_charges[-1]['rz'] = float(line[3])
								asec_charges[-1]['chg'] = ghost['chg'] / norm_factor
						
						for lp in lp_atoms:
							line = xyzfile.pop(0).split()
							if line[0] != dice_ghost_label:
								sys.exit("Error reading file {}".format(file))
							if dist < thickness[-1]:
								asec_charges.append({})
								asec_charges[-1]['rx'] = float(line[1])
								asec_charges[-1]['ry'] = float(line[2])
								asec_charges[-1]['rz'] = float(line[3])
								asec_charges[-1]['chg'] = lp['chg'] / norm_factor
			
			picked_mols.append(mol_count)
	
	fh.write("Done\n")
	
	string =  "In average, {:^7.2f} molecules ".format(sum(picked_mols)/norm_factor)
	string += "were selected from each of the {} configurations ".format(len(picked_mols))
	string += "of the production simulations to form the ASEC, comprising a shell with "
	string += "minimum thickness of {:>6.2f} Angstrom\n".format(sum(thickness)/norm_factor)
	
	fh.write(textwrap.fill(string, 86))
	fh.write("\n")
	
	otherfh = open("ASEC.dat", "w")
	for charge in asec_charges:
		otherfh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
								 charge['rx'], charge['ry'], charge['rz'], charge['chg']))
	otherfh.close()
	
	return asec_charges



# def principal_axes(inertia_tensor):
	
# 	try:
# 		evals, evecs = linalg.eigh(inertia_tensor)
# 	except:
# 		sys.exit("Error: diagonalization of inertia tensor did not converge")
	
# 	return evals, evecs



def print_geom(cycle, fh):
	
	fh.write("{}\n".format(len(molecules[0])))
	fh.write("Cycle # {}\n".format(cycle))
	for atom in molecules[0]:
		symbol = atomsymb[atom['na']]
		fh.write("{:<2s}    {:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(symbol, 
		                                              atom['rx'], atom['ry'], atom['rz']))
	
	return



def print_mol_info(molecule, fh):
	
	com = center_of_mass(molecule)
	fh.write("    Center of mass = ( {:>10.4f} , {:>10.4f} , {:>10.4f} )\n".format(com[0], 
	                                                                      com[1], com[2]))
	inertia = inertia_tensor(molecule)
	evals, evecs = principal_axes(inertia)
	
	fh.write("    Moments of inertia =  {:>9E}  {:>9E}  {:>9E}\n".format(evals[0], 
	                                                                  evals[1], evals[2]))
	
	fh.write("    Major principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
	                                                  evecs[0,0], evecs[1,0], evecs[2,0]))
	fh.write("    Inter principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
	                                                  evecs[0,1], evecs[1,1], evecs[2,1]))
	fh.write("    Minor principal axis = ( {:>10.6f} , {:>10.6f} , {:>10.6f} )\n".format(
	                                                  evecs[0,2], evecs[1,2], evecs[2,2]))
	
	sizes = sizes_of_molecule(molecule)
	fh.write("    Characteristic lengths = ( {:>6.2f} , {:>6.2f} , {:>6.2f} )\n".format(
	                                                        sizes[0], sizes[1], sizes[2]))
	mol_mass = total_mass(molecule)
	fh.write("    Total mass = {:>8.2f} au\n".format(mol_mass))
	
	chg_dip = charges_and_dipole(molecule)
	fh.write("    Total charge = {:>8.4f} e\n".format(chg_dip[0]))
	fh.write("    Dipole moment = ( {:>9.4f} , {:>9.4f} , {:>9.4f} )     Total = {:>9.4f} Debye\n\n".format(
	                                      chg_dip[1], chg_dip[2], chg_dip[3], chg_dip[4]))
	
	return



# def sizes_of_molecule(molecule):
	
# 	x_list = []
# 	y_list = []
# 	z_list = []
# 	for atom in molecule:
# 		if atom['na'] != ghost_number:
# 			x_list.append(atom['rx'])
# 			y_list.append(atom['ry'])
# 			z_list.append(atom['rz'])
	
# 	x_max = max(x_list)
# 	x_min = min(x_list)
# 	y_max = max(y_list)
# 	y_min = min(y_list)
# 	z_max = max(z_list)
# 	z_min = min(z_list)
	
# 	sizes = [x_max - x_min, y_max - y_min, z_max - z_min]
	
# 	return sizes



# def standard_orientation(molecule):
	
# 	center_of_mass_to_origin(molecule)
# 	tensor = inertia_tensor(molecule)
# 	evals, evecs = principal_axes(tensor)
# 	if round(linalg.det(evecs)) == -1:
# 		evecs[0,2] *= -1
# 		evecs[1,2] *= -1
# 		evecs[2,2] *= -1
# 	if round(linalg.det(evecs)) != 1:
# 		sys.exit("Error: could not make a rotation matrix while adopting the standard orientation")
	
# 	rot_matrix = evecs.T
# 	for atom in molecule:
# 		position = np.array([ atom['rx'], atom['ry'], atom['rz'] ])
# 		new_position = np.matmul(rot_matrix, position.T).T
# 		atom['rx'] = new_position[0]
# 		atom['ry'] = new_position[1]
# 		atom['rz'] = new_position[2]
	
# 	return



# def total_mass(molecule):
	
# 	mass = 0
# 	for atom in molecule:
# 		mass += atom['mass']
	
# 	return mass



# def translate(molecule, vector):
	
# 	new_molecule = deepcopy(molecule)
# 	for atom in new_molecule:
# 		atom['rx'] += vector[0]
# 		atom['ry'] += vector[1]
# 		atom['rz'] += vector[2]
	
# 	return new_molecule



def rmsd_fit(projecting_mol, reference_mol):
	
	if len(projecting_mol) != len(reference_mol):
		sys.exit("Error in RMSD fit procedure: molecules have different number of atoms")
	dim = len(projecting_mol)
	
	new_projecting_mol = deepcopy(projecting_mol)
	new_reference_mol = deepcopy(reference_mol)
	
	center_of_mass_to_origin(new_projecting_mol)
	center_of_mass_to_origin(new_reference_mol)
	
	x = []
	y = []
	
	for atom in new_projecting_mol:
		x.extend([ atom['rx'], atom['ry'], atom['rz'] ])
	
	for atom in new_reference_mol:
		y.extend([ atom['rx'], atom['ry'], atom['rz'] ])
	
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
		new_projecting_mol[i]['rx'] = x[i,0]
		new_projecting_mol[i]['ry'] = x[i,1]
		new_projecting_mol[i]['rz'] = x[i,2]
	
	tr_vector = center_of_mass(reference_mol)
	projected_mol = translate(new_projecting_mol, tr_vector)
	
	return rmsd, projected_mol



