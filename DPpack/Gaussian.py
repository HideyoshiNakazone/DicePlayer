import os, sys
import textwrap
import subprocess

from DPpack.PTable import *
from DPpack.SetGlobals import *
from DPpack.MolHandling import *
from DPpack.Misc import *

#######################################  functions  ######################################

def read_forces_fchk(file, fh):
	
	forces = []
	try:
		with open(file) as tmpfh:
			fchkfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = fchkfile.pop(0).strip()
	while  start.find("Cartesian Gradient") != 0:	##  expression in begining of line
		start = fchkfile.pop(0).strip()
	
	degrees = 3 * len(molecules[0])
	count = 0
	while True:
		values = fchkfile.pop(0).split()
		forces.extend([ float(x) for x in values ])
		count += len(values)
		if count >= degrees:
			forces = forces[:degrees]
			break
	
	gradient = np.array(forces)
	
	fh.write("\nGradient read from file {}:\n".format(file))
	fh.write("-----------------------------------------------------------------------\n"
			 "Center     Atomic                     Forces (Hartree/Bohr)\n"
			 "Number     Number              X                Y                Z\n"
			 "-----------------------------------------------------------------------\n")
	for i in range(len(molecules[0])):
		fh.write("  {:>5d}     {:>3d}        {:>14.9f}   {:>14.9f}   {:>14.9f}\n".format(
			   i + 1, molecules[0][i]['na'], forces.pop(0), forces.pop(0), forces.pop(0)))
	
	fh.write("-----------------------------------------------------------------------\n")
	
	force_max = np.amax(np.absolute(gradient))
	force_rms = np.sqrt(np.mean(np.square(gradient)))
	
	fh.write("  Max Force = {:>14.9f}      RMS Force = {:>14.9f}\n\n".format(
																	force_max, force_rms))
	
	return gradient



def read_hessian_fchk(file):
	
	force_const = []
	try:
		with open(file) as tmpfh:
			fchkfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = fchkfile.pop(0).strip()
	while  start.find("Cartesian Force Constants") != 0:
		start = fchkfile.pop(0).strip()
	
	degrees = 3 * len(molecules[0])
	last = round(degrees * (degrees + 1) / 2)
	count = 0
	while True:
		values = fchkfile.pop(0).split()
		force_const.extend([ float(x) for x in values ])
		count += len(values)
		if count >= last:
			force_const = force_const[:last]
			break
	
	hessian = np.zeros((degrees, degrees))
	for i in range(degrees):
		for j in range(i + 1):
			hessian[i,j] = force_const.pop(0)
			hessian[j,i] = hessian[i,j]
	
	return hessian



def read_hessian_log(file):
	
	try:
		with open(file) as tmpfh:
			logfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = logfile.pop(0).strip()
	while  start.find("The second derivative matrix:") != 0:
		start = logfile.pop(0).strip()
	
	degrees = 3 * len(molecules[0])
	hessian = np.zeros((degrees, degrees))
	
	k = 0
	while k < degrees:
		logfile.pop(0)
		for i in range(k, degrees):
			values = logfile.pop(0).split()[1:]
			for j in range(k, min(i + 1, k + 5)):
				hessian[i,j] = float(values.pop(0))
				hessian[j,i] = hessian[i,j]
		k += 5
	
	return hessian



def print_grad_hessian(cycle, cur_gradient, hessian):
	
	try:
		fh = open("grad_hessian.dat", "w")
	except:
		sys.exit("Error: cannot open file grad_hessian.dat")
	
	fh.write("Optimization cycle: {}\n".format(cycle))
	fh.write("Cartesian Gradient\n")
	degrees = 3 * len(molecules[0])
	for i in range(degrees):
		fh.write(" {:>11.8g}".format(cur_gradient[i]))
		if (i + 1) % 5 == 0 or i == degrees - 1:
			fh.write("\n")
	
	fh.write("Cartesian Force Constants\n")
	last = degrees * (degrees + 1) / 2
	count = 0
	for i in range(degrees):
		for j in range(i + 1):
			count += 1
			fh.write(" {:>11.8g}".format(hessian[i,j]))
			if count % 5 == 0 or count == last:
				fh.write("\n")
	
	fh.close()
	
	return



def make_force_input(cycle, asec_charges):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	file = path + os.sep + "asec.gjf"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("%Chk=asec.chk\n")
	if gaussian['mem'] != None:
		fh.write("%Mem={}MB\n".format(gaussian['mem']))
	fh.write("%Nprocs={}\n".format(player['nprocs'] * dice['ncores']))
	
	kword_line = "#P " + gaussian['level'] + " " + gaussian['keywords']
	kword_line += " Force Charge NoSymm"
	
	if cycle >= player['switchcyc']:
		kword_line += " Pop={} Density=Current".format(gaussian['pop'])
	
	if cycle > 1:
		kword_line += " Guess=Read"
	
	fh.write(textwrap.fill(kword_line, 90))
	fh.write("\n")
	
	fh.write("\nForce calculation - Cycle number {}\n".format(cycle))
	fh.write("\n")
	fh.write("{},{}\n".format(gaussian['chgmult'][0], gaussian['chgmult'][1]))
	
	for atom in molecules[0]:
		symbol = atomsymb[atom['na']]
		fh.write("{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(symbol, 
													  atom['rx'], atom['ry'], atom['rz']))
	
	## If also performing charge fit in the same calculation
	if cycle >= player['switchcyc']:
		for ghost in ghost_atoms:
			fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
												   ghost['rx'], ghost['ry'], ghost['rz']))
		
		for lp in lp_atoms:
			fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
															lp['rx'], lp['ry'], lp['rz']))
	
	fh.write("\n")
	
	## If gmiddle file was informed, write its contents in asec.gjf
	if gaussian['gmiddle'] != None:
		if not os.path.isfile(gaussian['gmiddle']):
			sys.exit("Error: cannot find file {} in main directory".format(
																	 gaussian['gmiddle']))
		try:
			with open(gaussian['gmiddle']) as gmiddlefile:
				gmiddle = gmiddlefile.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(gaussian['gmiddle']))
		
		for line in gmiddle:
			fh.write(line)
		
		fh.write("\n")
	
	## Write the ASEC:
	for charge in asec_charges:
		fh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
								 charge['rx'], charge['ry'], charge['rz'], charge['chg']))
	
	fh.write("\n")
	
	## If gbottom file was informed, write its contents in asec.gjf
	if gaussian['gbottom'] != None:
		if not os.path.isfile(gaussian['gbottom']):
			sys.exit("Error: cannot find file {} in main directory".format(
																	 gaussian['gbottom']))
		try:
			with open(gaussian['gbottom']) as gbottomfile:
				gbottom = gbottomfile.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(gaussian['gbottom']))
		
		for line in gbottom:
			fh.write(line)
		
		fh.write("\n")
	
	fh.close()
	
	return



def make_charge_input(cycle, asec_charges):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	file = path + os.sep + "asec2.gjf"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("%Chk=asec.chk\n")
	if gaussian['mem'] != None:
		fh.write("%Mem={}MB\n".format(gaussian['mem']))
	fh.write("%Nprocs={}\n".format(player['nprocs'] * dice['ncores']))
	
	kword_line = "#P " + gaussian['chglevel'] + " " + gaussian['keywords'] + " Charge NoSymm"
	
	if player['opt'] != "no" or cycle > 1:
		kword_line += " Guess=Read"
	
	kword_line += " Pop={} Density=Current\n".format(gaussian['pop'])
	
	fh.write(textwrap.fill(kword_line, 90))
	fh.write("\n")
	
	fh.write("\nCharge calculation - Cycle number {}\n".format(cycle))
	fh.write("\n")
	fh.write("{},{}\n".format(gaussian['chgmult'][0], gaussian['chgmult'][1]))
	
	for atom in molecules[0]:
		symbol = atomsymb[atom['na']]
		fh.write("{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(symbol, 
													  atom['rx'], atom['ry'], atom['rz']))
	
	if cycle >= player['switchcyc']:
		for ghost in ghost_atoms:
			fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
												   ghost['rx'], ghost['ry'], ghost['rz']))
		
		for lp in lp_atoms:
			fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
															lp['rx'], lp['ry'], lp['rz']))
	
	fh.write("\n")
	
	## If gmiddle file was informed, write its contents in asec.gjf
	if gaussian['gmiddle'] != None:
		if not os.path.isfile(gaussian['gmiddle']):
			sys.exit("Error: cannot find file {} in main directory".format(
																	 gaussian['gmiddle']))
		try:
			with open(gaussian['gmiddle']) as gmiddlefile:
				gmiddle = gmiddlefile.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(gaussian['gmiddle']))
		
		for line in gmiddle:
			fh.write(line)
		
		fh.write("\n")
	
	## Write the ASEC:
	for charge in asec_charges:
		fh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}     {:>11.8f}\n".format(
								 charge['rx'], charge['ry'], charge['rz'], charge['chg']))
	
	fh.write("\n")
	
	## If gbottom file was informed, write its contents in asec.gjf
	if gaussian['gbottom'] != None:
		if not os.path.isfile(gaussian['gbottom']):
			sys.exit("Error: cannot find file {} in main directory".format(
																	 gaussian['gbottom']))
		try:
			with open(gaussian['gbottom']) as gbottomfile:
				gbottom = gbottomfile.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(gaussian['gbottom']))
		
		for line in gbottom:
			fh.write(line)
		
		fh.write("\n")
	
	fh.close()
	
	return



def read_charges(file, fh):
	
	try:
		with open(file) as tmpfh:
			glogfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = glogfile.pop(0).strip()
	while  start != "Fitting point charges to electrostatic potential":
		start = glogfile.pop(0).strip()
	
	glogfile = glogfile[3:]		## Consume 3 more lines
	
	fh.write("\nAtomic charges:\n")
	fh.write("------------------------------------\n")
	for atom in molecules[0]:
		line = glogfile.pop(0).split()
		atom_str = line[1]
		charge = float(line[2])
		atom['chg'] = charge
		fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))
	
	if gaussian['pop'] == "chelpg":
		for ghost in ghost_atoms:
			line = glogfile.pop(0).split()
			atom_str = line[1]
			charge = float(line[2])
			ghost['chg'] = charge
			fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))
		
		for lp in lp_atoms:
			line = glogfile.pop(0).split()
			atom_str = line[1]
			charge = float(line[2])
			lp['chg'] = charge
			fh.write(" {:<2s}      {:>10.6f}\n".format(atom_str, charge))
	
	fh.write("------------------------------------\n")
		
	return



def run_gaussian(cycle, type, fh):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	work_dir = os.getcwd()
	os.chdir(path)
	
	if type == "force":
		infile = "asec.gjf"
	elif type == "charge":
		infile = "asec2.gjf"
	
	fh.write("\nCalculation of {}s initiated with Gaussian on {}\n".format(type, date_time()))
	
	exit_status = subprocess.call([player['qmprog'], infile])
	
	if exit_status != 0:
		sys.exit("Gaussian process did not exit properly")
	
	fh.write("Calculation of {}s finished on {}\n".format(type, date_time()))
	
	os.chdir(work_dir)
	
	return



def run_formchk(cycle, fh):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	work_dir = os.getcwd()
	os.chdir(path)
		
	fh.write("Formatting the checkpoint file... ")
	
	exit_status = subprocess.call(["formchk", "asec.chk"])
	
	fh.write("Done\n")
	
	os.chdir(work_dir)
	
	return



