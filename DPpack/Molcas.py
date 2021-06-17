import os, sys
import textwrap
import subprocess

from DPpack.PTable import *
from DPpack.SetGlobals import *
from DPpack.MolHandling import *
from DPpack.Misc import *

#######################################  functions  ######################################

def read_forces_log(file, fh):
	
	forces = []
	try:
		with open(file) as tmpfh:
			logfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = logfile.pop(0).strip()
	while  start.find("Molecular gradients") < 0:	##  expression not found
		start = logfile.pop(0).strip()
	
	logfile = logfile[7:]		##  skip next 7 lines
	
	for i in range(len(molecules[0])):
		values = logfile.pop(0).split()
		values = values[1:]
		forces.extend([ float(x) for x in values ])
	
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



def read_hessian_log(file):
	
	force_const = []
	try:
		with open(file) as tmpfh:
			logfile = tmpfh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	start = logfile.pop(0).strip()
	while  start.find("Force constant matrix") < 0:
		start = logfile.pop(0).strip()
	
	logfile = logfile[1:]		## skip next 1 line
	
	degrees = 3 * len(molecules[0])
	dim = degrees
	last = round(dim * dim)
	count = 0
	while True:
		values = logfile.pop(0).rstrip()
		while len(values) != 0:
			new_value = values[:16]
			values = values[16:]
			force_const.append(float(new_value))
			count += 1
		
		if count >= last:
			break
	
	hessian = np.array(force_const).reshape(dim, dim)
	hessian = hessian[:degrees, :degrees]		## remove degrees related to ghost atoms
	
	for i in range(degrees):
		for j in range(i + 1):
			hessian[j,i] = hessian[i,j]		## force the hessian to be symmetric
	
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



def make_asec_file(cycle, asec_charges):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	file = path + os.sep + "asec.xfield"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("{} Angstrom\n".format(len(asec_charges)))
	
	## Write the ASEC:
	for charge in asec_charges:
		fh.write("{:>10.5f}   {:>10.5f}   {:>10.5f}    {:>11.8f}    0.0  0.0  0.0\n".format(
								 charge['rx'], charge['ry'], charge['rz'], charge['chg']))
	
	fh.write("End of input\n")
	fh.close()
	
	return



def make_force_input(cycle, asec_charges):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	file = path + os.sep + "asec.input"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write(" &Gateway\n")
	fh.write("  Coord\n")
	
	nsites = len(molecules[0])
	if cycle >= player['switchcyc']:
		nsites += len(ghost_atoms) + len(lp_atoms)
	
	fh.write("{}\n".format(nsites))
	fh.write("\nForce calculation - Cycle number {}\n".format(cycle))
	
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
	
	fh.write("basis = {}\n".format(molcas['basis']))
	fh.write("group= nosym\n")
	fh.write(" XFIELD\n")
	fh.write(">>> Include asec.xfield\n")
	
	if not os.path.isfile(molcas['mbottom']):
		sys.exit("Error: cannot find file {} in main directory".format(molcas['mbottom']))
	try:
		with open(molcas['mbottom']) as mbottomfile:
			mbottom = mbottomfile.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(molcas['mbottom']))
	
	for line in mbottom:
		fh.write(line)
	
	fh.write(" &Alaska\nPNEW\n &SLAPAF\nCartesian\n")
	fh.close()
	
	return



def make_charge_input(cycle, asec_charges):
	
	path = "step{:02d}".format(cycle) + os.sep + "qm"
	file = path + os.sep + "asec.input"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write(" &Gateway\n")
	fh.write("  Coord\n")
	
	nsites = len(molecules[0])
	if cycle >= player['switchcyc']:
		nsites += len(ghost_atoms) + len(lp_atoms)
	
	fh.write("{}\n".format(nsites))
	fh.write("\nForce calculation - Cycle number {}\n".format(cycle))
	
	for atom in molecules[0]:
		symbol = atomsymb[atom['na']]
		fh.write("{:<2s}    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(symbol, 
													  atom['rx'], atom['ry'], atom['rz']))
	
	for ghost in ghost_atoms:
		fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
												   ghost['rx'], ghost['ry'], ghost['rz']))
	
	for lp in lp_atoms:
		fh.write("Bq    {:>10.5f}   {:>10.5f}   {:>10.5f}\n".format(
															lp['rx'], lp['ry'], lp['rz']))
	
	fh.write("basis = {}\n".format(molcas['basis']))
	fh.write("group= nosym\n")
	fh.write(" XFIELD\n")
	fh.write(">>> Include asec.xfield\n")
	
	if not os.path.isfile(molcas['mbottom']):
		sys.exit("Error: cannot find file {} in main directory".format(molcas['mbottom']))
	try:
		with open(molcas['mbottom']) as mbottomfile:
			mbottom = mbottomfile.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(molcas['mbottom']))
	
	for line in mbottom:
		fh.write(line)
	
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



