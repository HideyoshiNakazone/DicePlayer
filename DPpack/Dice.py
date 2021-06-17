import sys, os, time
import subprocess

from copy import deepcopy

from numpy import random

from DPpack.PTable import *
from DPpack.SetGlobals import *
from DPpack.Misc import *


dice_end_flag = "End of simulation"		## The normal end flag
dice_flag_line = -2    					## must be in the line before the last
umaAng3_to_gcm3 = 1.6605				## Conversion between uma/Ang3 to g/cm3

max_seed = 4294967295					## Maximum allowed value for a seed (numpy)

#######################################  functions  ######################################

def make_inputs(cycle, proc):
	
	step_dir = "step{:02d}".format(cycle)
	proc_dir = "p{:02d}".format(proc)
	path = step_dir + os.sep + proc_dir
	
	num = time.time()					##  Take the decimal places 7 to 12 of the 
	num = (num - int(num)) * 1e6		##  time in seconds as a floating point
	num = int((num - int(num)) * 1e6)	##  to make an integer in the range 1-1e6
	random.seed( (os.getpid() * num) % (max_seed + 1) )

	if not dice['randominit']:
		xyzfile = dice['outname'] + ".xyz.last-" + "p{:02d}".format(proc)
		make_init_file(path, xyzfile)
	
	if len(dice['nstep']) == 2:  ## Means NVT simulation
		
		make_nvt_ter(path)
		make_nvt_eq(path)
		
	elif len(dice['nstep']) == 3:  ## Means NPT simulation
		
		if dice['randominit']:
			make_nvt_ter(path)
		else:
			dice['dens'] = new_density(proc)
		
		make_npt_ter(path)
		make_npt_eq(path)
	
	else:
		sys.exit("Error: bad number of entries for 'nstep'")
	
	make_potential(path)
			
	return


def make_nvt_ter(path):
	
	file = path + os.sep + "NVT.ter"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("title = {} - NVT Thermalization\n".format(dice['title']))
	fh.write("ncores = {}\n".format(dice['ncores']))
	fh.write("ljname = {}\n".format(dice['ljname']))
	fh.write("outname = {}\n".format(dice['outname']))
	
	string = " ".join(str(x) for x in dice['nmol'])
	fh.write("nmol = {}\n".format(string))
	
	fh.write("dens = {}\n".format(dice['dens']))
	fh.write("temp = {}\n".format(dice['temp']))
	
	if dice['randominit']:
		fh.write("init = yes\n")
		fh.write("nstep = {}\n".format(dice['nstep'][0]))
	else:
		fh.write("init = yesreadxyz\n")
		fh.write("nstep = {}\n".format(player['altsteps']))
	
	fh.write("vstep = 0\n")
	fh.write("mstop = 1\n")
	fh.write("accum = no\n")
	fh.write("iprint = 1\n")
	fh.write("isave = 0\n")
	fh.write("irdf = 0\n")
	
	seed = int(1e6 * random.random())
	fh.write("seed = {}\n".format(seed))
	
	fh.close()
	
	return


def make_nvt_eq(path):
	
	file = path + os.sep + "NVT.eq" 	
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("title = {} - NVT Production\n".format(dice['title']))
	fh.write("ncores = {}\n".format(dice['ncores']))
	fh.write("ljname = {}\n".format(dice['ljname']))
	fh.write("outname = {}\n".format(dice['outname']))
	
	string = " ".join(str(x) for x in dice['nmol'])
	fh.write("nmol = {}\n".format(string))
	
	fh.write("dens = {}\n".format(dice['dens']))
	fh.write("temp = {}\n".format(dice['temp']))
	fh.write("init = no\n")
	fh.write("nstep = {}\n".format(dice['nstep'][1]))
	fh.write("vstep = 0\n")
	fh.write("mstop = 1\n")
	fh.write("accum = no\n")
	fh.write("iprint = 1\n")
	fh.write("isave = {}\n".format(dice['isave']))
	fh.write("irdf = {}\n".format(10 * player['nprocs']))
	
	seed = int(1e6 * random.random())
	fh.write("seed = {}\n".format(seed))
	
	fh.close()
	
	return


def make_npt_ter(path):
	
	file = path + os.sep + "NPT.ter"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("title = {} - NPT Thermalization\n".format(dice['title']))
	fh.write("ncores = {}\n".format(dice['ncores']))
	fh.write("ljname = {}\n".format(dice['ljname']))
	fh.write("outname = {}\n".format(dice['outname']))
	
	string = " ".join(str(x) for x in dice['nmol'])
	fh.write("nmol = {}\n".format(string))
	
	fh.write("press = {}\n".format(dice['press']))
	fh.write("temp = {}\n".format(dice['temp']))
	
	if dice['randominit']:
		fh.write("init = no\n")   ## Because there will be a previous NVT simulation
		fh.write("vstep = {}\n".format(int(dice['nstep'][1] / 5)))
	else:
		fh.write("init = yesreadxyz\n")
		fh.write("dens = {:<8.4f}\n".format(dice['dens']))
		fh.write("vstep = {}\n".format(int(player['altsteps'] / 5)))
	
	fh.write("nstep = 5\n")
	fh.write("mstop = 1\n")
	fh.write("accum = no\n")
	fh.write("iprint = 1\n")
	fh.write("isave = 0\n")
	fh.write("irdf = 0\n")
	
	seed = int(1e6 * random.random())
	fh.write("seed = {}\n".format(seed))
	
	fh.close()
	
	return


def make_npt_eq(path):
	
	file = path + os.sep + "NPT.eq"
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("title = {} - NPT Production\n".format(dice['title']))
	fh.write("ncores = {}\n".format(dice['ncores']))
	fh.write("ljname = {}\n".format(dice['ljname']))
	fh.write("outname = {}\n".format(dice['outname']))
	
	string = " ".join(str(x) for x in dice['nmol'])
	fh.write("nmol = {}\n".format(string))
	
	fh.write("press = {}\n".format(dice['press']))
	fh.write("temp = {}\n".format(dice['temp']))
	
	fh.write("nstep = 5\n")
	
	fh.write("vstep = {}\n".format(int(dice['nstep'][2] / 5)))
	fh.write("init = no\n")
	fh.write("mstop = 1\n")
	fh.write("accum = no\n")
	fh.write("iprint = 1\n")
	fh.write("isave = {}\n".format(dice['isave']))
	fh.write("irdf = {}\n".format(10 * player['nprocs']))
	
	seed = int(1e6 * random.random())
	fh.write("seed = {}\n".format(seed))
	
	fh.close()
	
	return


def make_init_file(path, file):
	
	if not os.path.isfile(file):
		sys.exit("Error: cannot find the xyz file {} in main directory".format(file))
	try:
		with open(file) as fh:
			xyzfile = fh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	nsites_mm = 0
	for i in range(1, len(dice['nmol'])):
		nsites_mm += dice['nmol'][i] * len(molecules[i])
	
	nsites_mm *= -1  ## Become an index to count from the end of xyzfile (list)
	xyzfile = xyzfile[nsites_mm :]  ## Only the MM atoms of the last configuration remains
	
	file = path + os.sep + dice['outname'] + ".xy"
	
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	for atom in molecules[0]:
		fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(atom['rx'], atom['ry'], 
																			  atom['rz']))
	
	for ghost in ghost_atoms:
		fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(ghost['rx'], ghost['ry'], 
																			 ghost['rz']))
	
	for lps in lp_atoms:
		fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(lps['rx'], lps['ry'], 
																			   lps['rz']))
	
	for line in xyzfile:
		atom = line.split()
		rx = float(atom[1])
		ry = float(atom[2])
		rz = float(atom[3])
		fh.write("{:>10.5f}  {:>10.5f}  {:>10.5f}\n".format(rx, ry, rz))
	
	fh.write("$end")
	
	fh.close()
	
	return


def make_potential(path):
	
	fstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f}\n"
	
	file = path + os.sep + dice['ljname']
	try:
		fh = open(file, "w")
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	fh.write("{}\n".format(dice['combrule']))
	fh.write("{}\n".format(len(dice['nmol'])))
	
	nsites_qm = len(molecules[0]) + len(ghost_atoms) + len(lp_atoms)
	
	## Print the sites of the QM molecule
	fh.write("{}\n".format(nsites_qm))
	for atom in molecules[0]:
		fh.write(fstr.format(atom['lbl'], atom['na'], atom['rx'], atom['ry'], atom['rz'], 
			             			   			   atom['chg'], atom['eps'], atom['sig']))
	
	ghost_label = molecules[0][-1]['lbl'] + 1
	for ghost in ghost_atoms:
		fh.write(fstr.format(ghost_label, ghost_number, ghost['rx'], ghost['ry'], 
			             			   					 ghost['rz'], ghost['chg'], 0, 0))
	
	ghost_label += 1
	for lp in lp_atoms:
		fh.write(fstr.format(ghost_label, ghost_number, lp['rx'], lp['ry'], lp['rz'], 
			             			   									 lp['chg'], 0, 0))
	
	## Print the sites of the other molecules
	for mol in molecules[1:]:
		fh.write("{}\n".format(len(mol)))
		for atom in mol:
			fh.write(fstr.format(atom['lbl'], atom['na'], atom['rx'], atom['ry'], 
			             			   atom['rz'], atom['chg'], atom['eps'], atom['sig']))
	
	return


def make_proc_dir(cycle, proc):
	
	step_dir = "step{:02d}".format(cycle)
	proc_dir = "p{:02d}".format(proc)
	path = step_dir + os.sep + proc_dir
	try:
		os.makedirs(path)
	except:
		sys.exit("Error: cannot make directory {}".format(path))
	
	return



def run_dice(cycle, proc, fh):
	
	step_dir = "step{:02d}".format(cycle)
	proc_dir = "p{:02d}".format(proc)
	path = step_dir + os.sep + proc_dir
	working_dir = os.getcwd()
	os.chdir(path)
	
	fh.write("Simulation process {} initiated with pid {}\n".format(proc_dir, os.getpid()))
	
	if len(dice['nstep']) == 2:  ## Means NVT simulation
		
		## NVT thermalization
		string = "(from " + ("random" if dice['randominit'] else "previous") + " configuration)"
		fh.write("p{:02d}> NVT thermalization initiated {} on {}\n".format(proc, string, 
																			 date_time()))
		
		infh = open("NVT.ter")
		outfh = open("NVT.ter.out", "w")
		
		exit_status = subprocess.call(dice['progname'], stdin=infh, stdout=outfh)
		infh.close()
		outfh.close()
		
		if os.getppid() == 1:	## Parent process is dead
			sys.exit()
		
		if exit_status != 0:
			sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		else:
			outfh = open("NVT.ter.out")          ## Open again to seek the normal end flag
			flag = outfh.readlines()[dice_flag_line].strip()
			outfh.close()
			if flag != dice_end_flag:
				sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		
		## NVT production
		fh.write("p{:02d}> NVT production initiated on {}\n".format(proc, date_time()))
		
		infh = open("NVT.eq")
		outfh = open("NVT.eq.out", "w")
		
		exit_status = subprocess.call(dice['progname'], stdin=infh, stdout=outfh)
		infh.close()
		outfh.close()
		
		if os.getppid() == 1:	## Parent process is dead
			sys.exit()
		
		if exit_status != 0:
			sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		else:
			outfh = open("NVT.eq.out")           ## Open again to seek the normal end flag
			flag = outfh.readlines()[dice_flag_line].strip()
			outfh.close()
			if flag != dice_end_flag:
				sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		
		fh.write("p{:02d}> ----- NVT production finished on {}\n".format(proc, 
																			 date_time()))
	
	elif len(dice['nstep']) == 3:  ## Means NPT simulation
		
		## NVT thermalization if randominit
		if dice['randominit']:
			string = "(from random configuration)"
			fh.write("p{:02d}> NVT thermalization initiated {} on {}\n".format(proc, 
																	 string, date_time()))
			infh = open("NVT.ter")
			outfh = open("NVT.ter.out", "w")
			
			exit_status = subprocess.call(dice['progname'], stdin=infh, stdout=outfh)
			infh.close()
			outfh.close()
			
			if os.getppid() == 1:	## Parent process is dead
				sys.exit()
		
			if exit_status != 0:
				sys.exit("Dice process p{:02d} did not exit properly".format(proc))
			else:
				outfh = open("NVT.ter.out")      ## Open again to seek the normal end flag
				flag = outfh.readlines()[dice_flag_line].strip()
				outfh.close()
				if flag != dice_end_flag:
					sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		
		## NPT thermalization
		string = (" (from previous configuration) " if not dice['randominit'] else " ")
		fh.write("p{:02d}> NPT thermalization initiated{}on {}\n".format(proc, string, 
																			 date_time()))
		infh = open("NPT.ter")
		outfh = open("NPT.ter.out", "w")
		
		exit_status = subprocess.call(dice['progname'], stdin=infh, stdout=outfh)
		infh.close()
		outfh.close()
		
		if os.getppid() == 1:	## Parent process is dead
			sys.exit()
		
		if exit_status != 0:
			sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		else:
			outfh = open("NPT.ter.out")          ## Open again to seek the normal end flag
			flag = outfh.readlines()[dice_flag_line].strip()
			outfh.close()
			if flag != dice_end_flag:
				sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		
		## NPT production
		fh.write("p{:02d}> NPT production initiated on {}\n".format(proc, date_time()))
		
		infh = open("NPT.eq")
		outfh = open("NPT.eq.out", "w")
		
		exit_status = subprocess.call(dice['progname'], stdin=infh, stdout=outfh)
		infh.close()
		outfh.close()
		
		if os.getppid() == 1:	## Parent process is dead
			sys.exit()
		
		if exit_status != 0:
			sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		else:
			outfh = open("NPT.eq.out")           ## Open again to seek the normal end flag
			flag = outfh.readlines()[dice_flag_line].strip()
			outfh.close()
			if flag != dice_end_flag:
				sys.exit("Dice process p{:02d} did not exit properly".format(proc))
		
		fh.write("p{:02d}> ----- NPT production finished on {}\n".format(proc, 
																			 date_time()))
		
	os.chdir(working_dir)
	
	return



def print_last_config(cycle, proc):
	
	step_dir = "step{:02d}".format(cycle)
	proc_dir = "p{:02d}".format(proc)
	path = step_dir + os.sep + proc_dir
	file = path + os.sep + dice['outname'] + ".xyz"
	if not os.path.isfile(file):
		sys.exit("Error: cannot find the xyz file {}".format(file))
	try:
		with open(file) as fh:
			xyzfile = fh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	nsites = ( len(molecules[0]) + len(ghost_atoms) + len(lp_atoms) ) * dice['nmol'][0]	
	for i in range(1, len(dice['nmol'])):
		nsites += dice['nmol'][i] * len(molecules[i])
	
	nsites += 2  ## To include the comment line and the number of atoms (xyz file format)
	
	nsites *= -1  ## Become an index to count from the end of xyzfile (list)
	xyzfile = xyzfile[nsites :]  ## Take the last configuration
	
	
	file = dice['outname'] + ".xyz.last-" + proc_dir
	fh = open(file, "w")
	for line in xyzfile:
		fh.write(line)
	
	fh.close()
	
	return



def new_density(proc):
	
	file = dice['outname'] + ".xyz.last-" + "p{:02d}".format(proc)
	if not os.path.isfile(file):
		sys.exit("Error: cannot find the xyz file {} in main directory".format(file))
	try:
		with open(file) as fh:
			xyzfile = fh.readlines()
	except:
		sys.exit("Error: cannot open file {}".format(file))
	
	box = xyzfile[1].split()
	volume = float(box[-3]) * float(box[-2]) * float(box[-1])
	
	total_mass = 0
	for i in range(len(molecules)):
		mol_mass = 0
		for atom in molecules[i]:
			mol_mass += atom['mass']
		total_mass += mol_mass * dice['nmol'][i]
	
	density = (total_mass / volume) * umaAng3_to_gcm3
	
	return density



def simulation_process(cycle, proc, logfh):
	
	try:
		make_proc_dir(cycle, proc)
		make_inputs(cycle, proc)
		run_dice(cycle, proc, logfh)
	except Exception as err:
		sys.exit(err)
	
	return



