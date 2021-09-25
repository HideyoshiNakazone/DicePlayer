#!/export/apps/python/361/bin/python3

import os, sys, time, signal
import argparse
import shutil
from multiprocessing import Process, connection

import DPpack.Dice as Dice
import DPpack.Gaussian as Gaussian
from DPpack.PTable import *
from DPpack.SetGlobals import *
from DPpack.MolHandling import *
from DPpack.Misc import *


if __name__ == '__main__':
####  Read and store the arguments passed to the program  ####
####  and set the usage and help messages                 ####

	parser = argparse.ArgumentParser(prog='Diceplayer')
	parser.add_argument('--continue', dest='opt_continue' , default=False, action='store_true')
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	parser.add_argument('-i', dest='infile', default='control.in', metavar='INFILE', 
    			                   help='input file of diceplayer [default = control.in]')
	parser.add_argument('-o', dest='outfile', default='run.log', metavar='OUTFILE', 
    				                 help='output file of diceplayer [default = run.log]')
	## Study the option of a parameter for continuing the last process via data from control.in and run.log files

	args = parser.parse_args()

####  Open OUTFILE for writing and print keywords and initial info

	try:

		if args.opt_continue and os.path.exists(args.outfile):
			
			outfile = open(args.outfile,'r')
			run_file = outfile.readlines()
			control_sequence = '                                         Step # '
			
			for line in run_file:
				if control_sequence in line:
					cyc = int(line[-2]) + 1

			outfile.close()
			os.rename(os.path.abspath(args.outfile),os.path.abspath(args.outfile)+".backup")
			outfile = open(args.outfile,'w')


		if os.path.exists(args.outfile):		
			os.rename(os.path.abspath(args.outfile),os.path.abspath(args.outfile)+".backup")
			outfile = open(args.outfile,'w')
		else:
			outfile = open(args.outfile,"w")

	except EnvironmentError as err:
		sys.exit(err)

	try:

		if os.path.exists(args.infile):
			infile = open(args.infile,"r")

	except EnvironmentError as err:
		sys.exit(err)

####  Read and check the keywords in INFILE

	internal = Internal(infile, outfile)

	internal.read_keywords()

	if args.opt_continue:
		internal.player.cyc = cyc

	internal.check_keywords()
	internal.print_keywords()

# ####  Check whether the executables are in the path

	# internal.check_executables()

# ####  Read the potential, store the info in 'molecules' and prints the info in OUTFILE

	internal.read_potential()

	# if internal.player.lps == "yes":
	# 	read_lps()

	# if internal.player.ghosts == "yes":
	# 	read_ghosts()

	internal.print_potential()

####  Bring the molecules to standard orientation and prints info about them

	for i in range(len(internal.system.molecule)):
		internal.outfile.write("\nMolecule type {}:\n\n".format(i + 1))
		internal.system.molecule[i].print_mol_info(internal.outfile)
		internal.outfile.write("    Translating and rotating molecule to standard orientation...")
		internal.system.molecule[i].standard_orientation()
		internal.outfile.write(" Done\n\n    New values:\n")
		internal.system.molecule[i].print_mol_info(internal.outfile)
	
	internal.outfile.write(90 * "=")
	internal.outfile.write("\n")

####  Open the geoms.xyz file and prints the initial geometry if starting from zero

	if internal.player.cyc == 1:
		try:
			geomsfh = open("geoms.xyz", "w", 1)
		except EnvironmentError as err:
			sys.exit(err)
		internal.system.print_geom(0, geomsfh)
	else:
		try:
			geomsfh = open("geoms.xyz", "A", 1)
		except EnvironmentError as err:
			sys.exit(err)

	# internal.outfile.write("\nStarting the iterative process.\n")
	
	# ## Initial position (in Bohr)
	# position = internal.system.molecule[0].read_position()
	
	# ## If restarting, read the last gradient and hessian
	# if internal.player.cyc > 1:
	# 	if internal.player.qmprog in ("g03", "g09", "g16"):
	# 		Gaussian.read_forces("grad_hessian.dat")
	# 		Gaussian.read_hessian_fchk("grad_hessian.dat")
	
	# 	#if player['qmprog'] == "molcas":
	# 		#Molcas.read_forces("grad_hessian.dat")
	# 		#Molcas.read_hessian("grad_hessian.dat")
	
	# ####
	# ####  Start the iterative process
	# ####
	
# 	for cycle in range(internal.player.cyc, internal.player.cyc + internal.player.maxcyc):
	
# 		internal.outfile.write("\n" + 90 * "-" + "\n")
# 		internal.outfile.write("{} Step # {}\n".format(40 * " ", cycle))
# 		internal.outfile.write(90 * "-" + "\n\n")
	
# 		make_step_dir(cycle)
	
# 		if internal.player.altsteps == 0 or cycle == 1:
# 			internal.dice.randominit = True
# 		else:
# 			internal.dice.randominit = False
	
# 		####
# 		####  Start block of parallel simulations
# 		####
		
# 		procs = []
# 		sentinels = []
# 		for proc in range(1, internal.player.nprocs + 1):
		
# 			p = Process(target=Dice.simulation_process, args=(cycle, proc, internal.outfile))
# 			p.start()
# 			procs.append(p)
# 			sentinels.append(p.sentinel)
			
# 		while procs:
# 			finished = connection.wait(sentinels)
# 			for proc_sentinel in finished:
# 				i = sentinels.index(proc_sentinel)
# 				status = procs[i].exitcode
# 				procs.pop(i)
# 				sentinels.pop(i)
# 				if status != 0:
# 					for p in procs:
# 						p.terminate()
# 					sys.exit(status)
		
# 		for proc in range(1, internal.player.nprocs + 1):
# 			Dice.print_last_config(cycle, proc)
		
# 		####
# 		####  End of parallel simulations block
# 		####	
		
# 		## Make ASEC
# 		internal.outfile.write("\nBuilding the ASEC and vdW meanfields... ")
# 		asec_charges = internal.populate_asec_vdw(cycle)
		
# 		## After ASEC is built, compress files bigger than 1MB
# 		for proc in range(1, internal.player.nprocs + 1):
# 			path = "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
# 			compress_files_1mb(path)
		
# 		####
# 		####  Start QM calculation
# 		####
		
# 		make_qm_dir(cycle)
		
# 		if internal.player.opt == "yes":
			
# 			##
# 			##  Gaussian block
# 			##
# 			if internal.player.qmprog in ("g03", "g09", "g16"):
				
# 				if cycle > 1:
# 					src = "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
# 					dst = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
# 					shutil.copyfile(src, dst)
				
# 				Gaussian.make_force_input(cycle, asec_charges)
# 				Gaussian.run_gaussian(cycle, "force", internal.outfile)
# 				Gaussian.run_formchk(cycle, internal.outfile)
				
# 				## Read the gradient
# 				file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.fchk"
# 				gradient = Gaussian.read_forces(file, internal.outfile)
# 				if len(cur_gradient) > 0:
# 					old_gradient = cur_gradient
# 				cur_gradient = gradient
				
# 				## If 1st step, read the hessian
# 				if cycle == 1:
# 					if internal.player.readhessian == "yes":
# 						file = "grad_hessian.dat"
# 						internal.outfile.write("\nReading the hessian matrix from file {}\n".format(file))
# 						hessian = Gaussian.read_hessian_fchk(file)
# 					else:
# 						file = "step01" + os.sep + "qm" + os.sep + "asec.fchk"
# 						internal.outfile.write("\nReading the hessian matrix from file {}\n".format(file))
# 						hessian = internal.gaussian.read_hessian(file)
				
# 				## From 2nd step on, update the hessian
# 				else:
# 					internal.outfile.write("\nUpdating the hessian matrix using the BFGS method... ")
# 					hessian = internal.system.molecule[0].update_hessian(step, cur_gradient, old_gradient, hessian)
# 					internal.outfile.write("Done\n")
				
# 				## Save gradient and hessian
# 				internal.gaussian.print_grad_hessian(cycle, cur_gradient, hessian)
				
# 				## Calculate the step and update the position
# 				step = internal.calculate_step(cur_gradient, hessian, internal.outfile)
# 				position += step
				
# 				## Update the geometry of the reference molecule
# 				internal.system.molecule[0].update_molecule(position, internal.outfile)
				
# 				## If needed, calculate the charges
# 				if cycle < internal.player.switchcyc:
					
# 					internal.gaussian.make_charge_input(cycle, asec_charges)
# 					internal.gaussian.run_gaussian(cycle, "charge", internal.outfile)
				
# 				## Read the new charges and update molecules[0]
# 				if cycle < internal.player.switchcyc:
# 					file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
# 					internal.gaussian.read_charges(file, internal.outfile)
# 				else:
# 					file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.log"
# 					internal.gaussian.read_charges(file, internal.outfile)
				
# 				## Print new info for molecule[0]
# 				internal.outfile.write("\nNew values for molecule type 1:\n\n")
# 				internal.system.molecule[0].print_mol_info()
				
# 				## Print new geometry in geoms.xyz
# 				internal.system.molecule[0].print_geom(cycle, geomsfh)
				
# 			##
# 			##  Molcas block
# 			##
# 			#if player['qmprog'] == "molcas":
		
		
# 		#elif player['opt'] == "ts":
			
# 			##
# 			##  Gaussian block
# 			##
# 			#if player['qmprog'] in ("g03", "g09", "g16"):
				
				
				
# 			##
# 			##  Molcas block
# 			##
# 			#if player['qmprog'] == "molcas":
				
		
# 		else:    ## Only relax the charge distribution
			
# 			if internal.player.qmprog in ("g03", "g09", "g16"):
				
# 				if cycle > 1:
# 					src = "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
# 					dst = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
# 					shutil.copyfile(src, dst)
				
# 				Gaussian.make_charge_input(cycle, asec_charges)
# 				Gaussian.run_gaussian(cycle, "charge", internal.outfile)
				
# 				file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
# 				Gaussian.read_charges(file)
				
# 				## Print new info for molecule[0]
# 				internal.outfile.write("\nNew values for molecule type 1:\n\n")
# 				internal.system.molecule[0].print_mol_info()
				
# 			#if player['qmprog'] == "molcas":
				
		
	
# 	####
# 	####  End of the iterative process
# 	####

# ## imprimir ultimas mensagens, criar um arquivo de potencial para ser usado em eventual
# ## continuacao, fechar arquivos (geoms.xyz, run.log, ...)

# 	internal.outfile.write("\nDiceplayer finished normally!\n")
# 	internal.outfile.close()
# ####
# ####  End of the program
# ####