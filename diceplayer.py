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
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')
	parser.add_argument('-i', dest='infile', default='control.in', metavar='INFILE', 
    			                   help='input file of diceplayer [default = control.in]')
	parser.add_argument('-o', dest='outfile', default='run.log', metavar='OUTFILE', 
    				                 help='output file of diceplayer [default = run.log]')
	## Study the option of a parameter for continuing the last process via data from control.in and run.log files

	args = parser.parse_args()

####  Read and check the keywords in INFILE

	read_keywords(args.infile)
	check_keywords(args.infile)

####  Open OUTFILE for writing and print keywords and initial info

	try:
		if player['initcyc'] > 1 and os.path.exists(args.outfile):
			oldname = args.outfile + ".old"
			os.replace(args.outfile, oldname)
		logfh = open(args.outfile, 'w', 1)
	except EnvironmentError as err:
		sys.exit(err)

	print_keywords(logfh)

####  Check whether the executables are in the path

	check_executables(logfh)

####  Read the potential, store the info in 'molecules' and prints the info in OUTFILE

	read_potential(args.infile)

	if player['lps'] == "yes":
		read_lps()

	if player['ghosts'] == "yes":
		read_ghosts()

	print_potential(logfh)

####  Bring the molecules to standard orientation and prints info about them

	for i in range(len(molecules)):
		logfh.write("\nMolecule type {}:\n\n".format(i + 1))
		print_mol_info(molecules[i], logfh)
		logfh.write("    Translating and rotating molecule to standard orientation...")
		standard_orientation(molecules[i])
		logfh.write(" Done\n\n    New values:\n")
		print_mol_info(molecules[i], logfh)
	
	logfh.write(90 * "=")
	logfh.write("\n")

####  Open the geoms.xyz file and prints the initial geometry if starting from zero

	if player['initcyc'] == 1:
		try:
			geomsfh = open("geoms.xyz", "w", 1)
		except EnvironmentError as err:
			sys.exit(err)
		print_geom(0, geomsfh)
	else:
		try:
			geomsfh = open("geoms.xyz", "A", 1)
		except EnvironmentError as err:
			sys.exit(err)


	logfh.write("\nStarting the iterative process.\n")
	
	## Initial position (in Bohr)
	position = read_position(molecules[0])
	
	## If restarting, read the last gradient and hessian
	if player['initcyc'] > 1:
		if player['qmprog'] in ("g03", "g09", "g16"):
			Gaussian.read_forces("grad_hessian.dat")
			Gaussian.read_hessian_fchk("grad_hessian.dat")
	
		#if player['qmprog'] == "molcas":
			#Molcas.read_forces("grad_hessian.dat")
			#Molcas.read_hessian("grad_hessian.dat")
	
	####
	####  Start the iterative process
	####
	
	for cycle in range(player['initcyc'], player['initcyc'] + player['maxcyc']):
	
		logfh.write("\n" + 90 * "-" + "\n")
		logfh.write("{} Step # {}\n".format(40 * " ", cycle))
		logfh.write(90 * "-" + "\n\n")
	
		make_step_dir(cycle)
	
		if player['altsteps'] == 0 or cycle == 1:
			dice['randominit'] = True
		else:
			dice['randominit'] = False
	
		####
		####  Start block of parallel simulations
		####
		
		procs = []
		sentinels = []
		for proc in range(1, player['nprocs'] + 1):
		
			p = Process(target=Dice.simulation_process, args=(cycle, proc, logfh))
			p.start()
			procs.append(p)
			sentinels.append(p.sentinel)
			
		while procs:
			finished = connection.wait(sentinels)
			for proc_sentinel in finished:
				i = sentinels.index(proc_sentinel)
				status = procs[i].exitcode
				procs.pop(i)
				sentinels.pop(i)
				if status != 0:
					for p in procs:
						p.terminate()
					sys.exit(status)
		
		for proc in range(1, player['nprocs'] + 1):
			Dice.print_last_config(cycle, proc)
		
		####
		####  End of parallel simulations block
		####	
		
		## Make ASEC
		logfh.write("\nBuilding the ASEC and vdW meanfields... ")
		asec_charges = populate_asec_vdw(cycle, logfh)
		
		## After ASEC is built, compress files bigger than 1MB
		for proc in range(1, player['nprocs'] + 1):
			path = "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
			compress_files_1mb(path)
		
		####
		####  Start QM calculation
		####
		
		make_qm_dir(cycle)
		
		if player['opt'] == "yes":
			
			##
			##  Gaussian block
			##
			if player['qmprog'] in ("g03", "g09", "g16"):
				
				if cycle > 1:
					src = "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
					dst = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
					shutil.copyfile(src, dst)
				
				Gaussian.make_force_input(cycle, asec_charges)
				Gaussian.run_gaussian(cycle, "force", logfh)
				Gaussian.run_formchk(cycle, logfh)
				
				## Read the gradient
				file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.fchk"
				gradient = Gaussian.read_forces(file, logfh)
				if len(cur_gradient) > 0:
					old_gradient = cur_gradient
				cur_gradient = gradient
				
				## If 1st step, read the hessian
				if cycle == 1:
					if player['readhessian'] == "yes":
						file = "grad_hessian.dat"
						logfh.write("\nReading the hessian matrix from file {}\n".format(file))
						hessian = Gaussian.read_hessian_fchk(file)
					else:
						file = "step01" + os.sep + "qm" + os.sep + "asec.fchk"
						logfh.write("\nReading the hessian matrix from file {}\n".format(file))
						hessian = Gaussian.read_hessian(file)
				
				## From 2nd step on, update the hessian
				else:
					logfh.write("\nUpdating the hessian matrix using the BFGS method... ")
					hessian = update_hessian(step, cur_gradient, old_gradient, hessian)
					logfh.write("Done\n")
				
				## Save gradient and hessian
				Gaussian.print_grad_hessian(cycle, cur_gradient, hessian)
				
				## Calculate the step and update the position
				step = calculate_step(cur_gradient, hessian, logfh)
				position += step
				
				## Update the geometry of the reference molecule
				update_molecule(position, logfh)
				
				## If needed, calculate the charges
				if cycle < player['switchcyc']:
					
					Gaussian.make_charge_input(cycle, asec_charges)
					Gaussian.run_gaussian(cycle, "charge", logfh)
				
				## Read the new charges and update molecules[0]
				if cycle < player['switchcyc']:
					file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
					Gaussian.read_charges(file, logfh)
				else:
					file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.log"
					Gaussian.read_charges(file, logfh)
				
				## Print new info for molecule[0]
				logfh.write("\nNew values for molecule type 1:\n\n")
				print_mol_info(molecules[0], logfh)
				
				## Print new geometry in geoms.xyz
				print_geom(cycle, geomsfh)
				
			##
			##  Molcas block
			##
			#if player['qmprog'] == "molcas":
		
		
		#elif player['opt'] == "ts":
			
			##
			##  Gaussian block
			##
			#if player['qmprog'] in ("g03", "g09", "g16"):
				
				
				
			##
			##  Molcas block
			##
			#if player['qmprog'] == "molcas":
				
		
		else:    ## Only relax the charge distribution
			
			if player['qmprog'] in ("g03", "g09", "g16"):
				
				if cycle > 1:
					src = "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
					dst = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
					shutil.copyfile(src, dst)
				
				Gaussian.make_charge_input(cycle, asec_charges)
				Gaussian.run_gaussian(cycle, "charge", logfh)
				
				file = "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
				Gaussian.read_charges(file)
				
				## Print new info for molecule[0]
				logfh.write("\nNew values for molecule type 1:\n\n")
				print_mol_info(molecules[0], logfh)
				
			#if player['qmprog'] == "molcas":
				
		
	
	####
	####  End of the iterative process
	####

## imprimir ultimas mensagens, criar um arquivo de potencial para ser usado em eventual
## continuacao, fechar arquivos (geoms.xyz, run.log, ...)

	logfh.write("\nDiceplayer finished normally!\n")
	logfh.close()
####
####  End of the program
####