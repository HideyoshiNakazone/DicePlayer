#!/usr/bin/python3

from multiprocessing import Process, connection
import os, sys, time, signal
import setproctitle
import numpy as np
import argparse
import shutil
import pickle

from diceplayer.DPpack.PTable import *
from diceplayer.DPpack.SetGlobals import *
from diceplayer.DPpack.MolHandling import *
from diceplayer.DPpack.Misc import *

__version = 'dev'
setproctitle.setproctitle("diceplayer-{}".format(__version))

if __name__ == '__main__':
####  Read and store the arguments passed to the program  ####
####  and set the usage and help messages                 ####

	parser = argparse.ArgumentParser(prog='Diceplayer')
	parser.add_argument('--continue', dest='opt_continue' , default=False, action='store_true')
	parser.add_argument('--version', action='version', version= "diceplayer-"+__version)
	parser.add_argument('-i', dest='infile', default='control.in', metavar='INFILE', 
								   help='input file of diceplayer [default = control.in]')
	parser.add_argument('-o', dest='outfile', default='run.log', metavar='OUTFILE', 
									 help='output file of diceplayer [default = run.log]')
	## Study the option of a parameter for continuing the last process via data from control.in and run.log files

	args = parser.parse_args()

####  Open OUTFILE for writing and print keywords and initial info

	try:

		if args.opt_continue and os.path.exists(args.outfile):
			
			save = pickle.load(open("latest-step.pkl","rb"))

			if os.path.isfile(args.outfile+".backup"):
				os.remove(args.outfile+".backup")

			os.rename(args.outfile,args.outfile+".backup")
			outfile = open(args.outfile,'w',1)

		elif os.path.exists(args.outfile):		
			os.rename(args.outfile, args.outfile+".backup")
			outfile = open(args.outfile,'w',1)
		else:
			outfile = open(args.outfile,"w",1)

	except Exception as err:
		sys.exit(err)

	try:

		if os.path.exists(args.infile):
			infile = open(args.infile,"r")

	except Exception as err:
		sys.exit(err)

####  Read and check the keywords in INFILE

	internal = Internal(infile, outfile)

	internal.read_keywords()

	internal.check_keywords()
	internal.print_keywords()

	if args.opt_continue:
		internal.player.initcyc = save[0] + 1
		internal.system = save[1]
	else:
		internal.read_potential()

####  Check whether the executables are in the path
####		and print potential to Log File

	internal.check_executables()

	internal.print_potential()

####  Bring the molecules to standard orientation and prints info about them

	for i in range(len(internal.system.molecule)):
		
		internal.outfile.write("\nMolecule type {} - {}:\n\n".format(i + 1, internal.system.molecule[i].molname))
		internal.system.molecule[i].print_mol_info(internal.outfile)
		internal.outfile.write("    Translating and rotating molecule to standard orientation...")
		internal.system.molecule[i].standard_orientation()
		internal.outfile.write(" Done\n\n    New values:\n")
		internal.system.molecule[i].print_mol_info(internal.outfile)
	
	internal.outfile.write(90 * "=")
	internal.outfile.write("\n")

	if not args.opt_continue:
		make_simulation_dir()
	else:
		simdir = "simfiles"
		stepdir = "step{:02d}".format(internal.player.initcyc)
		if os.path.exists(simdir+os.sep+stepdir):
			shutil.rmtree(simdir+os.sep+stepdir)

####  Open the geoms.xyz file and prints the initial geometry if starting from zero

	if internal.player.initcyc == 1:
		try:
			path = "geoms.xyz"
			geomsfh = open(path, "w", 1)
		except EnvironmentError as err:
			sys.exit(err)
		internal.system.print_geom(0, geomsfh)
		geomsfh.write(40 * "-"+"\n")
	else:
		try:
			path = "geoms.xyz"
			geomsfh = open(path, "a", 1)
		except EnvironmentError as err:
			sys.exit(err)

	internal.outfile.write("\nStarting the iterative process.\n")
	
	## Initial position (in Bohr)
	position = internal.system.molecule[0].read_position()
	
	## If restarting, read the last gradient and hessian
	# if internal.player.initcyc > 1:
	# 	if internal.player.qmprog in ("g03", "g09", "g16"):
	# 		Gaussian.read_forces("grad_hessian.dat")
	# 		Gaussian.read_hessian_fchk("grad_hessian.dat")
	
		#if player['qmprog'] == "molcas":
			#Molcas.read_forces("grad_hessian.dat")
			#Molcas.read_hessian("grad_hessian.dat")
	
	###
	###  Start the iterative process
	###

	internal.outfile.write("\n" + 90 * "-" + "\n")
	
	for cycle in range(internal.player.initcyc, internal.player.initcyc + internal.player.maxcyc):
	
		internal.outfile.write("{} Step # {}\n".format(40 * " ", cycle))
		internal.outfile.write(90 * "-" + "\n\n")

		make_step_dir(cycle)
	
		####
		####  Start block of parallel simulations
		####
		
		internal.dice_start(cycle)
		
		###
		###  End of parallel simulations block
		###	
		
		## Make ASEC
		# internal.outfile.write("\nBuilding the ASEC and vdW meanfields... ")
		# asec_charges = internal.populate_asec_vdw(cycle)
		
		# ## After ASEC is built, compress files bigger than 1MB
		# for proc in range(1, internal.player.nprocs + 1):
		# 	path = "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
		# 	compress_files_1mb(path)
		
		###
		###  Start QM calculation
		###
		
		# make_qm_dir(cycle)

		# if internal.player.qmprog in ("g03", "g09", "g16"):
			
		# 	if cycle > 1:

		# 		src = "simfiles" + os.sep + "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
		# 		dst = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
		# 		shutil.copyfile(src, dst)
			
		# 	internal.make_gaussian_input(cycle)
		# 	internal.gaussian.run_gaussian(cycle, "force", internal.outfile)
		# 	internal.gaussian.run_formchk(cycle, internal.outfile)
				
		# 	## Read the gradient
		# 	file = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.fchk"

		# 	try: 
		# 		gradient
		# 		old_gradient = gradient
		# 	except:
		# 		pass

		# 	gradient = internal.read_forces_fchk(file, internal.outfile)

		# 	# If 1st step, read the hessian
		# 	if cycle == 1:

		# 		if internal.player.readhessian == "yes":

		# 			file = "grad_hessian.dat"
		# 			internal.outfile.write("\nReading the hessian matrix from file {}\n".format(file))
		# 			hessian = internal.read_hessian_log(file)

		# 		else:

		# 			file = "simfiles" + os.sep + "step01" + os.sep + "qm" + os.sep + "asec.fchk"
		# 			internal.outfile.write("\nReading the hessian matrix from file {}\n".format(file))
		# 			hessian = internal.read_hessian_fchk(file)

		# 	# From 2nd step on, update the hessian
		# 	else:
		# 		internal.outfile.write("\nUpdating the hessian matrix using the BFGS method... ")
		# 		hessian = internal.system.molecule[0].update_hessian(step, gradient, old_gradient, hessian)
		# 		internal.outfile.write("Done\n")
				
		# 	# Save gradient and hessian
		# 	internal.print_grad_hessian(cycle, gradient, hessian)

		# 	# Calculate the step and update the position
		# 	step = internal.calculate_step(cycle, gradient, hessian)

		# 	position += step
				
		# 	# ## Update the geometry of the reference molecule
		# 	internal.system.update_molecule(position, internal.outfile)
				
# 				## If needed, calculate the charges
# 				if cycle < internal.player.switchcyc:
					
# 					# internal.gaussian.make_charge_input(cycle, asec_charges)
# 					internal.gaussian.run_gaussian(cycle, "charge", internal.outfile)
				
# 				## Read the new charges and update molecules[0]
# 				if cycle < internal.player.switchcyc:
# 					file = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
# 					internal.gaussian.read_charges(file, internal.outfile)
# 				else:
# 					file = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.log"
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
# 					src = "simfiles" + os.sep + "step{:02d}".format(cycle - 1) + os.sep + "qm" + os.sep + "asec.chk"
# 					dst = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec.chk"
# 					shutil.copyfile(src, dst)
				
# 				# internal.gaussian.make_charge_input(cycle, asec_charges)
# 				internal.gaussian.run_gaussian(cycle, "charge", internal.outfile)
				
# 				file = "simfiles" + os.sep + "step{:02d}".format(cycle) + os.sep + "qm" + os.sep + "asec2.log"
# 				internal.read_charges(file)
				
# 				## Print new info for molecule[0]
# 				internal.outfile.write("\nNew values for molecule type 1:\n\n")
# 				internal.system.molecule[0].print_mol_info()
				
# 			#if player['qmprog'] == "molcas":
				
		internal.system.print_geom(cycle, geomsfh)
		geomsfh.write(40 * "-"+"\n")

		internal.outfile.write("\n+" + 88 * "-" + "+\n")

		pickle.dump([cycle,internal.system], open("latest-step.pkl", "wb"))
	####
	####  End of the iterative process
	####

## imprimir ultimas mensagens, criar um arquivo de potencial para ser usado em eventual
## continuacao, fechar arquivos (geoms.xyz, run.log, ...)

	internal.outfile.write("\nDiceplayer finished normally!\n")
	internal.outfile.close()
####
####  End of the program
####
