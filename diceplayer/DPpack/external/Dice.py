from diceplayer.DPpack.MolHandling import *
from diceplayer.DPpack.PTable import *
from diceplayer.DPpack.Misc import *


from typing import IO, Tuple, List, TextIO, Union

from numpy.core.numeric import partition
from numpy import random

from multiprocessing import Process, connection
import setproctitle
import subprocess
import os
import sys
import shutil
import textwrap
import types


dice_end_flag = "End of simulation"  # The normal end flag
dice_flag_line = -2  # must be in the line before the last
umaAng3_to_gcm3 = 1.6605  # Conversion between uma/Ang3 to g/cm3

max_seed = 4294967295  # Maximum allowed value for a seed (numpy)


class Dice:

	def __init__(self, infile: TextIO, outfile: TextIO) -> None:

		self.title = "Diceplayer run"
		self.progname = "dice"

		self.infile = infile
		self.outfile = outfile

		self.path = None
		self.nprocs: int = None

		self.randominit = 'first'
		self.temp = 300.0
		self.press = 1.0
		self.isave = 1000         # ASEC construction will take this into account
		self.ncores = 1

		# Investigate the possibility of using 'box = Lx Ly Lz' instead.
		self.dens = None
		# self.box = None		# So 'geom' would be set by diceplayer and 'cutoff' would be
		# switched off. One of them must be given.
		self.combrule = "*"
		self.ljname = None
		self.outname = None
		# Up to 4 integer values related to up to 4 molecule types
		self.nmol: List[int] = []
		# 2 or 3 integer values related to 2 or 3 simulations
		self.nstep: List[int] = []
		# (NVT th + NVT eq) or (NVT th + NPT th + NPT eq).
		# This will control the 'nstep' keyword of Dice
		self.upbuf = 360

	def __new_density(self, cycle: int, proc: int) -> float:

		sim_dir = "simfiles"
		step_dir = "step{:02d}".format(cycle-1)
		proc_dir = "p{:02d}".format(proc)
		path = sim_dir + os.sep + step_dir + os.sep + proc_dir
		file = path + os.sep + "last.xyz"

		if not os.path.isfile(file):
			sys.exit(
				"Error: cannot find the xyz file {} in main directory".format(file))
		try:
			with open(file) as fh:
				xyzfile = fh.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(file))

		box = xyzfile[1].split()
		volume = float(box[-3]) * float(box[-2]) * float(box[-1])

		total_mass = 0
		for i in range(len(self.molecule)):

			total_mass += self.molecule[i].total_mass * self.nmol[i]

		density = (total_mass / volume) * umaAng3_to_gcm3

		return density

	def __print_last_config(self, cycle: int, proc: int) -> None:
		
		sim_dir = "simfiles"
		step_dir = "step{:02d}".format(cycle)
		proc_dir = "p{:02d}".format(proc)
		path = sim_dir + os.sep + step_dir + os.sep + proc_dir
		file = path + os.sep + "phb.xyz"
		if not os.path.isfile(file):
			sys.exit("Error: cannot find the xyz file {}".format(file))
		try:
			with open(file) as fh:
				xyzfile = fh.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(file))
		
		nsites = len(self.molecule[0].atom) * self.nmol[0]	
		for i in range(1, len(self.nmol)):
			nsites += self.nmol[i] * len(self.molecule[i].atom)
		
		nsites += 2  ## To include the comment line and the number of atoms (xyz file format)
		
		nsites *= -1  ## Become an index to count from the end of xyzfile (list)
		xyzfile = xyzfile[nsites :]  ## Take the last configuration
		
		
		file = path + os.sep + "last.xyz"
		fh = open(file, "w")
		for line in xyzfile:
			fh.write(line)

	def __make_dice_inputs(self, cycle: int, proc: int) -> None:

		sim_dir = "simfiles"
		step_dir = "step{:02d}".format(cycle)
		proc_dir = "p{:02d}".format(proc)
		path = sim_dir + os.sep + step_dir + os.sep + proc_dir

		num = time.time()  # Take the decimal places 7 to 12 of the
		num = (num - int(num)) * 1e6  # time in seconds as a floating point
		# to make an integer in the range 1-1e6
		num = int((num - int(num)) * 1e6)
		random.seed((os.getpid() * num) % (max_seed + 1))

		if (self.randominit == 'first' and cycle > self.initcyc):
			last_step_dir = "step{:02d}".format(cycle-1)
			last_path = sim_dir + os.sep + last_step_dir + os.sep + proc_dir
			xyzfile = last_path + os.sep + "last.xyz"
			self.__make_init_file(path, xyzfile)

		if len(self.nstep) == 2:  # Means NVT simulation

			self.__make_nvt_ter(cycle, path)
			self.__make_nvt_eq(path)

		elif len(self.nstep) == 3:  # Means NPT simulation

			if (self.randominit == 'first' and cycle > self.initcyc):
				self.dens = self.__new_density(cycle, proc)
			else:
				self.__make_nvt_ter(cycle, path)

			self.__make_npt_ter(cycle, path)
			self.__make_npt_eq(path)

		else:
			sys.exit("Error: bad number of entries for 'nstep'")

		self.__make_potential(path)

		# if (self.randominit == 'first' and cycle > self.initcyc):

		# 	last_path = sim_dir + os.sep + "step{:02d}".format(cycle-1) + os.sep + proc_dir
		# 	shutil.copyfile(last_path + os.sep + "phb.dat", path + os.sep + "phb.dat")

	def __make_nvt_ter(self, cycle: int, path: str) -> None:

		file = path + os.sep + "NVT.ter"
		try:
			fh = open(file, "w")
		except:
			sys.exit("Error: cannot open file {}".format(file))

		fh.write("title = {} - NVT Thermalization\n".format(self.title))
		fh.write("ncores = {}\n".format(self.ncores))
		fh.write("ljname = {}\n".format(self.ljname))
		fh.write("outname = {}\n".format(self.outname))

		string = " ".join(str(x) for x in self.nmol)
		fh.write("nmol = {}\n".format(string))

		fh.write("dens = {}\n".format(self.dens))
		fh.write("temp = {}\n".format(self.temp))

		if (self.randominit == 'first' and cycle > self.initcyc):
			fh.write("init = yesreadxyz\n")
			fh.write("nstep = {}\n".format(self.altsteps))
		else:
			fh.write("init = yes\n")
			fh.write("nstep = {}\n".format(self.nstep[0]))

		fh.write("vstep = 0\n")
		fh.write("mstop = 1\n")
		fh.write("accum = no\n")
		fh.write("iprint = 1\n")
		fh.write("isave = 0\n")
		fh.write("irdf = 0\n")

		seed = int(1e6 * random.random())
		fh.write("seed = {}\n".format(seed))
		fh.write("upbuf = {}".format(self.upbuf))

		fh.close()

	def __make_nvt_eq(self, path: str) -> None:

		file = path + os.sep + "NVT.eq"
		try:
			fh = open(file, "w")
		except:
			sys.exit("Error: cannot open file {}".format(file))

		fh.write("title = {} - NVT Production\n".format(self.title))
		fh.write("ncores = {}\n".format(self.ncores))
		fh.write("ljname = {}\n".format(self.ljname))
		fh.write("outname = {}\n".format(self.outname))

		string = " ".join(str(x) for x in self.nmol)
		fh.write("nmol = {}\n".format(string))

		fh.write("dens = {}\n".format(self.dens))
		fh.write("temp = {}\n".format(self.temp))
		fh.write("init = no\n")
		fh.write("nstep = {}\n".format(self.nstep[1]))
		fh.write("vstep = 0\n")
		fh.write("mstop = 1\n")
		fh.write("accum = no\n")
		fh.write("iprint = 1\n")
		fh.write("isave = {}\n".format(self.isave))
		fh.write("irdf = {}\n".format(10 * self.nprocs))

		seed = int(1e6 * random.random())
		fh.write("seed = {}\n".format(seed))

		fh.close()

	def __make_npt_ter(self, cycle: int, path: str) -> None:

		file = path + os.sep + "NPT.ter"
		try:
			fh = open(file, "w")
		except:
			sys.exit("Error: cannot open file {}".format(file))

		fh.write("title = {} - NPT Thermalization\n".format(self.title))
		fh.write("ncores = {}\n".format(self.ncores))
		fh.write("ljname = {}\n".format(self.ljname))
		fh.write("outname = {}\n".format(self.outname))

		string = " ".join(str(x) for x in self.nmol)
		fh.write("nmol = {}\n".format(string))

		fh.write("press = {}\n".format(self.press))
		fh.write("temp = {}\n".format(self.temp))

		if (self.randominit == 'first' and cycle > self.initcyc):
			fh.write("init = yesreadxyz\n")
			fh.write("dens = {:<8.4f}\n".format(self.dens))
			fh.write("vstep = {}\n".format(int(self.altsteps / 5)))
		else:
			# Because there will be a previous NVT simulation
			fh.write("init = no\n")
			fh.write("vstep = {}\n".format(int(self.nstep[1] / 5)))

		fh.write("nstep = 5\n")
		fh.write("mstop = 1\n")
		fh.write("accum = no\n")
		fh.write("iprint = 1\n")
		fh.write("isave = 0\n")
		fh.write("irdf = 0\n")

		seed = int(1e6 * random.random())
		fh.write("seed = {}\n".format(seed))

		fh.close()

	def __make_npt_eq(self, path: str) -> None:

		file = path + os.sep + "NPT.eq"
		try:
			fh = open(file, "w")
		except:
			sys.exit("Error: cannot open file {}".format(file))

		fh.write("title = {} - NPT Production\n".format(self.title))
		fh.write("ncores = {}\n".format(self.ncores))
		fh.write("ljname = {}\n".format(self.ljname))
		fh.write("outname = {}\n".format(self.outname))

		string = " ".join(str(x) for x in self.nmol)
		fh.write("nmol = {}\n".format(string))

		fh.write("press = {}\n".format(self.press))
		fh.write("temp = {}\n".format(self.temp))

		fh.write("nstep = 5\n")

		fh.write("vstep = {}\n".format(int(self.nstep[2] / 5)))
		fh.write("init = no\n")
		fh.write("mstop = 1\n")
		fh.write("accum = no\n")
		fh.write("iprint = 1\n")
		fh.write("isave = {}\n".format(self.isave))
		fh.write("irdf = {}\n".format(10 * self.nprocs))

		seed = int(1e6 * random.random())
		fh.write("seed = {}\n".format(seed))

		fh.close()

	def __make_init_file(self, path: str, file: TextIO) -> None:

		if not os.path.isfile(file):
			sys.exit(
				"Error: cannot find the xyz file {} in main directory".format(file))
		try:
			with open(file) as fh:
				xyzfile = fh.readlines()
		except:
			sys.exit("Error: cannot open file {}".format(file))

		nsites_mm = 0
		for i in range(1, len(self.nmol)):
			nsites_mm += self.nmol[i] * len(self.molecule[i].atom)

		# Become an index to count from the end of xyzfile (list)
		nsites_mm *= -1
		# Only the MM atoms of the last configuration remains
		xyzfile = xyzfile[nsites_mm:]

		file = path + os.sep + self.outname + ".xy"

		try:
			fh = open(file, "w", 1)
		except:
			sys.exit("Error: cannot open file {}".format(file))

		for atom in self.molecule[0].atom:
			fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(
				atom.rx, atom.ry, atom.rz))

		# for i in self.molecule[0].ghost_atoms:
		# 	with self.molecule[0].atom[i] as ghost:
		# 		fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(ghost.rx, ghost.ry, ghost.rz))

		# for i in self.molecule[0].lp_atoms:
		# 	with self.molecule[0].atom[i] as lp:
		# 		fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(lp.rx, lp.ry, lp.rz))

		for line in xyzfile:
			atom = line.split()
			rx = float(atom[1])
			ry = float(atom[2])
			rz = float(atom[3])
			fh.write("{:>10.6f}  {:>10.6f}  {:>10.6f}\n".format(rx, ry, rz))

		fh.write("$end")

		fh.close()

	def __make_potential(self, path: str) -> None:

		fstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f}\n"

		file = path + os.sep + self.ljname
		try:
			fh = open(file, "w")
		except:
			sys.exit("Error: cannot open file {}".format(file))

		fh.write("{}\n".format(self.combrule))
		fh.write("{}\n".format(len(self.nmol)))

		nsites_qm = len(self.molecule[0].atom) + len(
			self.molecule[0].ghost_atoms) + len(self.molecule[0].lp_atoms)

		# Print the sites of the QM self.molecule
		fh.write("{} {}\n".format(nsites_qm, self.molecule[0].molname))
		for atom in self.molecule[0].atom:
			fh.write(fstr.format(atom.lbl, atom.na, atom.rx, atom.ry, atom.rz,
								 atom.chg, atom.eps, atom.sig))

		ghost_label = self.molecule[0].atom[-1].lbl + 1
		for i in self.molecule[0].ghost_atoms:
			fh.write(fstr.format(ghost_label, ghost_number, self.molecule[0].atom[i].rx, self.molecule[0].atom[i].ry,
								 self.molecule[0].atom[i].rz, self.molecule[0].atom[i].chg, 0, 0))

		ghost_label += 1
		for lp in self.molecule[0].lp_atoms:
			fh.write(fstr.format(ghost_label, ghost_number, lp['rx'], lp['ry'], lp['rz'],
								 lp['chg'], 0, 0))

		# Print the sites of the other self.molecules
		for mol in self.molecule[1:]:
			fh.write("{} {}\n".format(len(mol.atom), mol.molname))
			for atom in mol.atom:
				fh.write(fstr.format(atom.lbl, atom.na, atom.rx, atom.ry,
									 atom.rz, atom.chg, atom.eps, atom.sig))

	def __make_proc_dir(self, cycle: int, proc: int) -> None:

		sim_dir = "simfiles"
		step_dir = "step{:02d}".format(cycle)
		proc_dir = "p{:02d}".format(proc)
		path = sim_dir + os.sep + step_dir + os.sep + proc_dir
		try:
			os.makedirs(path)
		except:
			sys.exit("Error: cannot make directory {}".format(path))

	def __run_dice(self, cycle: int, proc: int, fh: str) -> None:

		sim_dir = "simfiles"
		step_dir = "step{:02d}".format(cycle)
		proc_dir = "p{:02d}".format(proc)

		try:
			fh.write("Simulation process {} initiated with pid {}\n".format(
				sim_dir + os.sep + step_dir + os.sep + proc_dir, os.getpid()))

		except Exception as err:
			print("I/O error({0}): {1}".format(err))

		path = sim_dir + os.sep + step_dir + os.sep + proc_dir
		working_dir = os.getcwd()
		os.chdir(path)

		if len(self.nstep) == 2:  # Means NVT simulation

			if (self.randominit == 'first' and cycle > self.initcyc):
				string_tmp = 'previous'
			else:
				string_tmp = 'random'

			# NVT thermalization
			string = "(from " + string_tmp + " configuration)"
			fh.write("p{:02d}> NVT thermalization finished {} on {}\n".format(proc, string,
																			  date_time()))

			infh = open("NVT.ter")
			outfh = open("NVT.ter.out", "w")

			if shutil.which("bash") != None:
				exit_status = subprocess.call(
					["bash", "-c", "exec -a dice-step{}-p{} {} < {} > {}".format(cycle, proc, self.progname, infh.name, outfh.name)])
			else:
				exit_status = subprocess.call(
					self.progname, stin=infh.name, stout=outfh.name)

			infh.close()
			outfh.close()

			if os.getppid() == 1:  # Parent process is dead
				sys.exit()

			if exit_status != 0:
				sys.exit(
					"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))
			else:
				# Open again to seek the normal end flag
				outfh = open("NVT.ter.out")
				flag = outfh.readlines()[dice_flag_line].strip()
				outfh.close()
				if flag != dice_end_flag:
					sys.exit(
						"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))

			# NVT production
			fh.write("p{:02d}> NVT production initiated on {}\n".format(
				proc, date_time()))

			infh = open("NVT.eq")
			outfh = open("NVT.eq.out", "w")

			if shutil.which("bash") != None:
				exit_status = subprocess.call(
					["bash", "-c", "exec -a dice-step{}-p{} {} < {} > {}".format(cycle, proc, self.progname, infh.name, outfh.name)])
			else:
				exit_status = subprocess.call(
					self.progname, stin=infh.name, stout=outfh.name)

			infh.close()
			outfh.close()

			if os.getppid() == 1:  # Parent process is dead
				sys.exit()

			if exit_status != 0:
				sys.exit(
					"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))
			else:
				# Open again to seek the normal end flag
				outfh = open("NVT.eq.out")
				flag = outfh.readlines()[dice_flag_line].strip()
				outfh.close()
				if flag != dice_end_flag:
					sys.exit(
						"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))

			fh.write("p{:02d}> ----- NVT production finished on {}\n".format(proc,
																			 date_time()))

		elif len(self.nstep) == 3:  # Means NPT simulation

			# NVT thermalization if randominit
			if self.randominit == 'always' or (self.randominit == 'first' and cycle == 1) or self.continued:
				string = "(from random configuration)"
				fh.write("p{:02d}> NVT thermalization initiated {} on {}\n".format(proc,
																				   string, date_time()))
				infh = open("NVT.ter")
				outfh = open("NVT.ter.out", "w")

				if shutil.which("bash") != None:
					exit_status = subprocess.call(
						["bash", "-c", "exec -a dice-step{}-p{} {} < {} > {}".format(cycle, proc, self.progname, infh.name, outfh.name)])
				else:
					exit_status = subprocess.call(
						self.progname, stin=infh.name, stout=outfh.name)

				infh.close()
				outfh.close()

				if os.getppid() == 1:  # Parent process is dead
					sys.exit()

				if exit_status != 0:
					sys.exit(
						"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))
				else:
					# Open again to seek the normal end flag
					outfh = open("NVT.ter.out")
					flag = outfh.readlines()[dice_flag_line].strip()
					outfh.close()
					if flag != dice_end_flag:
						sys.exit(
							"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))

			# NPT thermalization
			if not self.randominit == 'always' or ((self.randominit == 'first' and cycle > self.initcyc)):
				string = " (from previous configuration) "
			else:
				string = " "
			fh.write("p{:02d}> NPT thermalization finished {} on {}\n".format(proc, string,
																			  date_time()))

			infh = open("NPT.ter")
			outfh = open("NPT.ter.out", "w")

			if shutil.which("bash") != None:
				exit_status = subprocess.call(
					["bash", "-c", "exec -a dice-step{}-p{} {} < {} > {}".format(cycle, proc, self.progname, infh.name, outfh.name)])
			else:
				exit_status = subprocess.call(
					self.progname, stin=infh.name, stout=outfh.name)

			infh.close()
			outfh.close()

			if os.getppid() == 1:  # Parent process is dead
				sys.exit()

			if exit_status != 0:
				sys.exit(
					"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))
			else:
				# Open again to seek the normal end flag
				outfh = open("NPT.ter.out")
				flag = outfh.readlines()[dice_flag_line].strip()
				outfh.close()
				if flag != dice_end_flag:
					sys.exit(
						"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))

			# NPT production
			fh.write("p{:02d}> NPT production initiated on {}\n".format(
				proc, date_time()))

			infh = open("NPT.eq")
			outfh = open("NPT.eq.out", "w")

			if shutil.which("bash") != None:
				exit_status = subprocess.call(
					["bash", "-c", "exec -a dice-step{}-p{} {} < {} > {}".format(cycle, proc, self.progname, infh.name, outfh.name)])
			else:
				exit_status = subprocess.call(
					self.progname, stin=infh.name, stout=outfh.name)

			infh.close()
			outfh.close()

			if os.getppid() == 1:  # Parent process is dead
				sys.exit()

			if exit_status != 0:
				sys.exit(
					"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))
			else:
				# Open again to seek the normal end flag
				outfh = open("NPT.eq.out")
				flag = outfh.readlines()[dice_flag_line].strip()
				outfh.close()
				if flag != dice_end_flag:
					sys.exit(
						"Dice process step{:02d}-p{:02d} did not exit properly".format(cycle, proc))

			fh.write("p{:02d}> ----- NPT production finished on {}\n".format(proc,
																			 date_time()))

		os.chdir(working_dir)

	def __simulation_process(self, cycle: int, proc: int):
		setproctitle.setproctitle(
			"diceplayer-step{:0d}-p{:0d}".format(cycle, proc))

		try:
			self.__make_proc_dir(cycle, proc)
			self.__make_dice_inputs(cycle, proc)
			self.__run_dice(cycle, proc, self.outfile)
		except Exception as err:
			sys.exit(err)

	def configure(self, initcyc: int, nprocs: int, altsteps: int, nmol: List[int], molecule: List[Molecule]):

		self.initcyc = initcyc

		self.nprocs = nprocs
		self.altsteps = altsteps

		self.molecule = molecule

	def start(self, cycle: int) -> None:

		procs = []
		sentinels = []

		for proc in range(1, self.nprocs + 1):

			p = Process(target=self.__simulation_process, args=(cycle, proc))
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

		for proc in range(1, self.nprocs + 1):
			self.__print_last_config(cycle, proc)

	def reset(self):

		del self.nprocs
		del self.altsteps
		del self.molecule
