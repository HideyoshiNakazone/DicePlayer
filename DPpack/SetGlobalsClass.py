import os, sys
import shutil
import textwrap

from DPpack.MolHandling import *
from DPpack.PTable import *
from DPpack.Misc import *

env = ["OMP_STACKSIZE"]

bohr2ang = 0.52917721092
ang2bohr = 1/bohr2ang

class Internal:

	def __init__(self, infile):

		self.infile = infile
		self.system = System()

		self.player = self.Player()
		self.player_keywords = [a for a in dir(self.player) if not a.startswith('__') and not callable(getattr(self.player, a))]

		self.dice = self.Dice()
		self.dice_keywords = [a for a in dir(self.dice) if not a.startswith('__') and not callable(getattr(self.dice, a))]

		self.gaussian = self.Gaussian()
		self.gaussian_keywords = [a for a in dir(self.gaussian) if not a.startswith('__') and not callable(getattr(self.gaussian, a))]

		# self.molcas = self.Molcas()
		# self.molcas_keywords = [a for a in dir(self.molcas) if not a.startswith('__') and not callable(getattr(self.molcas, a))]

		## Constanst that shall be set for global use

		self.tol_rms_force = 3e-4		# Hartree/Bohr
		self.tol_max_force = 4.5e-4		# Hartree/Bohr
		self.tol_rms_step = 1.2e-3		# Bohr
		self.tol_max_step = 1.8e-3		# Bohr
		self.trust_radius = None

		## Dice:
		self.combrule = None
		self.randominit = None

	def read_keywords(self):

		try:
			with open(self.infile) as fh:
				controlfile = fh.readlines()
		except EnvironmentError:
			sys.exit("Error: cannot open file {}".format(self.infile))
		
		for line in controlfile:
			
			key, value = line.partition("=")[::2]  # Discards the '='
			key = key.strip().lower()
			if key in ('title', 'keywords'):
				value = value.strip()
			else:
				value = value.split()
			
			####  Read the Diceplayer related keywords
			if key in self.player_keywords and len(value) != 0:  ##  'value' is not empty!
				
				if key == 'qmprog' and value[0].lower() in ("g03", "g09", "g16", "molcas"):
					setattr(self.player, key, value[0].lower())
				
				elif key == 'opt' and value[0].lower() in ("yes", "no", "ts"):
					setattr(self.player, key, value[0].lower())
				
				#elif key == 'zipprog' and value[0].lower() in ("zip", "gzip", "bzip"):
					#player[key] = value[0].lower()
				
				elif key in ('lps', 'ghosts') and value[0].lower() in ("yes", "no"):
					setattr(self.player, key, value[0].lower())
				
				elif key in ('readhessian', 'vdwforces') and value[0].lower() in ("yes", "no"):
					setattr(self.player, key, value[0].lower())
				
				elif key in ('maxcyc', 'initcyc', 'nprocs', 'altsteps', 'switchcyc'):
					err = "Error: expected a positive integer for keyword {} in file {}".format(key, self.infile)
					try:
						new_value = int(value[0])
						if new_value >= 1:
							setattr(self.player, key, new_value)
						elif key == 'altsteps' and new_value == 0:
							setattr(self.player, key, 0)
					except ValueError:
						sys.exit(err)
				
				elif key == 'maxstep':  # Cannot be less than 0.01
					err = "Error: expected a float greater than 0.01 for keyword {} in file {}".format(key, self.infile)
					try:
						new_value = float(value[0])
						if new_value < 0.01:
							sys.exit(err)
						else:
							setattr(self.player, key).append(new_value)
					except ValueError:
						sys.exit(err)
						
			####  Read the Dice related keywords
			elif key in self.dice_keywords and len(value) != 0:  ##  'value' is not empty!
				
				if key == 'title':
					setattr(self.dice, key, value)
				
				elif key in ('ljname', 'outname', 'progname'):
					setattr(self.dice, key, value[0])
				
				elif key in ('ncores', 'isave'):
					err = "Error: expected a positive integer for keyword {} in file {}".format(key, self.infile)
					if not value[0].isdigit():
						sys.exit(err)
					new_value = int(value[0])
					if new_value >= 1:
						setattr(self.dice, key, new_value)
				
				elif key in ('temp', 'press', 'dens'):  # Cannot be less than 1e-10
					err = "Error: expected a positive float for keyword {} in file {}".format(key, self.infile)
					try:
						new_value = float(value[0])
						if new_value < 1e-10:
							sys.exit(err)
						else:
							setattr(self.dice, key, new_value)
					except ValueError:
						sys.exit(err)
				
				elif key == 'nmol':  # If defined, must be well defined (only positive integer values)
					err = "Error: expected 1 to 4 positive integers for keyword {} in file {}".format(key, self.infile)
					args = min(4, len(value))
					for i in range(args):
						if value[i].isdigit():
							new_value = int(value[i])
							if new_value < 1:
								sys.exit(err)
							else:
								setattr(self.dice, key, new_value)
						elif i == 0:
							sys.exit(err)
						else:
							break
				
				elif key == 'nstep':  # If defined, must be well defined (only positive integer values)
					err = "Error: expected 2 or 3 positive integers for keyword {} in file {}".format(key, self.infile)
					if len(value) < 2:
						sys.exit(err)
					args = min(3, len(value))
					for i in range(args):
						if value[i].isdigit():
							new_value = int(value[i])
							if new_value < 1:
								sys.exit(err)
							else:
								setattr(self.dice, key, new_value)
						elif i < 2:
							sys.exit(err)
						else:
							break
			
			####  Read the Gaussian related keywords
			elif key in self.gaussian_keywords and len(value) != 0:  ##  'value' is not empty!
				
				if key == 'mem':  # Memory in MB (minimum of 100)
					err = "Error: expected a positive integer for keyword {} in file {}".format(key, self.infile)
					if not value[0].isdigit():
						sys.exit(err)
					new_value = int(value[0])
					if new_value >= 100:
						setattr(self.gaussian, key, new_value)
				
				elif key == 'keywords':
					setattr(self.gaussian, key, value)
				
				elif key == 'chgmult':  # If defined, must be well defined (2 integer values)
					err = "Error: expected 2 integers for keyword {} in file {}".format(key, self.infile)
					if len(value) < 2:
						sys.exit(err)
					for i in range (2):
						try:
							setattr(self.gaussian, key)[i] = int(value[i])
						except ValueError:
							sys.exit(err)
				
				elif key in ('level', 'chglevel'):
					setattr(self.gaussian, key, value[0])
				
				elif key in ('gmiddle', 'gbottom'):
					setattr(self.gaussian, key, value[0])
				
				elif key == 'pop' and value[0].lower() in ("chelpg", "mk", "nbo"):
					setattr(self.gaussian, key, value[0].lower())
			
			####  Read the Molcas related keywords
			elif key in self.molcas_keywords and len(value) != 0:  ##  'value' is not empty!
				
				if key == 'root': # If defined, must be well defined (only positive integer values)
					err = "Error: expected a positive integer for keyword {} in file {}".format(key, self.infile)
					if not value[0].isdigit():
						sys.exit(err)
					new_value = int(value[0])
					if new_value >= 1:
						setattr(self.molcas, key, new_value)
				
				elif key in ('mbottom', 'orbfile'):
					setattr(self.molcas, key, value[0])
				
				elif key == 'basis':
					setattr(self.molcas ,key, value[0])
			
			#### End

		
	def check_keywords(self):
		
		min_steps = 20000
		
		if self.dice.ljname == None:
			sys.exit("Error: 'ljname' keyword not specified in file {}".format(self.infile))
		
		if self.dice.outname == None:
			sys.exit("Error: 'outname' keyword not specified in file {}".format(self.infile))
		
		if self.dice.dens == None:
			sys.exit("Error: 'dens' keyword not specified in file {}".format(self.infile))
		
		if len(self.dice.nmol) == 0:
			sys.exit("Error: 'nmol' keyword not defined appropriately in file {}".format(self.infile))
		
		if len(self.dice.nstep) == 0:
			sys.exit("Error: 'nstep' keyword not defined appropriately in file {}".format(self.infile))
		
		## Check only if QM program is Gaussian:
		if self.player.qmprog in ("g03", "g09", "g16"):

			if self.gaussian.level == None:
				sys.exit("Error: 'level' keyword not specified in file {}".format(self.infile))
		
			if self.gaussian.gmiddle != None:
				if not os.path.isfile(self.gaussian.gmiddle):
					sys.exit("Error: file {} not found".format(self.gaussian.gmiddle))

			if self.gaussian.gbottom != None:
				if not os.path.isfile(self.gaussian.gbottom):
					sys.exit("Error: file {} not found".format(self.gaussian.gbottom))
			
			if self.gaussian.pop != "chelpg" and (self.player.ghosts == "yes" or self.player.lps == "yes"):
				sys.exit("Error: ghost atoms or lone pairs only available with 'pop = chelpg')")
			
			if self.gaussian.chglevel == None:
				self.gaussian.chglevel = self.gaussian.level
		
		## Check only if QM program is Molcas:
		# if self.player.qmprog == "molcas":
			
		# 	if self.molcas.mbottom == None:
		# 		sys.exit("Error: 'mbottom' keyword not specified in file {}".format(self.infile))
		# 	else:
		# 		if not os.path.isfile(self.molcas.mbottom):
		# 			sys.exit("Error: file {} not found".format(self.molcas.mbottom))
			
		# 	if self.molcas.basis == None:
		# 		sys.exit("Error: 'basis' keyword not specified in file {}".format(self.infile))
		
		
		if self.player.altsteps != 0:
			
			### Verifica se tem mais de 1 molecula QM
			### (No futuro usar o RMSD fit para poder substituir todas as moleculas QM
			### no arquivo outname.xy - Need to change the make_init_file!!)
			if self.dice.nmol[0] > 1:
				sys.exit("Error: altsteps > 0 only possible with 1 QM molecule (nmol = 1 n2 n3 n4)")
						
			# if not zero, altsteps cannot be less than min_steps
			self.player.altsteps = max(min_steps, self.player.altsteps)
			# altsteps value is always the nearest multiple of 1000
			self.player.altsteps = round(self.player.altsteps / 1000) * 1000
		
		
		for i in range(len(self.dice.nstep)):
			# nstep can never be less than min_steps
			self.dice.nstep[i] = max(min_steps, self.dice.nstep[i])
			# nstep values are always the nearest multiple of 1000
			self.dice.nstep[i] = round(self.dice.nstep[i] / 1000) * 1000
		
		# isave must be between 100 and 2000
		self.dice.isave = max(100, self.dice.isave)
		self.dice.isave = min(2000, self.dice.isave)
		# isave value is always the nearest multiple of 100
		self.dice.isave = round(self.dice.isave / 100) * 100

	def print_keywords(self, fh):
		
		fh.write("##########################################################################################\n"
				"#############               Welcome to DICEPLAYER version 1.0                #############\n"
				"##########################################################################################\n"
				"\n")
		fh.write("Your python version is {}\n".format(sys.version))
		fh.write("\n")
		fh.write("Program started on {}\n".format(weekday_date_time()))
		fh.write("\n")
		fh.write("Environment variables:\n")
		for var in env:
			fh.write("{} = {}\n".format(var, 
									(os.environ[var] if var in os.environ else "Not set")))
		
		fh.write("\n==========================================================================================\n"
				"                         CONTROL variables being used in this run:\n"
				"------------------------------------------------------------------------------------------\n"
				"\n")

		for key in sorted(self.player_keywords):
			if getattr(self.player,key) != None:
				if isinstance(getattr(self.player,key), list):
					string = " ".join(str(x) for x in getattr(self.player,key))
					fh.write("{} = {}\n".format(key, string))
				else:	
					fh.write("{} = {}\n".format(key, getattr(self.player,key)))
		
		fh.write("\n")

		fh.write("------------------------------------------------------------------------------------------\n"
				"                         DICE variables being used in this run:\n"
				"------------------------------------------------------------------------------------------\n"
				"\n")

		for key in sorted(self.dice_keywords):
			if getattr(self.dice,key) != None:
				if isinstance(getattr(self.dice,key), list):
					string = " ".join(str(x) for x in getattr(self.dice,key))
					fh.write("{} = {}\n".format(key, string))
				else:	
					fh.write("{} = {}\n".format(key, getattr(self.dice,key)))
		
		fh.write("\n")
		
		if self.player.qmprog in ("g03", "g09", "g16"):

			fh.write("------------------------------------------------------------------------------------------\n"
					"                         GAUSSIAN variables being used in this run:\n"
					"------------------------------------------------------------------------------------------\n"
					"\n")
			
			for key in sorted(self.gaussian_keywords):
				if getattr(self.gaussian,key) != None:
					if isinstance(getattr(self.gaussian,key), list):
						string = " ".join(str(x) for x in getattr(self.gaussian,key))
						fh.write("{} = {}\n".format(key, string))
					else:	
						fh.write("{} = {}\n".format(key, getattr(self.gaussian,key)))
			
			fh.write("\n")
		
		# elif self.player.qmprog == "molcas":

		# 	fh.write("------------------------------------------------------------------------------------------\n"
		# 			"                         MOLCAS variables being used in this run:\n"
		# 			"------------------------------------------------------------------------------------------\n"
		# 			"\n")
			
		# 	for key in sorted(molcas):
		# 		if molcas[key] != None:
		# 			if isinstance(molcas[key], list):
		# 				string = " ".join(str(x) for x in molcas[key])
		# 				fh.write("{} = {}\n".format(key, string))
		# 			else:	
		# 				fh.write("{} = {}\n".format(key, molcas[key]))
			
		# 	fh.write("\n")

	def read_potential(self): # Deve ser atualizado para o uso de 
		
		try:
			with open(self.dice.ljname) as file:
				ljfile = file.readlines()
		except EnvironmentError as err:
			sys.exit(err)
		
		combrule = ljfile.pop(0).split()[0]
		if combrule not in ("*", "+"):
			sys.exit("Error: expected a '*' or a '+' sign in 1st line of file {}".format(self.dice.ljname))
		self.dice.combrule = combrule
		
		ntypes = ljfile.pop(0).split()[0]
		if not ntypes.isdigit():
			sys.exit("Error: expected an integer in the 2nd line of file {}".format(self.dice.ljname))
		ntypes = int(ntypes)
		
		if ntypes != len(self.dice.nmol):
			sys.exit("Error: number of molecule types in file {} must match that of 'nmol' keyword in file {}".format(
																	self.dice.ljname, self.infile))
		
		line = 2
		for i in range(ntypes):

			line += 1
			nsites = ljfile.pop(0).split()[0]
			if not nsites.isdigit():
				sys.exit("Error: expected an integer in line {} of file {}".format(line, self.dice.ljname))
			
			nsites = int(nsites)
			self.system.add_type(Molecule())
			
			for j in range(nsites):

				line += 1
				new_atom = ljfile.pop(0).split()
				
				if len(new_atom) < 8:
					sys.exit("Error: expected at least 8 fields in line {} of file {}".format(line, dice['ljname']))
				
				self.system.molecule[i].add_atom()
				
				if not new_atom[0].isdigit():
					sys.exit("Error: expected an integer in field 1, line {} of file {}".format(line, dice['ljname']))
				lbl = int(new_atom[0])
				
				if not new_atom[1].isdigit():
					sys.exit("Error: expected an integer in field 2, line {} of file {}".format(line, dice['ljname']))
				
				atnumber = int(new_atom[1])
				if atnumber == ghost_number and i == 0:  # Ghost atom not allowed in the QM molecule
					sys.exit("Error: found a ghost atom in line {} of file {}".format(line, dice['ljname']))
				na = atnumber
				
				try:
					rx = float(new_atom[2])
				except:
					sys.exit("Error: expected a float in field 3, line {} of file {}".format(line, dice['ljname']))
				
				try:
					ry = float(new_atom[3])
				except:
					sys.exit("Error: expected a float in field 4, line {} of file {}".format(line, dice['ljname']))
				
				try:
					rz = float(new_atom[4])
				except:
					sys.exit("Error: expected a float in field 5, line {} of file {}".format(line, dice['ljname']))
				
				try:
					chg = float(new_atom[5])
				except:
					sys.exit("Error: expected a float in field 6, line {} of file {}".format(line, dice['ljname']))
				
				try:
					eps = float(new_atom[6])
				except:
					sys.exit("Error: expected a float in field 7, line {} of file {}".format(line, dice['ljname']))
				
				try:
					sig = float(new_atom[7])
				except:
					sys.exit("Error: expected a float in field 8, line {} of file {}".format(line, dice['ljname']))
				
				mass = atommass[na]
				
				if len(new_atom) > 8:
					masskey, mass = new_atom[8].partition("=")[::2]
					if masskey.lower() == 'mass' and len(mass) !=0:
						try:
							new_mass = float(mass)
							if new_mass > 0:
								mass = new_mass
						except:
							sys.exit(
							"Error: expected a positive float after 'mass=' in field 9, line {} of file {}".format(
																							line, dice['ljname']))

				self.system.molecule[i].add_atom(Atom(lbl,na,rx,ry,rz,chg,eps,sig,mass))
				
		to_delete = ['lbl','na','rx','ry','rz','chg','eps','sig','mass']
		for _var in to_delete:
			if _var in locals() or _var in globals():
				exec(f'del {_var}')

	class Player:

		def __init__(self):

			self.maxcyc = None
			# self.initcyc = 1 Eliminated
			self.nprocs = 1
			self.switchcyc = 3
			self.altsteps = 20000
			self.maxstep = .3
			self.qmprog = "g09"
			self.opt = "yes"
			self.freq = "no"
			self.readhessian = "no"
			self.lps = "no"
			self.ghosts = "no"
			self.vdwforces = "no"
			self.tol_factor = 1.2
		


	class Dice:

		def __init__(self):

			self.title = "Diceplayer run"
			self.progname = "dice"
			self.temp = 300.0
			self.press = 1.0
			self.isave = 1000         # ASEC construction will take this into account
			self.ncores = 1

			self.dens = None		# Investigate the possibility of using 'box = Lx Ly Lz' instead.
			# self.box = None		# So 'geom' would be set by diceplayer and 'cutoff' would be
									# switched off. One of them must be given.
			self.ljname = None
			self.outname = None
			self.nmol = []     # Up to 4 integer values related to up to 4 molecule types
			self.nstep = []    # 2 or 3 integer values related to 2 or 3 simulations
								# (NVT th + NVT eq) or (NVT th + NPT th + NPT eq).
								# This will control the 'nstep' keyword of Dice


	class Gaussian:

		def __init__(self):

			self.mem = None
			self.keywords = None
			self.chgmult = [0, 1]
			self.gmiddle = None   # In each case, if a filename is given, its content will be 
			self.gbottom = None   # inserted in the gaussian input
			self.pop = "chelpg"
			self.chglevel = None

			self.level = None

	# class Molcas:

	# 	def __init(self):

	# 		self.orbfile = "input.exporb"
	# 		self.root = 1

	# 		self.mbottom = None
	# 		self.basis = None