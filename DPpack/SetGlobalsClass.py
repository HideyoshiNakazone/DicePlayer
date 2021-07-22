import os, sys
import shutil
import textwrap

from DPpack.PTable import *
from DPpack.Misc import *

class Internal:

	def __init__(self, infile):

		self.infile = infile

		self.player = self.Player()
		self.player_keywords = [a for a in dir(self.player) if not a.startswith('__') and not callable(getattr(self.player, a))]

		self.dice = self.Dice()
		self.dice_keywords = [a for a in dir(self.dice) if not a.startswith('__') and not callable(getattr(self.dice, a))]

		self.gaussian = self.Gaussian()
		self.gaussian_keywords = [a for a in dir(self.gaussian) if not a.startswith('__') and not callable(getattr(self.gaussian, a))]

		self.molcas = self.Molcas()
		self.molcas_keywords = [a for a in dir(self.molcas) if not a.startswith('__') and not callable(getattr(self.molcas, a))]

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
					err = "Error: expected a positive integer for keyword {} in file {}".format(key, infile)
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

	class Player:

		def __init__(self):

			self.maxcyc = None
			self.initcyc = 1
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
			#dice['box'] = None		# So 'geom' would be set by diceplayer and 'cutoff' would be
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

	class Molcas:

		def __init(self):

			self.orbfile = "input.exporb"
			self.root = 1

			self.mbottom = None
			self.basis = None