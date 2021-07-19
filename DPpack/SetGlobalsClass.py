import os, sys
import shutil
import textwrap

from DPpack.PTable import *
from DPpack.Misc import *

class Internal:

	def __init__(self):

		self.player = self.Player()
		self.dice = self.Dice()
		self.gaussian = self.Gaussian()
		self.molca = self.Molca()

		## Constanst that shall be set for global use

		self.tol_rms_force = 3e-4		# Hartree/Bohr
		self.tol_max_force = 4.5e-4		# Hartree/Bohr
		self.tol_rms_step = 1.2e-3		# Bohr
		self.tol_max_step = 1.8e-3		# Bohr
		self.trust_radius = None

		## Dice:
		self.combrule = None
		self.randominit = None

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