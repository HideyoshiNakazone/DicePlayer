import os, sys
import shutil
import textwrap

from DPpack.PTable import *
from DPpack.Misc import *

####  Global hashes that control the behaviour of Diceplayer

player = {}
dice = {}
gaussian = {}
molcas = {}
internal = {}

#######################################################################
####        Global parameters that MAY be given by the user        ####
####    (If not given by the user, default values will be used)    ####
#######################################################################

## Diceplayer:
player['maxcyc'] = 1
player['initcyc'] = 1        # May restart an optimization (append to geoms.xyz from start)
player['nprocs'] = 1
player['switchcyc'] = 3		 # At which step start doing only one QM calculation (geom & chg)
player['altsteps'] = 20000   # Steps for thermalization when starting from previous step
player['maxstep'] = 0.3      # Maxstep for geometry relaxation in Bohr
player['qmprog'] = "g09"
player['opt'] = "yes"
player['freq'] = "no"
player['readhessian'] = "no"
player['lps'] = "no"
player['ghosts'] = "no"
player['vdwforces'] = "no"
player['tol_factor'] = 1.2		# Factor to multiply the default tolerance values

## Dice:
dice['title'] = "Diceplayer run"
dice['progname'] = "dice"
dice['temp'] = 300.0
dice['press'] = 1.0
dice['isave'] = 1000         # ASEC construction will take this into account
dice['ncores'] = 1

## Gaussian:
gaussian['mem'] = None
gaussian['keywords'] = None
gaussian['chgmult'] = [0, 1]
gaussian['gmiddle'] = None   # In each case, if a filename is given, its content will be 
gaussian['gbottom'] = None   # inserted in the gaussian input
gaussian['pop'] = "chelpg"
gaussian['chglevel'] = None

## Molcas:
molcas['orbfile'] = "input.exporb"
molcas['root'] = 1


########################################################################
####        Global parameters that MUST be given by the user        ####
########################################################################

## Dice:
dice['dens'] = None		# Investigate the possibility of using 'box = Lx Ly Lz' instead.
#dice['box'] = None		# So 'geom' would be set by diceplayer and 'cutoff' would be
						# switched off. One of them must be given.
dice['ljname'] = None
dice['outname'] = None
dice['nmol'] = []     # Up to 4 integer values related to up to 4 molecule types
dice['nstep'] = []    # 2 or 3 integer values related to 2 or 3 simulations
                      # (NVT th + NVT eq) or (NVT th + NPT th + NPT eq).
                      # This will control the 'nstep' keyword of Dice

## Gaussian:
gaussian['level'] = None

## Molcas:
molcas['mbottom'] = None
molcas['basis'] = None


## The following Dice keywords will be handled automatically by Diceplayer:
## * init ("yes" in the first thermalization and "yesreadxyz" for thermalizations  
##         starting from a previous step / "no" in subsequent simulations)
## * vstep ("0" for NVT simulations and 'nstep'/5 for NPT simulations)
## * nstep ('nstep' for NVT and "5" for NPT simulations )
## * irdf ("0" for thermalizations and "10*nprocs" for equilibrium)
## * seed (will be generated randomly each time a Dice input is created)

## The following Dice keywords will be set constant by Diceplayer for all simulations
## * mstop = 1 (So to guarantee that the ASEC will be correctly built)
## * accum = no (There is never a simulation continuation in Diceplayer)
## * iprint = 1 (Print energy info every step in Dice output)

## All the other Dice keywords will not be altered from their default values
## and therefore are not mentioned in Diceplayer


#####################################################################
####    Global parameters that are not accessible to the user    ####
####    (Intended to be used only internally by the program)     ####
#####################################################################

## Diceplayer:
internal['tol_rms_force'] = 3e-4		# Hartree/Bohr
internal['tol_max_force'] = 4.5e-4		# Hartree/Bohr
internal['tol_rms_step'] = 1.2e-3		# Bohr
internal['tol_max_step'] = 1.8e-3		# Bohr
internal['trust_radius'] = None

## Dice:
internal['combrule'] = None
internal['randominit'] = None

## Other global variables:
molecules = []  # Armazena todas as informacoes sobre cada tipo de molecula 
                # (lbl, na, rx, ry, rz, chg, eps, sig, mass)

internal['ghost_types'] = []
internal['ghost_atoms'] = []  # Store the ghost atoms (off-atom charge sites) in the QM molecule
        			          # (rx, ry, rz, chg)
internal['lp_types'] = []
internal['lp_atoms'] = []  # Store the lone pairs (special off-atom charge sites) in the QM molecule
                		   # (rx, ry, rz, chg)

## Numpy arrays:
step = []						## Values in Bohr
internal['position'] = []
internal['energy'] = []			## Values in Hartree
internal['gradient'] = []		## Values in Hartree/Bohr
internal['hessian'] = []		## Values in Hartree/Bohr^2

## Conversion factors:
bohr2ang = 0.52917721092
ang2bohr = 1/bohr2ang

######################################################################
####   Environment variables important for executing Diceplayer   ####
######################################################################

env = ["OMP_STACKSIZE"]

#######################################  functions  ######################################
##  Functions to process the input files and store the values in the global variables   ##
##########################################################################################

def read_keywords(infile):

	try:
		with open(infile) as fh:
			controlfile = fh.readlines()
	except EnvironmentError:
		sys.exit("Error: cannot open file {}".format(infile))
	
	for line in controlfile:
		
		key, value = line.partition("=")[::2]  # Discards the '='
		key = key.strip().lower()
		if key in ('title', 'keywords'):
			value = value.strip()
		else:
			value = value.split()
		
		####  Read the Diceplayer related keywords
		if key in player and len(value) != 0:  ##  'value' is not empty!
			
			if key == 'qmprog' and value[0].lower() in ("g03", "g09", "g16", "molcas"):
				player[key] = value[0].lower()
			
			elif key == 'opt' and value[0].lower() in ("yes", "no", "ts"):
				player[key] = value[0].lower()
			
			#elif key == 'zipprog' and value[0].lower() in ("zip", "gzip", "bzip"):
				#player[key] = value[0].lower()
			
			elif key in ('lps', 'ghosts') and value[0].lower() in ("yes", "no"):
				player[key] = value[0].lower()
			
			elif key in ('readhessian', 'vdwforces') and value[0].lower() in ("yes", "no"):
				player[key] = value[0].lower()
			
			elif key in ('maxcyc', 'initcyc', 'nprocs', 'altsteps', 'switchcyc'):
				err = "Error: expected a positive integer for keyword {} in file {}".format(key, infile)
				try:
					new_value = int(value[0])
					if new_value >= 1:
						player[key] = new_value
					elif key == 'altsteps' and new_value == 0:
						player[key] = 0
				except ValueError:
					sys.exit(err)
			
			elif key == 'maxstep':  # Cannot be less than 0.01
				err = "Error: expected a float greater than 0.01 for keyword {} in file {}".format(key, infile)
				try:
					new_value = float(value[0])
					if new_value < 0.01:
						sys.exit(err)
					else:
						player[key] = new_value
				except ValueError:
					sys.exit(err)
					
		####  Read the Dice related keywords
		elif key in dice and len(value) != 0:  ##  'value' is not empty!
			
			if key == 'title':
				dice[key] = value
			
			elif key in ('ljname', 'outname', 'progname'):
				dice[key] = value[0]
			
			elif key in ('ncores', 'isave'):
				err = "Error: expected a positive integer for keyword {} in file {}".format(key, infile)
				if not value[0].isdigit():
					sys.exit(err)
				new_value = int(value[0])
				if new_value >= 1:
					dice[key] = new_value
			
			elif key in ('temp', 'press', 'dens'):  # Cannot be less than 1e-10
				err = "Error: expected a positive float for keyword {} in file {}".format(key, infile)
				try:
					new_value = float(value[0])
					if new_value < 1e-10:
						sys.exit(err)
					else:
						dice[key] = new_value
				except ValueError:
					sys.exit(err)
			
			elif key == 'nmol':  # If defined, must be well defined (only positive integer values)
				err = "Error: expected 1 to 4 positive integers for keyword {} in file {}".format(key, infile)
				args = min(4, len(value))
				for i in range(args):
					if value[i].isdigit():
						new_value = int(value[i])
						if new_value < 1:
							sys.exit(err)
						else:
							dice[key].append(new_value)
					elif i == 0:
						sys.exit(err)
					else:
						break
			
			elif key == 'nstep':  # If defined, must be well defined (only positive integer values)
				err = "Error: expected 2 or 3 positive integers for keyword {} in file {}".format(key, infile)
				if len(value) < 2:
					sys.exit(err)
				args = min(3, len(value))
				for i in range(args):
					if value[i].isdigit():
						new_value = int(value[i])
						if new_value < 1:
							sys.exit(err)
						else:
							dice[key].append(new_value)
					elif i < 2:
						sys.exit(err)
					else:
						break
		
		####  Read the Gaussian related keywords
		elif key in gaussian and len(value) != 0:  ##  'value' is not empty!
			
			if key == 'mem':  # Memory in MB (minimum of 100)
				err = "Error: expected a positive integer for keyword {} in file {}".format(key, infile)
				if not value[0].isdigit():
					sys.exit(err)
				new_value = int(value[0])
				if new_value >= 100:
					gaussian[key] = new_value
			
			elif key == 'keywords':
				gaussian[key] = value
			
			elif key == 'chgmult':  # If defined, must be well defined (2 integer values)
				err = "Error: expected 2 integers for keyword {} in file {}".format(key, infile)
				if len(value) < 2:
					sys.exit(err)
				for i in range (2):
					try:
						gaussian[key][i] = int(value[i])
					except ValueError:
						sys.exit(err)
			
			elif key in ('level', 'chglevel'):
				gaussian[key] = value[0]
			
			elif key in ('gmiddle', 'gbottom'):
				gaussian[key] = value[0]
			
			elif key == 'pop' and value[0].lower() in ("chelpg", "mk", "nbo"):
				gaussian[key] = value[0].lower()
		
		####  Read the Molcas related keywords
		elif key in molcas and len(value) != 0:  ##  'value' is not empty!
			
			if key == 'root': # If defined, must be well defined (only positive integer values)
				err = "Error: expected a positive integer for keyword {} in file {}".format(key, infile)
				if not value[0].isdigit():
					sys.exit(err)
				new_value = int(value[0])
				if new_value >= 1:
					molcas[key] = new_value
			
			elif key in ('mbottom', 'orbfile'):
				molcas[key] = value[0]
			
			elif key == 'basis':
				molcas[key] = value[0]
		
		#### End
	return



def check_keywords(infile):
	
	min_steps = 20000
	
	if dice['ljname'] == None:
		sys.exit("Error: 'ljname' keyword not specified in file {}".format(infile))
	
	if dice['outname'] == None:
		sys.exit("Error: 'outname' keyword not specified in file {}".format(infile))
	
	if dice['dens'] == None:
		sys.exit("Error: 'dens' keyword not specified in file {}".format(infile))
	
	if len(dice['nmol']) == 0:
		sys.exit("Error: 'nmol' keyword not defined appropriately in file {}".format(infile))
	
	if len(dice['nstep']) == 0:
		sys.exit("Error: 'nstep' keyword not defined appropriately in file {}".format(infile))
	
	## Check only if QM program is Gaussian:
	if player['qmprog'] in ("g03", "g09", "g16"):
		if gaussian['level'] == None:
			sys.exit("Error: 'level' keyword not specified in file {}".format(infile))
	
		if gaussian['gmiddle'] != None:
			if not os.path.isfile(gaussian['gmiddle']):
				sys.exit("Error: file {} not found".format(gaussian['gmiddle']))

		if gaussian['gbottom'] != None:
			if not os.path.isfile(gaussian['gbottom']):
				sys.exit("Error: file {} not found".format(gaussian['gbottom']))
		
		if gaussian['pop'] != "chelpg" and (player['ghosts'] == "yes" or player['lps'] == "yes"):
			sys.exit("Error: ghost atoms or lone pairs only available with 'pop = chelpg')")
		
		if gaussian['chglevel'] == None:
			gaussian['chglevel'] = gaussian['level']
	
	## Check only if QM program is Molcas:
	if player['qmprog'] == "molcas":
		
		if molcas['mbottom'] == None:
			sys.exit("Error: 'mbottom' keyword not specified in file {}".format(infile))
		else:
			if not os.path.isfile(molcas['mbottom']):
				sys.exit("Error: file {} not found".format(molcas['mbottom']))
		
		if molcas['basis'] == None:
			sys.exit("Error: 'basis' keyword not specified in file {}".format(infile))
	
	
	if player['altsteps'] != 0:
		
		### Verifica se tem mais de 1 molecula QM
		### (No futuro usar o RMSD fit para poder substituir todas as moleculas QM
		### no arquivo outname.xy - Need to change the make_init_file!!)
		if dice['nmol'][0] > 1:
			sys.exit("Error: altsteps > 0 only possible with 1 QM molecule (nmol = 1 n2 n3 n4)")
					
		# if not zero, altsteps cannot be less than min_steps
		player['altsteps'] = max(min_steps, player['altsteps'])
		# altsteps value is always the nearest multiple of 1000
		player['altsteps'] = round(player['altsteps'] / 1000) * 1000
	
	
	for i in range(len(dice['nstep'])):
		# nstep can never be less than min_steps
		dice['nstep'][i] = max(min_steps, dice['nstep'][i])
		# nstep values are always the nearest multiple of 1000
		dice['nstep'][i] = round(dice['nstep'][i] / 1000) * 1000
	
	# isave must be between 100 and 2000
	dice['isave'] = max(100, dice['isave'])
	dice['isave'] = min(2000, dice['isave'])
	# isave value is always the nearest multiple of 100
	dice['isave'] = round(dice['isave'] / 100) * 100
		
	return



def print_keywords(fh):
	
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

	for key in sorted(player):
		if player[key] != None:
			if isinstance(player[key], list):
				string = " ".join(str(x) for x in player[key])
				fh.write("{} = {}\n".format(key, string))
			else:	
				fh.write("{} = {}\n".format(key, player[key]))
	
	fh.write("\n")

	fh.write("------------------------------------------------------------------------------------------\n"
	         "                         DICE variables being used in this run:\n"
	         "------------------------------------------------------------------------------------------\n"
	         "\n")

	for key in sorted(dice):
		if dice[key] != None:
			if isinstance(dice[key], list):
				string = " ".join(str(x) for x in dice[key])
				fh.write("{} = {}\n".format(key, string))
			else:	
				fh.write("{} = {}\n".format(key, dice[key]))
	
	fh.write("\n")
	
	if player['qmprog'] in ("g03", "g09", "g16"):

		fh.write("------------------------------------------------------------------------------------------\n"
	             "                         GAUSSIAN variables being used in this run:\n"
	             "------------------------------------------------------------------------------------------\n"
	             "\n")
		
		for key in sorted(gaussian):
			if gaussian[key] != None:
				if isinstance(gaussian[key], list):
					string = " ".join(str(x) for x in gaussian[key])
					fh.write("{} = {}\n".format(key, string))
				else:	
					fh.write("{} = {}\n".format(key, gaussian[key]))
		
		fh.write("\n")
	
	elif player['qmprog'] == "molcas":

		fh.write("------------------------------------------------------------------------------------------\n"
	             "                         MOLCAS variables being used in this run:\n"
	             "------------------------------------------------------------------------------------------\n"
	             "\n")
		
		for key in sorted(molcas):
			if molcas[key] != None:
				if isinstance(molcas[key], list):
					string = " ".join(str(x) for x in molcas[key])
					fh.write("{} = {}\n".format(key, string))
				else:	
					fh.write("{} = {}\n".format(key, molcas[key]))
		
		fh.write("\n")
		
	return



def read_potential(infile): # Deve ser atualizado para o uso de 
	
	try:
		with open(dice['ljname']) as file:
			ljfile = file.readlines()
	except EnvironmentError as err:
		sys.exit(err)
	
	combrule = ljfile.pop(0).split()[0]
	if combrule not in ("*", "+"):
		sys.exit("Error: expected a '*' or a '+' sign in 1st line of file {}".format(dice['ljname']))
	dice['combrule'] = combrule
	
	ntypes = ljfile.pop(0).split()[0]
	if not ntypes.isdigit():
		sys.exit("Error: expected an integer in the 2nd line of file {}".format(dice['ljname']))
	ntypes = int(ntypes)
	
	if ntypes != len(dice['nmol']):
		sys.exit("Error: number of molecule types in file {} must match that of 'nmol' keyword in file {}".format(
		                                                          dice['ljname'], infile))
	
	line = 2
	for i in range(ntypes):
		line += 1
		nsites = ljfile.pop(0).split()[0]
		if not nsites.isdigit():
			sys.exit("Error: expected an integer in line {} of file {}".format(line, dice['ljname']))
		
		nsites = int(nsites)
		molecules.append([])
		
		for j in range(nsites):
			line += 1
			new_atom = ljfile.pop(0).split()
			
			if len(new_atom) < 8:
				sys.exit("Error: expected at least 8 fields in line {} of file {}".format(line, dice['ljname']))
			
			molecules[i].append({})
			
			if not new_atom[0].isdigit():
				sys.exit("Error: expected an integer in field 1, line {} of file {}".format(line, dice['ljname']))
			molecules[i][j]['lbl'] = int(new_atom[0])
			
			if not new_atom[1].isdigit():
				sys.exit("Error: expected an integer in field 2, line {} of file {}".format(line, dice['ljname']))
			
			atnumber = int(new_atom[1])
			if atnumber == ghost_number and i == 0:  # Ghost atom not allowed in the QM molecule
				sys.exit("Error: found a ghost atom in line {} of file {}".format(line, dice['ljname']))
			molecules[i][j]['na'] = atnumber
			
			try:
				molecules[i][j]['rx'] = float(new_atom[2])
			except:
				sys.exit("Error: expected a float in field 3, line {} of file {}".format(line, dice['ljname']))
			
			try:
				molecules[i][j]['ry'] = float(new_atom[3])
			except:
				sys.exit("Error: expected a float in field 4, line {} of file {}".format(line, dice['ljname']))
			
			try:
				molecules[i][j]['rz'] = float(new_atom[4])
			except:
				sys.exit("Error: expected a float in field 5, line {} of file {}".format(line, dice['ljname']))
			
			try:
				molecules[i][j]['chg'] = float(new_atom[5])
			except:
				sys.exit("Error: expected a float in field 6, line {} of file {}".format(line, dice['ljname']))
			
			try:
				molecules[i][j]['eps'] = float(new_atom[6])
			except:
				sys.exit("Error: expected a float in field 7, line {} of file {}".format(line, dice['ljname']))
			
			try:
				molecules[i][j]['sig'] = float(new_atom[7])
			except:
				sys.exit("Error: expected a float in field 8, line {} of file {}".format(line, dice['ljname']))
			
			molecules[i][j]['mass'] = atommass[molecules[i][j]['na']]
			
			if len(new_atom) > 8:
				masskey, mass = new_atom[8].partition("=")[::2]
				if masskey.lower() == 'mass' and len(mass) !=0:
					try:
						new_mass = float(mass)
						if new_mass > 0:
							molecules[i][j]['mass'] = new_mass
					except:
						sys.exit(
						"Error: expected a positive float after 'mass=' in field 9, line {} of file {}".format(
						                                                                  line, dice['ljname']))
	
	return



def read_ghosts():
	
	max_atom_number = len(molecules[0])

	try:
		with open("ghosts.in") as fh:
			ghostfile = fh.readlines()
	except EnvironmentError:
		sys.exit("Error: cannot open file ghosts.in")
	
	for line in ghostfile:
		
		if len(line.split()) > 1:  # Discard lines with less than 2 fields
			
			key, *atom_numbers = line.split()
			key = key.lower()
			
			if key in ("g", "m", "z"):  # Discard lines that do not start with g|m|z
				ghost_types.append({})
				ghost_types[-1]['type'] = key
				ghost_types[-1]['numbers'] = []
				
				for num in atom_numbers:
					if not num.isdigit():
						sys.exit("Error: in file ghosts.in: only positive integers allowed after letter g|m|z")
					new_num = int(num)
					if new_num > max_atom_number:
						sys.exit("Error: in file ghosts.in: there is no atom number {}".format(new_num))
					else:
						ghost_types[-1]['numbers'].append(new_num)
				
				if len(ghost_types[-1]['numbers']) < 2:
					sys.exit("Error: in file ghosts.in: at least 2 atoms are necessary to make a ghost")
	
	if len(ghost_types) == 0:
		sys.exit("Error: no ghost atom found in ghosts.in")
	
	return



def read_lps():
	
	lp_alpha = 104.0  # Default values
	lp_dist = 0.7     #
	max_lp_type = 2
	min_alpha = 90.0
	max_alpha = 150.0
	min_dist = 0.5
	max_dist = 1.5
	max_atom_number = len(molecules[0])
	
	try:
		with open("lps.in") as fh:
			lpfile = fh.readlines()
	except EnvironmentError:
		sys.exit("Error: cannot open file lps.in")
	
	for line in lpfile:
		
		if len(line.split()) > 1:  # Discard lines with less than 2 fields
			
			type, *atom_numbers = line.split()
		
			if type.isdigit():  # Discard lines that do not start with an integer
				new_type = int(type)
				if new_type > max_lp_type:
					sys.exit("Error: in file lps.in: allowed LP types from 1 to {}".format(max_lp_type))
				lp_types.append({})
				lp_types[-1]['type'] = new_type
				lp_types[-1]['numbers'] = []
				
				# Read types 1 and 2
				if new_type in (1, 2):
				
					if len(atom_numbers) < 3:
						sys.exit("Error: in file lps.in: at least 3 atoms are necessary to make LPs type 1 and 2")
					for i in range(3):
						num = atom_numbers.pop(0)
						if not num.isdigit():
							sys.exit("Error: in file lps.in: expected 3 atom numbers after LPs type 1 and 2")
						new_num = int(num)
						if new_num > max_atom_number or new_num < 1:
							sys.exit("Error: in file lps.in: there is no atom number {}".format(new_num))
						else:
							lp_types[-1]['numbers'].append(new_num)
					
					lp_types[-1]['alpha'] = lp_alpha
					lp_types[-1]['dist'] = lp_dist
					
					if len(atom_numbers) != 0:
						try:
							alpha = float(atom_numbers.pop(0))
							if alpha > min_alpha and alpha < max_alpha:
								lp_types[-1]['alpha'] = alpha
							else:
								atom_numbers = []
						except:
							atom_numbers = []
					
					if len(atom_numbers) != 0:
						try:
							dist = float(atom_numbers.pop(0))
							if dist > min_dist and dist < max_dist:
								lp_types[-1]['dist'] = dist
						except:
							None
				# End of types 1 and 2
					
	if len(lp_types) == 0:
		sys.exit("Error: no lone pair found in lps.in")

	return



def print_potential(fh):
	
	formatstr = "{:<3d} {:>3d}  {:>10.5f} {:>10.5f} {:>10.5f}  {:>10.6f} {:>9.5f} {:>7.4f} {:>9.4f}\n"
	fh.write("\n"
	         "==========================================================================================\n")
	fh.write("                    Potential parameters from file {}:\n".format(dice['ljname']))
	fh.write("------------------------------------------------------------------------------------------\n"
	         "\n")
	
	fh.write("Combination rule: {}\n".format(dice['combrule']))
	fh.write("Types of molecules: {}\n\n".format(len(molecules)))
	
	i = 0
	for mol in molecules:
		i += 1
		fh.write("{} atoms in molecule type {}:\n".format(len(mol), i))
		fh.write("---------------------------------------------------------------------------------\n"
		         "Lbl  AN       X          Y          Z         Charge    Epsilon   Sigma     Mass\n")
		fh.write("---------------------------------------------------------------------------------\n")

		for atom in mol:
			fh.write(formatstr.format(atom['lbl'], atom['na'], atom['rx'], atom['ry'], atom['rz'],
			                          atom['chg'], atom['eps'], atom['sig'], atom['mass']))
		
		fh.write("\n")
		
	if player['ghosts'] == "yes" or player['lps'] == "yes":
		fh.write("\n"
	         "------------------------------------------------------------------------------------------\n"
	         "                    Aditional potential parameters:\n"
	         "------------------------------------------------------------------------------------------\n")
	
	if player['ghosts'] == "yes":
		
		fh.write("\n")
		fh.write("{} ghost atoms appended to molecule type 1 at:\n".format(len(ghost_types)))
		fh.write("---------------------------------------------------------------------------------\n")
		
		atoms_string = ""
		for ghost in ghost_types:
			for atom in ghost['numbers']:
				atom_sym = atomsymb[ molecules[0][atom - 1]['na'] ].strip()
				atoms_string += "{}{} ".format(atom_sym,atom)
			
			if ghost['type'] == "g":
				fh.write(textwrap.fill("* Geometric center of atoms {}".format(atoms_string), 80))
			elif ghost['type'] == "m":
				fh.write(textwrap.fill("* Center of mass of atoms {}".format(atoms_string), 80))
			elif ghost['type'] == "z":
				fh.write(textwrap.fill("* Center of atomic number of atoms {}".format(atoms_string), 80))
			
			fh.write("\n")
	
	if player['lps'] == 'yes':
		
		fh.write("\n")
		fh.write("{} lone pairs appended to molecule type 1:\n".format(len(lp_types)))
		fh.write("---------------------------------------------------------------------------------\n")
		
		for lp in lp_types:
			# LP type 1 or 2
			if lp['type'] in (1, 2):
				atom1_num = lp['numbers'][0]
				atom1_sym = atomsymb[ molecules[0][atom1_num - 1]['na'] ].strip()
				atom2_num = lp['numbers'][1]
				atom2_sym = atomsymb[ molecules[0][atom2_num - 1]['na'] ].strip()
				atom3_num = lp['numbers'][2]
				atom3_sym = atomsymb[ molecules[0][atom3_num - 1]['na'] ].strip()
				
				fh.write(textwrap.fill(
				"* Type {} on atom {}{} with {}{} {}{}. Alpha = {:<5.1f} Deg and D = {:<4.2f} Angs".format(
					lp['type'], atom1_sym, atom1_num, atom2_sym, atom2_num, atom3_sym, atom3_num, lp['alpha'],
					lp['dist']), 86))
				fh.write("\n")
				
			# Other LP types
			
	fh.write("\n"
	         "==========================================================================================\n")
	
	return

##  Creation of continue_function

def check_executables(fh):
	
	fh.write("\n")
	fh.write(90 * "=")
	fh.write("\n\n")
	
	dice_path = shutil.which(dice['progname'])
	if dice_path != None:
		fh.write("Program {} found at {}\n".format(dice['progname'], dice_path))
	else:
		sys.exit("Error: cannot find dice executable")

	qmprog_path = shutil.which(player['qmprog'])
	if qmprog_path != None:
		fh.write("Program {} found at {}\n".format(player['qmprog'], qmprog_path))
	else:
		sys.exit("Error: cannot find {} executable".format(player['qmprog']))
	
	if player['qmprog'] in ("g03", "g09", "g16"):
		formchk_path = shutil.which("formchk")
		if formchk_path != None:
			fh.write("Program formchk found at {}\n".format(formchk_path))
		else:
			sys.exit("Error: cannot find formchk executable")

	
	return



