import os, sys, time
import shutil, gzip

#######################################  functions  ######################################

def weekday_date_time():
	
	return time.strftime("%A, %d %b %Y at %H:%M:%S")


def date_time():
	
	return time.strftime("%d %b %Y at %H:%M:%S")


def compress_files_1mb(path):
	
	working_dir = os.getcwd()
	os.chdir(path)
	
	files = filter(os.path.isfile, os.listdir(os.curdir))
	for file in files:
		if os.path.getsize(file) > 1024 * 1024:  ## If bigger than 1MB
			filegz = file + ".gz"
			try:
				with open(file, 'rb') as f_in:
					with gzip.open(filegz, 'wb') as f_out:
						shutil.copyfileobj(f_in, f_out)
			except:
				sys.exit("Error: cannot compress file {}".format(file))
	
	os.chdir(working_dir)
	
	return


def make_step_dir(cycle):
	
	step_dir = "step{:02d}".format(cycle)
	if os.path.exists(step_dir):
		sys.exit("Error: a file or directory {} already exists".format(step_dir))
	try:
		os.makedirs(step_dir)
	except:
		sys.exit("Error: cannot make directory {}".format(step_dir))
	
	return



def make_qm_dir(cycle):
	
	step_dir = "step{:02d}".format(cycle)
	path = step_dir + os.sep + "qm"
	try:
		os.makedirs(path)
	except:
		sys.exit("Error: cannot make directory {}".format(path))
	
	return



