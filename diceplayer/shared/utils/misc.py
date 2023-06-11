import gzip
import os
import shutil
import sys
import time
from typing import Final

#######################################  constants  ######################################


BOHR2ANG: Final[float] = 0.52917721092
ANG2BOHR: Final[float] = 1 / BOHR2ANG


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
                with open(file, "rb") as f_in:
                    with gzip.open(filegz, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
            except:
                sys.exit("Error: cannot compress file {}".format(file))

    os.chdir(working_dir)

    return


def make_step_dir(cycle):
    sim_dir = "simfiles"
    step_dir = "step{:02d}".format(cycle)
    path = sim_dir + os.sep + step_dir
    if os.path.exists(path):
        sys.exit("Error: a file or directory {} already exists".format(step_dir))
    try:
        os.makedirs(path)
    except:
        sys.exit("Error: cannot make directory {}".format(step_dir))


def make_qm_dir(cycle):
    sim_dir = "simfiles"
    step_dir = "step{:02d}".format(cycle)
    path = sim_dir + os.sep + step_dir + os.sep + "qm"
    try:
        os.makedirs(path)
    except:
        sys.exit("Error: cannot make directory {}".format(path))
