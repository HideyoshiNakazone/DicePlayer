#!/usr/bin/python3

import argparse
import os
import pickle
import shutil
import signal
import sys
import time
from multiprocessing import Process, connection

import numpy as np
import setproctitle

from diceplayer.DPpack.Environment.Atom import Atom
from diceplayer.DPpack.Environment.Molecule import Molecule
from diceplayer.DPpack.Player import Player
from diceplayer.DPpack.Utils.Misc import *

__version = "dev"
setproctitle.setproctitle("diceplayer-{}".format(__version))

if __name__ == "__main__":
    ####  Read and store the arguments passed to the program  ####
    ####  and set the usage and help messages                 ####

    parser = argparse.ArgumentParser(prog="Diceplayer")
    parser.add_argument(
        "--continue", dest="opt_continue", default=False, action="store_true"
    )
    parser.add_argument(
        "--version", action="version", version="diceplayer-" + __version
    )
    parser.add_argument(
        "-i",
        dest="infile",
        default="control.in",
        metavar="INFILE",
        help="input file of diceplayer [default = control.in]",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        default="run.log",
        metavar="OUTFILE",
        help="output file of diceplayer [default = run.log]",
    )
    ## Study the option of a parameter for continuing the last process via data from control.in and run.log files

    args = parser.parse_args()

    ####  Open OUTFILE for writing and print keywords and initial info

    try:

        if args.opt_continue and os.path.exists(args.outfile):

            save = pickle.load(open("latest-step.pkl", "rb"))

            if os.path.isfile(args.outfile + ".backup"):
                os.remove(args.outfile + ".backup")

            os.rename(args.outfile, args.outfile + ".backup")
            outfile = open(args.outfile, "w", 1)

        elif os.path.exists(args.outfile):
            os.rename(args.outfile, args.outfile + ".backup")
            outfile = open(args.outfile, "w", 1)
        else:
            outfile = open(args.outfile, "w", 1)

    except Exception as err:
        sys.exit(err)

    try:

        if os.path.exists(args.infile):
            infile = open(args.infile, "r")

    except Exception as err:
        sys.exit(err)

    ####  Read and check the keywords in INFILE

    player = Player(infile, outfile)

    player.read_keywords()

    player.check_keywords()
    player.print_keywords()

    if args.opt_continue:
        player.player.initcyc = save[0] + 1
        player.system = save[1]
    else:
        player.read_potential()

    ####  Check whether the executables are in the path
    ####		and print potential to Log File

    player.check_executables()

    player.print_potential()

    ####  Bring the molecules to standard orientation and prints info about them

    for i in range(len(player.system.molecule)):

        player.outfile.write(
            "\nMolecule type {} - {}:\n\n".format(
                i + 1, player.system.molecule[i].molname
            )
        )
        player.system.molecule[i].print_mol_info(player.outfile)
        player.outfile.write(
            "    Translating and rotating molecule to standard orientation..."
        )
        player.system.molecule[i].standard_orientation()
        player.outfile.write(" Done\n\n    New values:\n")
        player.system.molecule[i].print_mol_info(player.outfile)

    player.outfile.write(90 * "=")
    player.outfile.write("\n")

    if not args.opt_continue:
        make_simulation_dir()
    else:
        simdir = "simfiles"
        stepdir = "step{:02d}".format(player.player.initcyc)
        if os.path.exists(simdir + os.sep + stepdir):
            shutil.rmtree(simdir + os.sep + stepdir)

    ####  Open the geoms.xyz file and prints the initial geometry if starting from zero

    if player.player.initcyc == 1:
        try:
            path = "geoms.xyz"
            geomsfh = open(path, "w", 1)
        except EnvironmentError as err:
            sys.exit(err)
        player.system.print_geom(0, geomsfh)
        geomsfh.write(40 * "-" + "\n")
    else:
        try:
            path = "geoms.xyz"
            geomsfh = open(path, "a", 1)
        except EnvironmentError as err:
            sys.exit(err)

    player.outfile.write("\nStarting the iterative process.\n")

    ## Initial position (in Bohr)
    position = player.system.molecule[0].read_position()

    ## If restarting, read the last gradient and hessian
    # if player.player.initcyc > 1:
    # 	if player.player.qmprog in ("g03", "g09", "g16"):
    # 		Gaussian.read_forces("grad_hessian.dat")
    # 		Gaussian.read_hessian_fchk("grad_hessian.dat")

    # if player['qmprog'] == "molcas":
    # Molcas.read_forces("grad_hessian.dat")
    # Molcas.read_hessian("grad_hessian.dat")

    ###
    ###  Start the iterative process
    ###

    player.outfile.write("\n" + 90 * "-" + "\n")

    for cycle in range(
        player.player.initcyc, player.player.initcyc + player.player.maxcyc
    ):

        player.outfile.write("{} Step # {}\n".format(40 * " ", cycle))
        player.outfile.write(90 * "-" + "\n\n")

        make_step_dir(cycle)

        ####
        ####  Start block of parallel simulations
        ####

        player.dice_start(cycle)

        ###
        ###  End of parallel simulations block
        ###

        # Make ASEC
        player.outfile.write("\nBuilding the ASEC and vdW meanfields... ")
        asec_charges = player.populate_asec_vdw(cycle)

        ## After ASEC is built, compress files bigger than 1MB
        for proc in range(1, player.player.nprocs + 1):
            path = "step{:02d}".format(cycle) + os.sep + "p{:02d}".format(proc)
            compress_files_1mb(path)

        ###
        ###  Start QM calculation
        ###

        player.gaussian_start(cycle, geomsfh)

        player.system.print_geom(cycle, geomsfh)
        geomsfh.write(40 * "-" + "\n")

        player.outfile.write("\n+" + 88 * "-" + "+\n")

        pickle.dump([cycle, player.system], open("latest-step.pkl", "wb"))
    ####
    ####  End of the iterative process
    ####

    ## imprimir ultimas mensagens, criar um arquivo de potencial para ser usado em eventual
    ## continuacao, fechar arquivos (geoms.xyz, run.log, ...)

    player.outfile.write("\nDiceplayer finished normally!\n")
    player.outfile.close()
####
####  End of the program
####
