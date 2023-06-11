from diceplayer import logger
from diceplayer.player import Player
from diceplayer.shared.interface.dice_interface import DiceInterface

import argparse
import logging
import pickle
import sys
from pathlib import Path

__VERSION = "v0.0.1"


def main():
    """
    Read and store the arguments passed to the program
    and set the usage and help messages
    """

    parser = argparse.ArgumentParser(prog="Diceplayer")
    parser.add_argument(
        "-c", "--continue", dest="opt_continue", default=False, action="store_true"
    )
    parser.add_argument(
        "-v", "--version", action="version", version="diceplayer-" + __VERSION
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="infile",
        default="control.yml",
        metavar="INFILE",
        help="input file of diceplayer [default = control.in]",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="outfile",
        default="run.log",
        metavar="OUTFILE",
        help="output file of diceplayer [default = run.log]",
    )
    args = parser.parse_args()

    # Open OUTFILE for writing and print keywords and initial info
    logger.set_logger(args.outfile, logging.INFO)

    if args.opt_continue:
        player = Player.from_save()
    else:
        player = Player.from_file(args.infile)

        player.read_potentials()

        player.create_simulation_dir()
        player.create_geoms_file()

    player.print_keywords()

    player.print_potentials()

    player.prepare_system()

    player.start()

    logger.info("\n+" + 88 * "-" + "+\n")

    player.print_results()

    logger.info("\n+" + 88 * "-" + "+\n")

    logger.info("Diceplayer finished successfully \n")


if __name__ == "__main__":
    main()
