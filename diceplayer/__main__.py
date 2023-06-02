from diceplayer.shared.interface.dice_interface import DiceInterface
from diceplayer.player import Player
from diceplayer import logger

from pathlib import Path
import argparse
import logging
import pickle
import sys

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
        "-i", "--input",
        dest="infile",
        default="control.yml",
        metavar="INFILE",
        help="input file of diceplayer [default = control.in]"
    )
    parser.add_argument(
        "-o", "--output",
        dest="outfile",
        default="run.log",
        metavar="OUTFILE",
        help="output file of diceplayer [default = run.log]"
    )
    args = parser.parse_args()

    # Open OUTFILE for writing and print keywords and initial info
    logger.set_logger(args.outfile, logging.INFO)

    try:

        pickle_path = Path("latest-step.pkl")
        if args.opt_continue and pickle_path.exists():
            with open(pickle_path) as pickle_file:
                save = pickle.load(pickle_file)

    except Exception as err:
        sys.exit(err)

    player = Player(args.infile)

    player.print_keywords()

    player.create_simulation_dir()

    player.read_potentials()
    player.print_potentials()

    player.prepare_system()

    player.start()

    logger.info("\n+" + 88 * "-" + "+\n")

    player.print_results()

    logger.info("\n+" + 88 * "-" + "+\n")

    logger.info("Diceplayer finished successfully \n")


if __name__ == "__main__":
    main()
