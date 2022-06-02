import argparse
from distutils.command.clean import clean
import os
import shutil
import sys

import PyInstaller.__main__

name = "diceplayer"

parser = argparse.ArgumentParser(prog="Diceplayer Setup")

parser.add_argument(
    "-b", "--build",
    dest="build",
    default=False,
    action="store_true",
    help="Builds the Diceplayer Binary",
)
parser.add_argument(
    "-c", "--clean",
    dest="clean",
    default=False,
    action="store_true",
    help="Cleans the Development Environment",
)

args = parser.parse_args()


def __build():

    PyInstaller.__main__.run(
        ["diceplayer/__main__.py", "--onefile", "-n{}".format(name)]
    )


def __clean():

    shutil.rmtree("build")
    shutil.rmtree("dist")
    os.remove("diceplayer.spec")


if __name__ == "__main__":

    if args.build:
        __build()
    elif args.clean:
        __clean()
    else:
        parser.print_help(sys.stderr)
        sys.exit(1)
