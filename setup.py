import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
# "packages": ["os"] is used as example only
build_exe_options = {"packages": ["os","setproctitle","nomkl","numpy"],
                     "bin_includes": ["/usr/bin/python3"],
                     "excludes": ["tkinter"],
                     "optimize": 0}

name = "diceplayer"

setup(
    name = name,
    version = "0.1",
    executables = [Executable("diceplayer/__main__.py", target_name=name)]
)