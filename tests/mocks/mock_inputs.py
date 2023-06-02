from unittest import mock


def get_config_example():
    return """
diceplayer:
    opt: no
    mem: 12
    maxcyc: 3
    ncores: 4
    nprocs: 4
    qmprog: 'g16'
    lps: no
    ghosts: no
    altsteps: 20000

    dice:
        nmol: [1, 50]
        dens: 0.75
        nstep: [2000, 3000, 4000]
        isave: 1000
        outname: 'phb'
        progname: '~/.local/bin/dice'
        ljname: 'phb.ljc'
        randominit: 'first'

    gaussian:
        qmprog: 'g16'
        level: 'MP2/aug-cc-pVDZ'
        keywords: 'freq'
"""


def get_potentials_exemple():
    return """\
*
2
1 TEST
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def get_potentials_error_combrule():
    return """\
.
2
1 TEST
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def get_potentials_error_ntypes():
    return """\
*
a
1 TEST
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def get_potentials_error_ntypes_config():
    return """\
*
3
1 TEST
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def get_potentials_error_nsites():
    return """\
*
2
. TEST
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def get_potentials_error_molname():
    return """\
*
2
1
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
1 PLACEHOLDER
 1   1   0.000000   0.000000   0.000000   0.000000  0.0000  0.0000
"""


def mock_open(file, *args, **kwargs):
    values = {
        "control.test.yml": get_config_example(),
        "phb.ljc": get_potentials_exemple(),
        "phb.error.combrule.ljc": get_potentials_error_combrule(),
        "phb.error.ntypes.ljc": get_potentials_error_ntypes(),
        "phb.error.ntypes.config.ljc": get_potentials_error_ntypes_config(),
        "phb.error.nsites.ljc": get_potentials_error_nsites(),
        "phb.error.molname.ljc": get_potentials_error_molname(),
    }
    mock_file = mock.mock_open(read_data=values[file])
    return mock_file()