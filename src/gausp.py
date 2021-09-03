#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""

import sys
import numpy as np
#import matplotlib.pyplot as plt
import csv
import json
import argparse
import subprocess
import logging

from cogef.constants import *
from cogef import modstruct
from cogef import gaussian
from cogef import molecule
from cogef import cogef_logging

from cogef._version import __version__ 



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-xyztrj", required=True, help="xyz trajectory file", type=argparse.FileType('r'))
    parser.add_argument("-charge",help="list of charges",type=int, nargs="+", default=[0])
    parser.add_argument("-multiplicity",help="list of multiplicities",type=int, nargs="+", default=[1])
    parser.add_argument("-nproc", help="nproc", type=str, default="12")
    parser.add_argument("-mem", help="memory", type=str, default="12GB")
    parser.add_argument("-method", help="UwB97XD/def2TZVPP", type=str, default=None)
    parser.add_argument("-dispersion", help="GD3/GD3BJ/GD3", type=str, choices=["GD3","GD3BJ"], default=None)

    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO"], default="INFO")

    args = parser.parse_args()

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file="gausp.log",log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("gausp")
    logger.info("Starting Gaussian single point calculations .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))

    data = gaussian.GaussianInput(mem=args.mem, nproc=args.nproc, charge=args.charge, multiplicity=args.multiplicity)
    if args.method:
        data.route["level_of_theory"] = args.method 
    if args.dispersion:
        data.route["empiricalDispersion"] = args.dispersion

    xyztrj = molecule.Molecule.read_xyz_trj(args.xyztrj)
    for nn,xyz in enumerate(xyztrj):
        filename = "sp_{:03d}".format(nn)
        data.molecule.read_xyz(xyz)
        with open(filename+".com", "w") as of:
            data.write_sp_input(of)
        rungauss = subprocess.run(["cogef_rung16",filename + ".com"], check=True, capture_output=True, text=True)

    # STOP LOGGING HERE
    logger.info("Finished gaussian single point calculations. " )