#!/usr/bin/env python3
""" gausp.py - read xyz files and write gaussian input for sp calculations
 
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
    parser.add_argument("-cm",help="charge/multiplicity line e.g '0,1'",type=str, default="0 1")
    parser.add_argument("-nproc", help="nproc", type=str, default="12")
    parser.add_argument("-mem", help="memory", type=str, default="12GB")
    parser.add_argument("-method", help="UwB97XD/def2TZVPP", type=str, default=None)
    parser.add_argument("-chk", help="use common chk file e.g. 'guess.chk'", type=str, default=None)
    parser.add_argument("-reverse", help="reverse calculation starting at last point", action="store_true", default=False)

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
    # start calculation    
    ginp = gaussian.GaussianInput(mem=args.mem, nproc=args.nproc, charge_multi = args.cm )
    glog = gaussian.CheckGaussianLogfile()
    route = gaussian.classGaussianRoute()
    route.route = ["#P", args.method.strip()]
    xyztrj = molecule.Molecule.read_xyz_trj(args.xyztrj)
    if args.reverse: 
        last = len(xyztrj) - 1
        xyztrj.reverse()
    else: last = 0
    for nn,xyz in enumerate(xyztrj):
        cycle = abs(last - ( nn  )) + 1
        filename = "sp_{:03d}".format(cycle)
        if args.chk: ginp.link0.chk=args.chk
        else: ginp.link0.chk=filename+".chk"
        ginp.molecule.read_xyz(xyz)
        with open(filename+".com", "w") as of:
            logger.info(f"Starting cycle {cycle}")
            ginp.of = of
            ginp.write_inputfile(link1=False, route=route, geom=True, modredundant=None)
        rungauss = subprocess.run(["cogef_rung16",filename + ".com"], check=True, capture_output=True, text=True)
        glog.filename = filename+".log"
        glog.read_log()
        ginp.molecule.write_xyz("gausp_trj.xyz", comment=glog.comment_line(point=cycle), write_mode="a")
        logger.info(f"Finished gaussian sp calculation at cycle {cycle}")

    # STOP LOGGING HERE
    logger.info("Finished gaussian single point calculations. " )