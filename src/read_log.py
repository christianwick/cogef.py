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
    parser.add_argument("glogs", help="gaussian log files", type=str, nargs="+")

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
    glog = gaussian.CheckGaussianLogfile()
    for cycle, glogfile in enumerate(args.glogs):
        glog.filename = glogfile
        glog.read_log()
        logger.info(f"imaginary frequencies: {glog.imag_frequencies} {len(glog.imag_frequencies)}")
        logger.info(f"lowest freq {glog.lowest_freq}")
        glog.molecule.write_xyz("glog_trj.xyz", comment=glog.comment_line(point=cycle+1), write_mode="a")

    # STOP LOGGING HERE
    logger.info("Finished reading gaussian log files " )