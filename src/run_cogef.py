#!/usr/bin/env python3
""" targeted_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""

import sys
import os
import numpy as np
import logging 

import csv
import json
import argparse


from cogef.constants import *
from cogef import modstruct
from cogef import gaussian
from cogef import cogef_logging
from cogef import driver as driver
from cogef._version import __version__ 

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-at", required=True, help="atom1 atom2 starting at 1", type=str)
    
    group_molecule = parser.add_argument_group("molecule definition")
    group_molecule.add_argument("-xyz", required=True, help="xyz file",type=argparse.FileType('r'))
    group_molecule.add_argument("-cm",help="charge and multiplictiy line",type=str, default="0 1")
    



    group_gaussian = parser.add_argument_group("Gaussian options")
    group_gaussian.add_argument("-nproc", help="nproc", type=str, default="12")
    group_gaussian.add_argument("-mem", help="memory", type=str, default="12GB")
    group_gaussian.add_argument("-method", help="B3LYP/6-31G*", type=str, default=None)
    group_gaussian.add_argument("-read_guess", help="read in guess from previous checkpoint file", action="store_true", default=False)

    group_output = parser.add_argument_group("Output options")
    group_output.add_argument("-trajectory", help="write trajectory to file", type=argparse.FileType("a"), default=None)
    group_output.add_argument("-checkpoint", help="name of checkpoint file. None if no checkpointing",
            type=str, default="checkpoint.xyz")

    group_tweaks = parser.add_argument_group("cogef")
    group_tweaks.add_argument("-runtype", help="restricted/unrestricted cogef", 
            choices=["restricted"],type=str,default="restricted")
    group_tweaks.add_argument("-dx", help="increment", default=0.02,type=float)
    group_tweaks.add_argument("-dp", help="add perturbation for fragment atoms", default=1.0,type=float)
    group_tweaks.add_argument("-fragment",help="list of atoms in fragment 1 as string, e.g. '1 2 3 4'",type=str, default="")
    group_tweaks.add_argument("-cycles", help="number n of cycles to run", type=int, default=1)
    
    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-logfile", help="name of the logfile", type=str, default="job_cogef.log")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO"], default="INFO")

    args = parser.parse_args()
    atom1 ,atom2 = [ int(x)-1 for x in args.at.split() ]

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file=args.logfile, log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("run_cogef")
    logger.info("Starting Cogef calculation .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))

    args.fragment = [ int(x) - 1 for x in args.fragment.split() ]

    driver_args = { "runtype" : args.runtype, "atom1" : atom1 , "atom2" : atom2, "dx" : 0.1, 
                    "level_of_theory" : args.method , "read_guess" : args.read_guess, "cycles" : args.cycles, 
                    "reverse" : None, "restart" : 0, "modredundant" : None, "symmetric" : True, 
                    "dp" : 1.0, "fragment" : [], "max_error_cycles"  : 5, 
                    "trajectory" : args.trajectory, "checkpoint" : args.checkpoint }
    data = driver.cogef_loop(**driver_args)
    data.ginp.molecule.read_xyz(args.xyz)
    data.run()
    


    logger.info("Finished Manual Cogef. " )