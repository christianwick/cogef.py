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
    # setup argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-xyztrj", required=True, help="xyz trajectory file", type=argparse.FileType('r'))
    parser.add_argument("-template", required=True, help="orca template", type=argparse.FileType('r'))
    parser.add_argument("-use_last_gbw", help="use last gbw file as input", action="store_true", default=False)
    parser.add_argument("-guess_orbs", help="use in guess orbs 'guess_orbs.gbw'", type=str, default="guess_orbs.gbw")
    parser.add_argument("-start_at", help="start at number n (-1 = last, default)", type=int, default=-1)
    parser.add_argument("-direction", help="run in reverse direction (starting at last xyz) or not", type=str, 
        choices=["REVERSE","FORWARD"], default="REVERSE")

    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO"], default="INFO")

    args = parser.parse_args()


    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file="simple_orca_sp.log", log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("simple_orca_sp")
    logger.info("Starting orca single point calculations .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))
    


    # read in xyz trajectory 
    xyztrj = molecule.Molecule.read_xyz_trj(args.xyztrj)
    
    # we set the last structure
    if args.start_at == -1 :
        args.start_at = len(xyztrj) -1
        print(len(xyztrj) ,args.start_at)
    # check the direction and strip trailing coordinates
    if args.direction == "REVERSE":
        xyztrj = xyztrj[:args.start_at+1]
        print(len(xyztrj))
        xyztrj.reverse()   
    else:
        xyztrj = xyztrj[args.start_at:]
        
    current = args.start_at
    # read the template file and start the main loop
    template = args.template.read()
    for nn,xyz in enumerate(xyztrj):
        # set last sp number and check direction
        last = current 
        if args.direction == "REVERSE":
            current = args.start_at - nn  
        else:
            current = args.start_at + nn
        logger.info(f"starting cycle {current:03d}")
        # setup file name, read xyz coordinates and store them in a single file 
        filename = "sp_{:03d}".format(current)
        mol = molecule.Molecule()
        mol.read_xyz(xyz)
        mol.write_xyz(filename + ".xyz")
        # replace STRINGS in orca template
        orca_inp = template.replace("REPLACE",filename)
        if args.use_last_gbw:
            if nn == 0:
                gbw_guess = args.guess_orbs
            else:
                gbw_guess = "sp_{:03d}.gbw".format(last)
            orca_inp = orca_inp.replace("LASTGBW",gbw_guess)
        # write orca input file
        with open(filename+".inp", "w") as of:
            of.write(orca_inp)
        # run orca 
        runorca = subprocess.run(["runorca_4_2",filename + ".inp"], check=True, capture_output=True, text=True)

    # STOP LOGGING HERE
    logger.info("Finished orca single point calculations. " )