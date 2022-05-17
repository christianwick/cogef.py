#!/usr/bin/env python3
""" run_cogef.py - perform COGEF calculations at various levels of theory

 
"""

import sys
import os
import numpy as np
import logging 

import csv
import json
import argparse

from cogef import misc
from cogef import cogef_logging
from cogef import driver as driver
from cogef._version import __version__ 

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', action='version',version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument("-at", required=True, help="atom1 atom2 starting at 1 as list eg. '1,2'", type=str)
    
    group_molecule = parser.add_argument_group("molecule definition")
    group_molecule.add_argument("-xyz", required=True, help="xyz file",type=argparse.FileType('r'))
    group_molecule.add_argument("-cm",help="charge and multiplictiy line",type=str, default="0 1")
    
    group_oniom = parser.add_argument_group("ONIOM/AMBER calculations")
    group_oniom.add_argument("-oniom", help="oniom template file",type=argparse.FileType('r'))
    group_oniom.add_argument("-amber", help="amber template file",type=argparse.FileType('r'))
    group_oniom.add_argument("-oniomopt", help="oniom opt method: hybrid,me or ee_sp",choices=["hybrid","me","ee","ee_sp"],type=str)
    group_oniom.add_argument("-constraint", help="oniom constraints: modredundant or opt_flag",
            choices=["modredundant","opt_flag"], type=str, default="modredundant")
    group_oniom.add_argument("-no_micro", help="do not use microiterations in ONIOM", action="store_true", default=False)
    group_oniom.add_argument("-quadmac", help="do use the QuadMac option during microiterations in ONIOM",
        action="store_true", default=False)


    group_gaussian = parser.add_argument_group("Gaussian options")
    group_gaussian.add_argument("-nproc", help="nproc", type=str, default=None)
    group_gaussian.add_argument("-mem", help="memory", type=str, default=None)
    group_gaussian.add_argument("-method", help="e.g. B3LYP/6-31G*", type=str, default=None)
    group_gaussian.add_argument("-startchk", help="use initial guess from guess.chk", action="store_true", default=False)
    group_gaussian.add_argument("-maxcyc", help="Sets the maximum number of optimization steps to N", type=int, default=50)
    group_gaussian.add_argument("-maxconv", help="Sets the maximum number conventional SCF cycles", type=int, default=75)
    group_gaussian.add_argument("-modredundant", help="add additional modredundant sections separated by ',' eg. '1 2 F, 2 3 F'", 
        type=str, default="")
    group_gaussian.add_argument("-readfc", help="use initial ForceConstants from guess.chk", action="store_true", default=False)

    group_output = parser.add_argument_group("Output options")
    group_output.add_argument("-trajectory", help="write trajectory to file", type=str, default=None)
    group_output.add_argument("-checkpoint", help="name of the checkpoint file.",
            type=str, default="checkpoint.xyz")

    group_tweaks = parser.add_argument_group("cogef")
    group_tweaks.add_argument("-runtype", help="restricted/unrestricted/rTS/uTS cogef", 
            choices=["restricted","unrestricted","rTS","uTS"], type=str, default="restricted")
    group_tweaks.add_argument("-cguess", help="None/dx/strain", 
            choices=["dx","strain","new_strain"], type=str, default=None)
    group_tweaks.add_argument("-dx", help="step size.", default=0.02,type=float)
    group_tweaks.add_argument("-alpha", help="add perturbation for fragment atoms", default=1.0,type=float)
    group_tweaks.add_argument("-beta", help="damp the additive strain perturbation", default=1.0,type=float)
    group_tweaks.add_argument("-nbonds", help="damp the additive strain perturbation according to the number of bonds nbonds",
                                default=1.0,type=float)
    group_tweaks.add_argument("-frag1",help="list of atoms in fragment 1 as string, e.g. '1,2,3,4' or '1-4'",type=str, default="")
    group_tweaks.add_argument("-frag2",help="list of atoms in fragment 2 as string, e.g. '1,2,3,4' or '1-4'",type=str, default="")
    group_tweaks.add_argument("-exclude",help="list of atoms to exclude from fragments as string, e.g. '1,2,3,4' or '1-4'",type=str, default="")
    group_tweaks.add_argument("-cycles", help="number n of cycles to run", type=int, default=1)
    group_tweaks.add_argument("-restart", help="restart from cylce n. needs guess.chk and checkpoint.xyz", type=int, default=1)
    group_tweaks.add_argument("-restart_xyz", help="restart coordinates from checkpoint.xyz", type=argparse.FileType('r'), default=None)
    group_tweaks.add_argument("-reverse", help="reverse from cylce n. needs guess.chk and start geometry as .xyz", type=int, default=None)
    group_tweaks.add_argument("-max_error_cycles", help="number n of error cycles to run", type=int, default=5)
    group_tweaks.add_argument("-symm_stretch", 
        help ="perform symmetric stretch moving both atoms (framgents). asymmetric otherwise.", action="store_true",
        default = False)
    group_tweaks.add_argument("-mulliken_h", 
        help ="use mulliken charges with Hydrogens summed into heavy atoms. Otherwise use standard mulliken charges",
        action="store_true", default = False)
    group_tweaks.add_argument("-no_mix", help ="Turn of mixing of HOMO LUMO orbitals completely.",
        action="store_true", default = False)
    
    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-logfile", help="name of the logfile", type=str, default="job_cogef.log")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")

    args = parser.parse_args()
    atom1 ,atom2 = [ int(x)-1 for x in args.at.split(",") ]
    args.frag1 = misc.convert_inputstr(args.frag1)
    args.frag2 = misc.convert_inputstr(args.frag2)
    args.exclude = misc.convert_inputstr(args.exclude)
    args.modredundant = misc.convert_modredundant(args.modredundant) 

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file=args.logfile, log_level = args.log_level, stream_level = args.log_level)
    logger = logging.getLogger("run_cogef")
    logger.info("Starting Cogef calculation .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        # we print out all options. internal numbering starts at 0, input at 1!
        # Therefore, we fix the printout for the logfile.
        if arg == "fragment": logger.info("    -{:10s} : {}".format(arg,[ x+1 for x in value]))
        elif arg == "exclude": logger.info("    -{:10s} : {}".format(arg,[ x+1 for x in value]))
        else: logger.info("    -{:10s} : {}".format(arg,value))


    driver_args = { "xyz" : args.xyz, "runtype" : args.runtype, "atom1" : atom1 , "atom2" : atom2, "dx" : args.dx, 
                    "level_of_theory" : args.method , "mem" : args.mem, "nproc" : args.nproc,
                    "maxcyc" : args.maxcyc, "maxconv" : args.maxconv, "cm" : args.cm ,
                    "startchk" : args.startchk, "cycles" : args.cycles, 
                    "reverse" : args.reverse, "restart" : args.restart, "restart_xyz" : args.restart_xyz,
                    "modredundant" : args.modredundant, "readfc" : args.readfc,
                    "symm_stretch" : args.symm_stretch, "alpha" : args.alpha, "frag1" : args.frag1, "frag2" : args.frag2, "exclude" : args.exclude,
                    "max_error_cycles"  : args.max_error_cycles, "mulliken_h" : args.mulliken_h, "trajectory" : args.trajectory,
                    "checkpoint" : args.checkpoint, "oniomtemplate" : args.oniom, "ambertemplate" : args.amber,
                    "oniomopt" : args.oniomopt, "constraint" : args.constraint ,
                    "no_mix" : args.no_mix, "no_micro" : args.no_micro, "quadmac" : args.quadmac ,
                    "beta" : args.beta, "nbonds" : args.nbonds,
                    "cguess" : args.cguess }
    if args.oniom:
        logger.debug("Starting ONIOM driver")
        data = driver.oniom_cogef_loop(**driver_args)
    elif args.amber:
        logger.debug("Starting AMBER driver")
        data = driver.amber_cogef_loop(**driver_args)
    else:
        logger.debug("Starting driver")
        data = driver.cogef_loop(**driver_args)
    data.run()
    
    logger.info("Finished cogef calculation. " )