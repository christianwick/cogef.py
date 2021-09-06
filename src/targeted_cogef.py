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
import subprocess

from cogef.constants import *
from cogef import modstruct
from cogef import gaussian
from cogef import cogef_logging
from cogef._version import __version__ 

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("atoms", help="atom1 atom2 starting at 1",nargs='+', type=int)
    
    group_molecule = parser.add_argument_group("molecule definition")
    group_molecule .add_argument("-xyz", required=True, help="xyz file",type=argparse.FileType('r'))
    group_molecule.add_argument("-charge",help="list of charges",type=int, nargs="+", default=[0])
    group_molecule.add_argument("-multiplicity",help="list of multiplicities",type=int, nargs="+", default=[1])

    group_gaussian = parser.add_argument_group("Gaussian options")
    group_gaussian.add_argument("-nproc", help="nproc", type=str, default="12")
    group_gaussian.add_argument("-mem", help="memory", type=str, default="12GB")
    group_gaussian.add_argument("-method", help="UwB97XD/def2TZVPP", type=str, default=None)
    group_gaussian.add_argument("-dispersion", help="GD3/GD3BJ/GD3", type=str, choices=["GD3","GD3BJ"], default=None)
    group_gaussian.add_argument("-calcfc", help="calculate Force Constants for difficult cases", choices=["CalcFC","RecalcFC=5","CalcAll", "None"], default="CalcFC")
    group_gaussian.add_argument("-rfo", help="use rfo",  choices=["RFO", "None"], default="RFO")

    group_output = parser.add_argument_group("Output options")
    group_output.add_argument("-trajectory", help="write trajectory to file", type=argparse.FileType("a"), default=None)

    group_tweaks = parser.add_argument_group("cogef")
    group_tweaks.add_argument("-dx", help="increment", default=0.02,type=float)
    group_tweaks.add_argument("-dp", help="add perturbation for fragment atoms", default=1.0,type=float)
    group_tweaks.add_argument("-fragment",help="list of atoms in fragment 1",type=int, nargs="+", default=[])
    group_tweaks.add_argument("-cycles", help="number n of cycles to run", type=int, default=1)
    group_tweaks.add_argument("-reverse", help="reverse direction start at n", type=int, default=0)
    group_tweaks.add_argument("-max_instab_cycles", help="number n of cycles to rerun after instability has been found", type=int, default=5)
    group_tweaks.add_argument("-max_stationary_cycles", help="number n of cycles to rerun after no stationary point has been found", type=int, default=5)
    group_tweaks.add_argument("-restart", help="restart at cycle n. (Note:"+
        "If you want to restart in reverse direction use the reverse option!)", type=int, default=0)
    group_tweaks.add_argument("-symm", help="move atoms symmetrically", action="store_true", default=False)
    group_tweaks.add_argument("-modred", help="modredundant section, separated by ';' ", default=None, type=str)

    group_deprecated = parser.add_argument_group("deprecated")
    group_deprecated.add_argument("-no_opt", help="!DEPRECATED! do not use optimized geometry in guess structure", action="store_true", default=False)
    group_deprecated.add_argument("-fragment_guess", help="always generate a fragment guess", action="store_true", default=False)

    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO"], default="INFO")

    parser.add_argument('-version', action='version', version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    atom1 = args.atoms[0]-1
    atom2 = args.atoms[1]-1

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("targeted_cogef")
    logger.info("Starting Targeted Cogef .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))

    args.fragment = [ int(x) - 1 for x in args.fragment ]
    data = gaussian.GaussianInputWithFragments(mem=args.mem, nproc=args.nproc, charge=args.charge, multiplicity=args.multiplicity)
    if args.method:
        data.route["level_of_theory"] = args.method 
    if args.dispersion:
        data.route["empiricalDispersion"] = args.dispersion

    # set opt options
    opt_args = { "opt": ["MaxCyc=100"] } 
    if args.calcfc != "None":
        opt_args["opt"].append(args.calcfc)
    if args.rfo != "None":
        opt_args["opt"].append("RFO")
    logger.info("Setting opt options: {}".format(opt_args))

    # set mod redundant input
    modredundant = ["{} {} F".format(args.atoms[0],args.atoms[1])]
    if args.modred:
        for section in args.modred.split(";"):
            modredundant.append(section)

    data.molecule.read_xyz(args.xyz)
    logger.info("Initial Coordinates: \n" + str(data.molecule))

    # to check for spin contamination 
    found_spin_contamination = False

    instab=False
    stationary=False
    # cogef main loop
    if args.reverse:
        args.dx = - abs(args.dx)
        # we loop from reverse n cycles backwards. smallest point is 000
        loop_range = range(args.reverse, max(-1, args.reverse - args.cycles -1 ), -1)
    else:
        loop_range = range(args.restart, args.cycles +1 )
    for ii in loop_range:
        logger.info("Starting cycle {}".format(ii))
        stationary_cycles = 0
        while stationary_cycles <= args.max_stationary_cycles:
            instab_cycles = 0
            while instab_cycles <= args.max_instab_cycles:
                filename="cogef_{:03d}".format(int(ii))
                with open(filename+".com", "w") as of:
                    data.write_input(of,modredundant=modredundant, 
                          initial_stab_opt = not stationary,
                          instability = instab,
                          route_args = opt_args)
                rungauss = subprocess.run(["cogef_rung16",filename+".com"], check=True, capture_output=True, text=True)
                # read gaussian log file
                glog = gaussian.CheckGaussianLogfile(filename+".log")
                glog.read_log()
                if not args.no_opt:
                    data.molecule.elements = glog.molecule.elements
                    data.molecule.coordinates = glog.molecule.coordinates
                    #args.no_opt = False
                # check for instability or break out 
                if glog.instability:
                    logger.warning("Instability detected at cycle {} {} {}".format(ii, stationary_cycles, instab_cycles, ))
                    os.rename(filename+".log", "{}_{}_{}.log".format(filename, stationary_cycles, instab_cycles))
                    instab_cycles +=1
                elif glog.error:
                    # at the moment, we will treat error terminations as instability and try another instab cycle.
                    logger.warning("Error termination detected at cycle {} {} {}".format(ii, stationary_cycles, instab_cycles, ))
                    os.rename(filename+".log", "{}_{}_{}.log".format(filename, stationary_cycles, instab_cycles))
                    instab_cycles += 1
                else:
                    logger.info("Wavefunction is stable.")
                    break
            # check for stationary point or break out 
            if not glog.stationary:
                logger.warning("No Stationary point found at cycle {} {} {}".format(ii, stationary_cycles, instab_cycles, ))
                os.rename(filename+".log", "{}_{}_{}.log".format(filename, stationary_cycles, instab_cycles))
                stationary_cycles += 1
            else:
                logger.info("Optimisation completed.")
                logger.debug("Optimised Coordinates: \n" + str(data.molecule))
                # write actual structure to disk.
                logger.info(f"Writing checkpoint.xyz")
                data.molecule.write_xyz(f"checkpoint.xyz",comment=glog.comment_line(point=ii))
                if args.trajectory:
                    logger.info(f"Adding structure to trajectory {args.trajectory}")
                    data.molecule.write_xyz(args.trajectory, comment=glog.comment_line(point=ii))
                # we check for spin contamination and write the first structure with spin contamination to disk.
                if not found_spin_contamination and glog.spin > 0.5 and not args.reverse: 
                    logger.info("Detected homolytic bond scission. Saving structure to disk...")
                    data.molecule.write_xyz(f"start_reverse_{ii:03d}.xyz", comment=glog.comment_line(point=ii))
                    found_spin_contamination = True
                break           
        data.molecule.coordinates = modstruct.mod_fragments(data.molecule.coordinates, atom1, atom2, args.dx, symmetric=args.symm, dp=args.dp, fragment=args.fragment)
        #data.molecule.write_xyz(sys.stdout)
        logger.debug("New Coordinates: \n" + str(data.molecule))

    # STOP LOGGING HERE
    logger.info("Summary of cogef calculation:" )
    logger.info("Total cycles: {}".format(abs(ii - args.reverse - args.restart )))
    logger.info("Finished Manual Cogef. " )