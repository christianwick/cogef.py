#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""

import sys
import os
import numpy as np
#import matplotlib.pyplot as plt
import csv
import json
import argparse
import subprocess

from cogef.constants import *
from cogef import modstruct
from cogef import gaussian


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-xyz", required=True, help="xyz file",type=argparse.FileType('r'))
    parser.add_argument("-dx", help="increment", default=0.02,type=float)
    parser.add_argument("-dp", help="add perturbation for fragment atoms", default=1.0,type=float)
    parser.add_argument("atoms", help="atom1 atom2 starting at 1",nargs='+', type=int)
    parser.add_argument("-cycles", help="number n of cycles to run", type=int, default=1)
    parser.add_argument("-max_instab_cycles", help="number n of cycles to rerun after instability has been found", type=int, default=5)
    parser.add_argument("-max_stationary_cycles", help="number n of cycles to rerun after no stationary point has been found", type=int, default=5)
    parser.add_argument("-restart", help="restart at cycle n", type=int, default=0)
    parser.add_argument("-symm", help="move atoms symmetrically", action="store_true", default=False)
    parser.add_argument("-new_guess", help="always generate a new guess", action="store_true", default=False)
    parser.add_argument("-no_opt", help="do not use optimized geometry in guess structure", action="store_true", default=False)
    parser.add_argument("-fragment",help="list of atoms in fragment 1",type=int, nargs="+", default=[1,2,3,4])
    parser.add_argument("-fragment_guess", help="always generate a fragment guess", action="store_true", default=False)
    parser.add_argument("-charge",help="list of charges",type=int, nargs="+", default=[0])
    parser.add_argument("-multiplicity",help="list of multiplicities",type=int, nargs="+", default=[1])
    parser.add_argument("-nproc", help="nproc", type=str, default="12")
    parser.add_argument("-mem", help="memory", type=str, default="12GB")
    parser.add_argument("-method", help="UwB97XD/def2TZVPP", type=str, default=None)
    parser.add_argument("-calcfc", help="calculate Force Constants for difficult cases", choices=["CalcFC","RecalcFC=5","CalcAll", "None"], default="CalcFC")
    parser.add_argument("-rfo", help="use rfo",  choices=["RFO", "None"], default="RFO")
    parser.add_argument("-dispersion", help="GD3/GD3BJ/GD3", type=str, choices=["GD3","GD3BJ"], default=None)

    args = parser.parse_args()
    atom1 = args.atoms[0]-1
    atom2 = args.atoms[1]-1

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

    data.molecule.read_xyz(args.xyz)
    instab=False
    stationary=False
    for ii in range(args.restart, args.cycles):
        stationary_cycles = 0
        while stationary_cycles <= args.max_stationary_cycles:
            instab_cycles = 0
            while instab_cycles <= args.max_instab_cycles:
                filename="cogef_{:03d}".format(int(ii))
                with open(filename+".com", "w") as of:
                    data.write_input(of,modredundant=["{} {} F".format(args.atoms[0],args.atoms[1])], 
                          initial_stab_opt = not stationary,
                          instability = instab,
                          route_args = opt_args)
                subprocess.run(["cogef_rung16",filename+".com"],check=True)
                xyz = subprocess.run(["gxyz",filename+".log","-l","--input"], capture_output=True, text=True)
                instab = gaussian.read_instability_from_log(filename+".log")
                stationary = gaussian.read_stationary_from_log(filename+".log")
                if not args.no_opt:
                    data.molecule.read_xyz(xyz.stdout)
                    #args.no_opt = False
                # check for instability or break out 
                if instab:
                    os.rename(filename+".log", "{}_{}_{}.log".format(filename, stationary_cycles, instab_cycles))
                    instab_cycles +=1
                else:
                    break
            # check for stationary point or break out 
            if not stationary:
                os.rename(filename+".log", "{}_{}_{}.log".format(filename, stationary_cycles, instab_cycles))
                stationary_cycles += 1
            else:
                break           
        data.molecule.coordinates = modstruct.mod_fragments(data.molecule.coordinates, atom1, atom2, args.dx, symmetric=args.symm, dp=args.dp, fragment=args.fragment)
        data.molecule.write_xyz(sys.stdout)