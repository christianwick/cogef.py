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

from cogef.constants import *
from cogef import modstruct
from cogef import gaussian
from cogef import molecule





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-xyztrj", required=True, help="xyz trajectory file", type=argparse.FileType('r'))
    parser.add_argument("-charge",help="list of charges",type=int, nargs="+", default=[0])
    parser.add_argument("-multiplicity",help="list of multiplicities",type=int, nargs="+", default=[1])
    parser.add_argument("-nproc", help="nproc", type=str, default="12")
    parser.add_argument("-mem", help="memory", type=str, default="12GB")
    parser.add_argument("-method", help="UwB97XD/def2TZVPP", type=str, default=None)
    parser.add_argument("-dispersion", help="GD3/GD3BJ/GD3", type=str, choices=["GD3","GD3BJ"], default=None)

    args = parser.parse_args()

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
        subprocess.run(["cogef_rung16",filename + ".com"], check=True)
