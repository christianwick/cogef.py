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
    parser.add_argument("-template", required=True, help="orca template", type=argparse.FileType('r'))
    parser.add_argument("-use_last_gbw", help="use last gbw file as input", action="store_true", default=False)
    parser.add_argument("-guess_orbs", help="use in guess orbs 'guess_orbs.gbw'", type=str, default="guess_orbs.gbw")

    args = parser.parse_args()

    xyztrj = molecule.Molecule.read_xyz_trj(args.xyztrj)
    xyztrj.reverse()
    num_mols = len(xyztrj)
    current = num_mols - 1
    template = args.template.read()
    for nn,xyz in enumerate(xyztrj):
        last = current 
        current = num_mols - nn - 1
        filename = "sp_{:03d}".format(current)
        mol = molecule.Molecule()
        mol.read_xyz(xyz)
        mol.write_xyz(filename + ".xyz")
        orca_inp = template.replace("REPLACE",filename)
        if args.use_last_gbw:
            if nn == 0:
                gbw_guess = args.guess_orbs
            else:
                gbw_guess = "sp_{:03d}.gbw".format(last)
            orca_inp = orca_inp.replace("LASTGBW",gbw_guess)
        with open(filename+".inp", "w") as of:
            of.write(orca_inp)
        subprocess.run(["runorca_4_2",filename + ".inp"], check=True)
