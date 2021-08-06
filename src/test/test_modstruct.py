#!/usr/bin/env python3
""" Test all classes in modstruct
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import numpy as np 
from cogef.gaussian import GaussianInput
from cogef.modstruct import mod_fragments, mod_single_atom, mod_atoms

xyz_file="""11
E(RB3LYP) = -119.164345992 on  Stat Point 2
H      -2.179023     0.366411     0.000000 
C      -1.279656    -0.261298     0.000000 
H      -1.321768    -0.908272     0.885321 
H      -1.321768    -0.908272    -0.885321 
C       0.000048     0.589806     0.000000 
H      -0.000057     1.248529     0.878918 
H      -0.000057     1.248529    -0.878918 
C       1.279627    -0.261297     0.000000 
H       1.321758    -0.908281     0.885302 
H       1.321758    -0.908281    -0.885302 
H       2.179042     0.366372     0.000000 
"""

if __name__ == "__main__":
    coords = GaussianInput()
    coords.molecule.read_xyz(xyz_file)
    counter = 0
    while counter < 10:
        coords.molecule.coordinates = mod_fragments(coords = coords.molecule.coordinates, atom1 = 0, atom2 = 10, dx = 0.5, dp = 1.1, symmetric=True)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1

    coords.molecule.read_xyz(xyz_file)
    counter = 0
    while counter < 10:
        coords.molecule.coordinates = mod_single_atom(coords = coords.molecule.coordinates, atom1 = 0, atom2 = 10, dx = 0.5, symmetric=True)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1
