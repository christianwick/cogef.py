#!/usr/bin/env python3
""" Test all classes in molecule
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import numpy as np 
from cogef.molecule import Molecule

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
    # test mod_fragments in forward direction
    print("Testing mod_fragments in forward direction")
    mol = Molecule()
    mol.read_xyz(xyz_file)
    print(mol)
    print(f"Distance between H1 and H11 : {mol.distance(0,10)}")
    print(f"Centroid: {mol.centroid(mol.coordinates)}")
    print("centering Molecule")
    mol.center_coordinates()
    print(mol)
    print(f"Distance between H1 and H11 : {mol.distance(0,10)}")
    print(f"Centroid: {mol.centroid(mol.coordinates)}")
    mol2 = Molecule()
    mol2.read_xyz(xyz_file)
    print(f"RMSD : {mol.rmsd(mol2.coordinates,mol.coordinates)}")


