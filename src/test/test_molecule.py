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

rotated_xyz="""11

H       1.47894      1.28058     -1.02406
C       1.15875      0.39713     -0.45848
H       1.42985     -0.49164     -1.04221
H       1.73517      0.36890      0.47483
C      -0.35101      0.42917     -0.17280
H      -0.90059      0.48889     -1.12191
H      -0.59748      1.34321      0.38416
C      -0.83193     -0.79674      0.61939
H      -0.62635     -1.72480      0.07116
H      -0.32103     -0.86428      1.58817
H      -1.91086     -0.75242      0.81140
"""


if __name__ == "__main__":
    print("Testing mod_fragments in forward direction")
    mol = Molecule()
    mol.read_xyz(xyz_file)
    print(mol)
    print(f"Distance between H1 and H11 : {mol.distance(0,10)}")
    print(f"Centroid: {mol.centroid(mol.coordinates)}")
    print(f"RMSD : {mol.rmsd(mol.coordinates,mol.coordinates)}")
    print("centering Molecule")

    mol.center_coordinates()
    print(mol)
    print(f"Distance between H1 and H11 : {mol.distance(0,10)}")
    print(f"Centroid: {mol.centroid(mol.coordinates)}")

    mol2 = Molecule()
    mol2.read_xyz(xyz_file)
    print(f"RMSD : {mol.rmsd(mol2.coordinates,mol.coordinates)}")

    mol3 = Molecule()
    mol3.read_xyz(rotated_xyz)
    print(f"RMSD ROTATED : {mol.rmsd(mol2.coordinates,mol3.coordinates)}")

    mol3.coordinates = Molecule.kabsch(mol3.coordinates,mol2.coordinates)
    print(mol3)
    print(f"RMSD KABSCH : {mol.rmsd(mol2.coordinates,mol3.coordinates)}")

    print("Distance Matrix:")
    mol3.distance_matrix(broadcast=True)
    print(mol3.d_mat)
    mol3.distance_matrix(broadcast=False)
    print(mol3.d_mat)



