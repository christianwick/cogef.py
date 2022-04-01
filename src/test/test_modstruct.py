#!/usr/bin/env python3
""" Test all classes in modstruct
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import numpy as np 
from cogef.gaussian import GaussianInput
from cogef.modstruct import mod_fragments, mod_single_atom

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

xyz_file_separated="""11
separated
H      -4.67902300     0.36643337     0.00000000 
C      -3.77965600    -0.26127563     0.00000000 
H      -3.82176800    -0.90824963     0.88532100 
H      -3.82176800    -0.90824963    -0.88532100 
C       2.50004800     0.58978363     0.00000000 
H       2.49994300     1.24850663     0.87891800 
H       2.49994300     1.24850663    -0.87891800 
C       3.77962700    -0.26131937     0.00000000 
H       3.82175800    -0.90830337     0.88530200 
H       3.82175800    -0.90830337    -0.88530200 
H       4.67904200     0.36634963     0.00000000 
"""


if __name__ == "__main__":
    coords = GaussianInput()
    # test mod_fragments in forward direction
    print("Testing mod_fragments in forward direction")
    coords.molecule.read_xyz(xyz_file)
    counter = 0
    while counter < 2:
        coords.molecule.coordinates = mod_fragments(coords = coords.molecule.coordinates,
                atom1 = 0, atom2 = 10,
                frag1=[0,1,2,3],
                frag2=[4,5,6,7,8,9,10],  
                exclude= [5,6],
                dx = 1.0, 
                alpha = 0.5, 
                symmetric=True)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1

    # test mod_fragments in reverse direction
    print("Testing mod_fragments in reverse direction")
    coords.molecule.read_xyz(xyz_file_separated)
    counter = 0
    while counter < 2:
        coords.molecule.coordinates = mod_fragments(coords = coords.molecule.coordinates,
                atom1 = 0, atom2 = 10,
                frag1=[0,1,2,3],
                frag2=[4,5,6,7,8,9,10], 
                exclude = [5,6],
                dx = -0.5, 
                alpha = 0.5, 
                symmetric=True)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1


    print("Test mod_single_atom")
    coords.molecule.read_xyz(xyz_file)
    counter = 0
    while counter < 2:
        coords.molecule.coordinates = mod_single_atom(coords = coords.molecule.coordinates,
                atom1 = 0, atom2 = 10,
                dx = 0.5, 
                symmetric=True)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1

    print("Testing mod_fragments for single atom movement")
    coords.molecule.read_xyz(xyz_file)
    counter = 0
    while counter < 2:
        coords.molecule.coordinates = mod_fragments(coords = coords.molecule.coordinates,
                atom1 = 0, atom2 = 10,
                frag1=[], 
                frag2=[],
                dx = 5.0, 
                alpha = 1.0, 
                symmetric=False)
        print("Cycle: {}".format(counter))
        print(str(coords.molecule))
        counter += 1

    # test mod_fragments in forward direction with strain
    print("Testing mod_fragments in forward direction")
    coords.molecule.read_xyz(xyz_file)
    counter = 0
    dx = 1.0
    while counter < 10:

        current_strain = counter * dx / 4.35807
        coords.molecule.coordinates = mod_fragments(coords = coords.molecule.coordinates,
                atom1 = 0, atom2 = 10,
                frag1=[0,1,2,3],
                frag2=[4,5,6,7,8,9,10],  
                exclude= [5,6],
                dx = dx, 
                alpha = 0.0,
                current_strain=current_strain, 
                symmetric=True)
        print("Cycle: {} Strain: {}".format(counter,current_strain))
        print(str(coords.molecule))
        counter += 1