#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import cogef.gaussian as g 
import sys 


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
    print()
    print("Testing Class GaussianInput:")
    print("----------------------------")
    print()
    test = g.GaussianInput()
    print()
    print("_write_link0:") 
    test._write_link0(of=sys.stdout)
    test._write_link0(of=sys.stdout,args={"%rwf": "test.rwg"})
    print()
    print("_write_route:")
    test._write_route(of=sys.stdout)
    test._write_route(of=sys.stdout, args={"opt": ["maxcyc=100","calcfc","rfo"] })
    print()
    print("_write_title:")
    test._write_title(of=sys.stdout)
    print("\n\n")
    print("---------------------")
    print("\n\n")
    print("Read write molecule")
    test.molecule.read_xyz(xyz_file)
    print(test.molecule)
    print("\n\n")
    print("---------------------")
    print("\n\n")
    print("Test write_input():")
    test.write_input(of=sys.stdout,modredundant=["1 11 F","2 8 F"])
    print("\n\n")
    print("---------------------")
    print("\n\n")
    print("Test write_sp_input():")
    test.write_sp_input(of=sys.stdout)
