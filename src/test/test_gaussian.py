#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""
# Modify path so that cogef is visible
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


mulliken = """
 Mulliken charges and spin densities with hydrogens summed into heavy atoms:
               1          2
     2  C    0.001520   0.012112
     5  C   -0.011722   0.963967
     8  C    0.025860  -0.776637
    11  N   -0.171410  -0.179924
    13  C    0.155753  -0.019517
 Electronic spatial extent (au):  <R**2>=           1543.5704
 """

if __name__ == "__main__":
    print()
    print("Testing Class GaussianInput:")
    print("----------------------------")
    print()
    test = g.GaussianInput()
    print()
    test.molecule.read_xyz(xyz_file)
    test.link0.chk="tst.chk"
    route="#P PM6 IOP(2/15=3) OPT=Modredundant"
    test.write_inputfile(link1=False,route=route,geom=True, modredundant=["1 2 F","1 3 F"])
    route="#P PM6 IOP(2/15=3) Freq guess=read geom=allcheck"
    test.write_inputfile(link1=True,route=route,geom=False)
    route="#P PM6 IOP(2/15=3) OPT guess=read geom=allcheck"
    test.write_inputfile(link1=True,route=route,geom=False )
    route="#P PM6 IOP(2/15=3) FREQ guess=read geom=allcheck"
    test.write_inputfile(link1=True,route=route,geom=False )

    # CheckGaussianLogfile
    mylog = g.CheckGaussianLogfile("example.log")
    mylog.read_log()
    print(f"Reading File {mylog.filename}")

    print(f"Instability: {mylog.instability}")
    print(f"Error termination: {mylog.error}")
    print(f"Stationary Point: {mylog.stationary}")
    print(f"Structure: {mylog.molecule}")

    print("Mulliken Data")
    print(mylog.mulliken)

