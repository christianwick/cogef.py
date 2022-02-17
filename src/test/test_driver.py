#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""
# Modify path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import cogef.gaussian as g 
import cogef.driver as driver 
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
    print("Testing Class cogef_dirver:")
    print("----------------------------")
    print()
    print()
    cogef = driver.cogef_loop(atom1=0,atom2=10,xyz=xyz_file,cycles=2,dx=5.0,dp=1.0,fragment=[0,1,2,3,4,5,6],runtype="restricted")
    cogef.run()

