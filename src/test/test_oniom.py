#!/usr/bin/env python3
""" Test all classes in molecule
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import numpy as np 
from cogef.molecule import Molecule
from cogef.molecule import OniomMolecule



template="""C-c3i--0.1171910   0     2.713000    -0.107000    -0.001000 H 
H-hci-0.0604880    0     2.781000     0.519000    -0.881000 H 
H-hci-0.0604880    0     2.781000     0.519000     0.879000 H 
H-hci-0.0604880    0     3.550000    -0.796000    -0.001000 H 
C-c3i--0.0852380   0     1.408000    -0.904000    -0.000000 H 
H-h1i-0.0584280    0     1.327000    -1.526000    -0.883000 H 
H-h1i-0.0584280    0     1.328000    -1.527000     0.882000 H 
S-s6i-0.9496850    0     0.000000     0.192000     0.000000 H 
O-oi--0.5707340    0    -0.001000     0.918000    -1.244000 H 
O-oi--0.5707340    0     0.001000     0.918000     1.244000 H 
C-c3i--0.0852380   0    -1.408000    -0.904000     0.002000 H 
H-h1i-0.0584280    0    -1.327000    -1.529000    -0.880000 H 
H-h1i-0.0584280    0    -1.328000    -1.525000     0.885000 H 
C-c3i--0.1171910   0    -2.713000    -0.107000    -0.001000 H 
H-hci-0.0604880    0    -2.780000     0.517000    -0.883000 H 
H-hci-0.0604880    0    -2.781000     0.521000     0.877000 H 
H-hci-0.0604880    0    -3.550000    -0.796000     0.000000 H 

     1      2 1.0      3 1.0      4 1.0      5 1.0 
     2      1 1.0 
     3      1 1.0 
     4      1 1.0 
     5      1 1.0      6 1.0      7 1.0      8 1.0 
     6      5 1.0 
     7      5 1.0 
     8      5 1.0      9 1.0     10 1.0     11 1.0 
     9      8 1.0 
    10      8 1.0 
    11      8 1.0     12 1.0     13 1.0     14 1.0 
    12     11 1.0 
    13     11 1.0 
    14     11 1.0     15 1.0     16 1.0     17 1.0 
    15     14 1.0 
    16     14 1.0 
    17     14 1.0 

NonBon 3 1 0 0 0.0 0.0 -2.0 0.0 0.0 -1.2 
VDW c3i     1.9080  0.1094
VDW hci     1.4870  0.0157
VDW h1i     1.3870  0.0157
VDW s6i     2.0000  0.2500
VDW oi     1.6612  0.2100
HrmStr1 c3i c3i  300.900    1.5375
HrmStr1 c3i s6i  233.500    1.8075
HrmStr1 s6i oi   512.700    1.4533
HrmStr1 c3i hci  330.600    1.0969
HrmStr1 c3i h1i  330.600    1.0969
HrmBnd1  c3i c3i s6i    62.07   110.220
HrmBnd1  c3i s6i oi     65.37   108.610
HrmBnd1  c3i s6i c3i    59.92   103.830
HrmBnd1  oi  s6i oi     73.59   120.050
HrmBnd1  c3i c3i h1i    46.39   109.560
HrmBnd1  hci c3i hci    39.40   107.580
HrmBnd1  hci c3i c3i    46.34   109.800
HrmBnd1  h1i c3i h1i    39.24   108.460
HrmBnd1  h1i c3i s6i    43.19   107.150
AmbTrs  c3i c3i s6i oi     0   0   0   0  0.0000000  0.0000000  0.1444444  0.0000000 1.0
AmbTrs  c3i c3i s6i c3i    0   0   0   0  0.0000000  0.0000000  0.1444444  0.0000000 1.0
AmbTrs  hci c3i c3i h1i    0   0   0   0  0.0000000  0.0000000  0.1555556  0.0000000 1.0
AmbTrs  hci c3i c3i s6i    0   0   0   0  0.0000000  0.0000000  0.1555556  0.0000000 1.0
AmbTrs  c3i s6i c3i h1i    0   0   0   0  0.0000000  0.0000000  0.1444444  0.0000000 1.0
AmbTrs  h1i c3i s6i oi     0   0   0   0  0.0000000  0.0000000  0.1444444  0.0000000 1.0
"""


if __name__ == "__main__":
    print("Testing OniomMolecule")
    mol = OniomMolecule()
    mol.read_oniom_template(template)
    print(mol)

