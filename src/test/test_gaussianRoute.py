#!/usr/bin/env python3
""" fragement_cogef.py - read xyz files and write gaussian input with modifed coords
 
"""
# Modify path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))

import cogef.gaussian as g 
import sys 



if __name__ == "__main__":
    print()
    print("Testing Class GaussianInput:")
    print("----------------------------")
    print()
    route = g.classGaussianRoute()
    route.route = ["#P","B3lyp/6-31G*","nosymm","stable=opt"]
    route.geom = ["connect"]
    route.opt = ["maxcyc=50"]
    route.scf = ["QC"]
    route.guess = ["mix"]
    print(str(route))
    route = g.classGaussianRoute()
    print(str(route))
