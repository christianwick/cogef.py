#!/usr/bin/env python3
""" Test all classes in mulliken.py
 
"""
# Modif path so that cogef is visible
import sys, os
sys.path.append(os.path.abspath(os.path.join('..')))


from cogef.mulliken import Mulliken

raw_str = """     2  C    0.001520   0.012112
     5  C   -0.011722   0.963967
     8  C    0.025860  -0.776637
    11  N   -0.171410  -0.179924
    13  C    0.155753  -0.019517 """
gaussian_raw_data = [ x for x in raw_str.split("\n") ]

if __name__ == "__main__":
    mull = Mulliken()
    mull._read_gaussian_raw_data(gaussian_raw_data)
    print(mull)
    mull._read_gaussian_raw_data(gaussian_raw_data)
    print(mull)

    mull.sort_by("spin")
    print(mull)
    mull.sort_by("atom_number")
    print(mull)
    mull.sort_by("atom_name")
    print(mull)
    mull.sort_by("charge")
    print(mull)

    print(mull.find_spin_larger_than())