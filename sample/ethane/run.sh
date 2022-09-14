#!/bin/bash

# run cogef on opt.xyz stretching atoms 2 3 by 0.05 angstroms
$COGEF_BASEDIR/run_cogef.py -xyz opt.xyz  -at "3,7"  -cycles 10 -dx 0.05 -method "UB3LYP/6-31G*" -runtype unrestricted -trajectory ethane_unrestricted.xyz

mkdir -p data             # make output dir for logfiles
mv cogef* data            # move log files to output dir




