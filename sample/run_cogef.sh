#!/bin/bash

# setting the path to the cogef script to be used
cogef=$COGEF_BASEDIR/manual_cogef.py

# run $COGEF_BASEDIR/manual_cogef.py -h for a list of options.

######################
# scan h h
# use standard asymmetric method
$cogef -xyz opt.xyz  2 3  -cycles 150  -dx 0.05 -method "UB3LYP/6-31G*"  -nproc=4 -mem=3GB  -calcfc None
mkdir -p h_h             # make output dir for logfiles
mv cogef* h_h            # move log files to output dir
# generate final xyz trajectory file
for i in h_h/cogef_???.log ; do
          gxyz --input  $i
done > h_h.xyz


