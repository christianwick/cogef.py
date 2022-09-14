#!/bin/bash


$COGEF_BASEDIR/run_cogef.py -xyz large.xyz  -at "2,31" -cycles 5 -dx 0.05 \
	-method "ONIOM(RB3LYP/6-31G*:AMBER=SOFTONLY)" -quadmac \
	-cm "0 1 0 1 0 1" -oniom large.template -oniomopt me -constraint opt_flag \
       	-nproc=12 -mem=30GB -trajectory traj_cogef.xyz 
mkdir -p data/
mv cogef* data/

