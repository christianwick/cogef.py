# COGEF.PY - Tools to run COnstraint Geometries simulate External Force calculations


## Description
COGEF.py is a set of python scripts that interface to external DFT codes to produce and 
analyse COGEF calculations.

Currently supported DFT codes: Gaussian 16 (other codes might be added in the future.

The functionality includes:
  - Standard COGEF calculations with restricted and unrestricted methodologies
  - AMBER COGEF calculations utilising the AMBER force field (implemented in Gaussian)
  - ONIOM COGEF calculations utilising the ONIOM formalism (as implemented in Gaussian)

The following tools are included in cogef.py:
  - src/run_cogef.py : actual cogef driver to run cogef calculations
  - src/analyse.py   : analyse cogef xyz trajectory files produced by run_cogef.py
  - src/gausp.py : conduct single point calculations on a cogef xyz trajectory file and 
               produce a new cogef trajectory file which can be analysed be analyse.py
 
helper tools which are considered to be beta and subject to drastic changes:
  - src/simple_orca_sp.py : conduct single point calculations on a cogef xyz trajectory file with orca
  - src/read_log.py : parse a series of gaussian log files and produce a cogef xyz trajectory file
  - src/min_xyz_traj.py : read two cogef xyz trajectory files and produce a new trajectory file containing
                      the structure with the lowest electronic energy at each point.
use these tools with care.


## Installation

currently only installation from source is possible. 

### install from source 

1. clone this git repository

2. run the installer script
```bash
# this will install pip and all dependecies into your current python enviroment
# it is highly recommended to use a virtual python environment (pyenv, conda, ..) 
cd cogef.py
./install.sh 
```
3. adjust the sample/cogef_rung16 file to your needs
sample/cogef_rung16: 
```bash
#!/bin/bash
###                                           ###
###           g16 run script                  ###
###         modfied for cogef                 ###
# --------------------------------------------- #
#
# Problems: Christian.wick@fau.de
# Version: 1.0.0
# --------------------------------------------- #

# run_cogef.py will call cogef_rung16 (this file) with the gaussian.com file
# as first argument. it expects a gaussian.log file as output.
# make sure that run_cogef.py is visible in your path when executing run_cogef.sh.

# adjust these lines to your needs (e.g. you could also call another script 
# that takes care about the gaussian calculation - helpful on some computing clusters etc.)

# RunGauss
logfile="${1%.com}.log"
$g16root/g16/g16 < $1 > $logfile 

```
 - run_cogef.py and gausp.py will assume that an executable "cogef_rung16" is in your path, which 
   is used to execute the gaussian job. Currently, this file needs to be adjusted to your local 
   gaussian installation. 

 - on ubunutu 21.04 using bash you could just copy your modified version to ~/bin and add this directory
   to your PATH in your .bashrc file, e.g.
  ```bash
  mkdir -p ~/bin
  cp sample/cogef_rung16 ~/bin
  ```
  and add to your .bashrc file the following line:
  ```bash
  PATH=$PATH:~/bin
  ```

## Usage

in bash: add the following lines to your .bashrc to access the cogef tools in your 
user session. 
```bash
# COGEF.py
export COGEF_BASEDIR=your-path-to-cogef.py/src
```
Now you can see all available options like this:
```bash
# run_cogef.py 
$COGEF_BASEDIR/run_cogef.py -h
# or gausp.py
$COGEF_BASEDOR/gasup.py -h
```

you can access the cogef tools in your scripts using
```bash
$COGEF_BASEDIR/run_cogef.py
```


## Example 
we provide two example inputs for standard cogef and 
oniom type cogef calculations in the ./sample directory

### Ethane 
In this example we will perform a COGEF calculation with the unrestricted procedure as described in [1]. 
The initial optimised coordinates for ethane can be found in the opt.xyz script.
We assume that you followed the install instructions to set up cogef.py. You should have adjusted the
file sample/cogef_rung16 and added it to your PATH. We also assume that you did set the COGEF_BASEDIR variable
correctly.

```bash
# copy the contents of the sample directory to a new directory outside the cogef
# base directory. e.g. your home directory 
cp -r $COGEF_BASEDIR/../sample/ ~/ 
cd ~/sample/ethane
```
In this directory you will find two files, opt.xyz and run.sh
opt.xyz:
```text
8
E(RB3LYP) = -79.8304198606 ETHANE B3LYP/6-31G*
H       1.021431     1.162936     0.000000 
C       0.000084     0.765388     0.000000 
H      -0.510417     1.163580     0.884397 
H      -0.510417     1.163580    -0.884397 
C       0.000084    -0.765157    -0.000000 
H       0.509724    -1.163738     0.884663 
H       0.509724    -1.163738    -0.884663 
H      -1.021053    -1.164002    -0.000000 
```
run.sh:
```text
#!/bin/bash

# run cogef on opt.xyz stretching atoms 2 3 by 0.05 angstroms
$COGEF_BASEDIR/run_cogef.py -xyz opt.xyz  -at "3,7"  -cycles 10 -dx 0.05 -method "UB3LYP/6-31G*" -runtype unrestricted -trajectory ethane_unrestricted.xyz

mkdir -p data             # make output dir for logfiles
mv cogef* data            # move log files to output dir
```

opt.xyz contains the optimised coordinates (optained previously by a geometry optimisation 
in gaussian 16) of an ethane molecule.
run.sh is a bashscript that will conduct the cogef caclulations with the following options:
```text
-xyz opt.xyz : starting coordinates in file opt.xyz (in xyz file format)
-at "3,7" : increasing the distance between atoms 3 and 7 (counting starts at 1)
-cycles 10 : increase the distance between the two atoms for 10 cycles
-dx 0.05 : increment in angstrom
-method "UB3LYP/6-31G*" : use this functional (important: If you want to conduct unrestricted calculations, you need to specify UHF/UKS
 formalism in gaussian (**U**B3LYP)
-runtype unrestricted : use the unrestricted procedure as described in [1] 
-trajectory ethane_unrestricted.xyz : write the cogef trajectory xyz file to ethane_unrestricted.xyz
```

You can start the calculation by executing the run.sh script
```bash
./run.sh
```

This will produce a lot of logging output. this is also saved to job_cogef.log for your convenience.
```text
~/sample/ethane $ ll
total 1136
drwxr-xr-x 3 christian smith    4096 Sep 14 12:00 ./
drwxr-xr-x 3 christian smith    4096 Sep 14 11:31 ../
-rw-r--r-- 1 christian smith     541 Sep 14 12:00 checkpoint.xyz
-rw-r--r-- 1 christian smith    5410 Sep 14 12:00 checkpoint.xyz_trj.xyz
drwxr-xr-x 2 christian smith    4096 Sep 14 12:00 data/
-rw-r--r-- 1 christian smith    5410 Sep 14 12:00 ethane_unrestricted.xyz
-rw-r--r-- 1 christian smith 1105920 Sep 14 12:00 guess.chk
-rw-r--r-- 1 christian smith   12720 Sep 14 12:00 job_cogef.log
-rw-r--r-- 1 christian smith     401 Sep 14 11:32 opt.xyz
-rwxr-xr-x 1 christian smith     349 Sep 14 11:49 run.sh*
```
The individual gaussian log files are stored in the directory ./data/ if the job did not get interrupted. 
The files checkpoint.xyz and checkpoint.xyz_trj.xyz contain the **unoptimized** input structures for each cycle.
The optimised final cogef trajectory is stored in ethane_unrestricted.xyz.

The comment line for individial data points looks like this:
```text
E(UB3LYP) = -79.8113126590 | S**2 = 0.000 | point 009 
```
a csv file containing the energies and distances for plotting can be obtained by
```bash
$COGEF_BASEDIR/analyse.py ethane_unrestricted.xyz  -d "3 7"
```

### ONIOM / AMBER
In the ./sample/oniom directory you can find an example to run cogef calculations with ONIOM support. 
the example run file is a bash script called ./run.sh. In addition to the large.xyz input file, the 
large.template file is necessary.
```bash
./run.sh
```
This section will be expanded in future releases.

## License / Permission of use

The code is provided as-is without any warranty for its correctness under the MIT license (see `LICENSE` file)

