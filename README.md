# COGEF.PY - Tools to run COnstraint Geometries simulate External Force calculations

## Usage


in bash: add the following lines to your .bashrc
```bash
# COGEF.py
export COGEF_BASEDIR=your-path-to-cogef.py/src
```

you can access the cogef tools in your scripts using:
```bash
$COGEF_BASEDIR/manual_cogef.py
```
see example run_script in ./sample


Note: COGEF.py works best in cooperation with the chem_tools provided by 
christian.wick@fau.de


## Scripts

src/ 

manual_cogef.py : run manual cogef calculations

targeted_cogef.py : run targeted cogef calculations

analyse_cogef.py : script to analyse cogef xyz trajectories 
                   produced with gxyz.py

simple_orca_sp.py : script to run single points on xyz trajectories 
                   with orca (needs additional configurations, details follow)

gausp.py : script to run single points on xyz trajectories
           with gaussian (needs additional configuration, details follow)

sample/ 

cogef_rung16  : script to run gaussian. Needs to be in your PATH!

run_cogef.sh : example run script for cogef calculations 




