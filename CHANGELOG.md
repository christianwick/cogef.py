# Changelog

## 1.5.1 Cogef.py 21.10.21
### added
- added the analysis tool analysis.py
- analysis.py will replace the old tool analyse_cogef.py in future versions

## 1.4.1 Cogef.py 11.09.21
### fixed
-simple_orca_sp.py: reverse direction works now as expected

## 1.4.0 Cogef.py 10.09.21
### added
#### simple_orca_sp.py: 
  - start_at a certain point and go in FORWARD or REVERSE direction
#### min_xyz_traj.py
  - kabsch fitting of final trajectories
#### molecule.py
  - functions to compute the distance matrix
  - functions to do KABSCH fitting
  - functions to compute RMSD between structures

## 1.3.0 Cogef.py 07.09.21
### added
- mulliken spin density is parsed and used to monitor which bond breaks
- summary of spin density and distance between the atoms with largest spin density
  is printed to the log file
- log file name can be changed in targeted_cogef
- targeted_cogef can read additional modredundant sections from the input

### deprecated
- manual_cogef is no longer maintained. use target cogef without fragment instead for
  identical results
- fragment_cogef is removed from the repository and no longer needed.

### changed
- modstruct: mod_single_atom will now move atom1 (if symm = False)

## 1.2.0 Cogef.py 03.09.21
### added
- looging in gausp.py
- logging in simple_orca_sp.py
- updated manual_cogef to use same modules as targeted_cogef.py
- added spin checking to manual_cogef

## 1.1.0 Cogef.py 27.08.21
### added
- parsing of gaussian log files in gaussian module
- reverse keyword in targeted cogef
- xyz data is extracted from gaussian log files
- min_xyz_traj.py to analyse reverse and forward cogef trajectories
### missing features:
- manual cogef needs to be updated to use the new gaussian module.

## 1.0.2 Cogef.py 19.08.21

### Fixed
- error terminations are now detected in manual_cogef and targeted_cogef. 

## 1.0.1 Cogef.py 11.08.21

### Fixed
- correct version is shown now

## 1.0.0 Cogef.py 10.08.21

### Added
- logging functionality
- grouping of command line arguments
- version control

  
