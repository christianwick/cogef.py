#!/usr/bin/env python3


import argparse
import logging
import sys

from cogef import cogef_logging
from cogef.molecule import Molecule
from cogef._version import __version__ 


def get_molecule_dict(trajectory):
    """
    input:  xyz trajectory = list Molecule.read_xyz_trj 
    
    returns: dict("point" = Molecule() )
    """
    molecules = dict()
    for mol in trajectory:
        data = Molecule()
        data.read_xyz(mol, process_comment=True)
        molecules[data.point] = data
    return(molecules)

def get_trajectories(xyzfiles):
    """ read in all trajectories and save them as dictionaries
     with the actual point as key. this way we can compare the
     energies per point (if it exists in the trajectory)
     we assume a comment line as follows:
     E(UB3LYP) = -118.996452168 | S**2 = 1.008 | point 069

     input:  xyzfiles : TextIOWrapper

     returns: trajectories : list
    """ 
    trajectories = []
    for fil in xyzfiles: 
        temp_trj = Molecule.read_xyz_trj(fil)
        trajectories.append(get_molecule_dict(temp_trj))
    return(trajectories)

def find_lowest_point(point,trajectories):
    """
    find the structure with the lowest energy at a defined point

    input:  point : int
            trajectories : list ( dict ( Molecule() ))
    
    returns: data_point_with_lowest_en = Molecule()
    """
    lowest_en = 0.0
    data_point_with_lowest_en = None
    for trj in trajectories:
        if point in trj.keys():
            if trj[point].energy <= lowest_en:
                data_point_with_lowest_en = trj[point]
                lowest_en = trj[point].energy
    return(data_point_with_lowest_en)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("xyzfiles",help="list of xyz files to compare", nargs="+", type=argparse.FileType("r"))
    parser.add_argument("-oxyz",help="output xyz trajextory", default=sys.stdout)
    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO"], default="INFO")
    
    args = parser.parse_args()

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file="min_xyz_traj.log",log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("min_xyz_traj")
    logger.info("Starting min_xyz_traj .... " )
    logger.info("Version: {}".format(__version__))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))

    # extract the molecules as list of trajectories from the
    # xyzfiles. The elements of the list are dictionaries (key = point)
    # and elements are of type Molecule()
    trajectories = get_trajectories(args.xyzfiles)

    # now we obtain a unique set of all the points in all
    # trajectories
    point_set = set()
    for trj in trajectories:
        point_set.update(trj.keys())
    
    # the final_xyz_trj will contain the structure with the 
    # lowest energy per point. we will sort the point set for this
    final_xyz_trj = []
    for point in sorted(point_set):
        data_point_with_lowest_en = find_lowest_point(point,trajectories)
        final_xyz_trj.append(data_point_with_lowest_en)

    # write the final trajectory
    for point in final_xyz_trj:
        point.write_xyz(args.oxyz,comment=point.comment,write_mode="a")          

