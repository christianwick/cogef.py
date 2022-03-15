#!/usr/bin/env python3
""" analyse.py - read xyz files with energy line as comment,
     compute distance between atoms
     plot energies
     write energy files and relative energy files for processing.
"""

import sys
import numpy as np
import csv
import logging
import argparse

from cogef import cogef_logging
from cogef import constants
from cogef.molecule import Molecule
from cogef._version import __version__ 


def get_molecules(trajectory):
    """
    input:  xyz trajectory = list Molecule.read_xyz_trj 
    
    returns: list( Molecule() )
    """
    molecules = []
    for mol in trajectory:
        data = Molecule()
        data.read_xyz(mol, process_comment=True)
        molecules.append(data)
    return(molecules)

def get_trajectories(xyzfiles):
    """ read in all trajectories and save them as a single list.
     we assume a comment line as follows:
     E(UB3LYP) = -118.996452168 | S**2 = 1.008 | point 069

     input:  xyzfiles : TextIOWrapper

     returns: trajectories : list 
    """ 
    trajectories = []
    for fil in xyzfiles: 
        temp_trj = Molecule.read_xyz_trj(fil)
        trajectories.extend(get_molecules(temp_trj))
    return(trajectories)

def get_dist_tuple(raw):
    """ 
    input raw = str

    returns dist_tuple = list ( tuple )
    """
    input_tuple = []
    temp = raw.split(";")
    for pair in temp:
        #input_tuple.append((int(pair.split()[0])-1,int(pair.split()[1])-1))
        input_tuple.append(tuple( [int(x) - 1  for x in pair.split() ] ))
    return(input_tuple)



class analyse_cogef():
    def __init__(self):
        self.energies = []
        self.spin = []
        self.distances = [] # this is a list of lists!
        self.distances_label = []
        self.rel_distances = []
        self.rel_distances_label = []
        self.rel_en_kj_mol = []
        self.molecules = []
        self.num_struc = []
        self.v_eff = []
        self.v_eff_label = []
        self.dist_tupl = []
        self.forces = []
        self.max_force = None
        self.eps = [] # difference of distances
        self.eps_label = [] # difference of distances

    def fill_with_molecules(self,trajectories, filter_en = False,
                            use_first_minimum = False, use_first_point = False, ref_en=None):
        for nn, mol in enumerate(trajectories):
            self.energies.append(mol.energy)
            self.spin.append(mol.spin)
            self.molecules.append(mol)
            self.num_struc.append(nn+1)
        self.energies = np.array(self.energies)
        self.rel_en_kj_mol = self._compute_rel_energy(self.energies,use_first_point=use_first_point,
            use_first_minimum=use_first_minimum, reference_energy=ref_en)
        if isinstance(filter_en, float):
            self._filter_data_by_energy(filter_en)
        en_max = self._find_maxima(self.energies)
        for val in en_max:
            logger.info("Found Energy maximum at point {:3} = {:7.3f} kJ/mol".format(val, self.rel_en_kj_mol[val]))

    def _compute_rel_energy(self,energies,conv_factor=constants.hartree_to_kJ_mol, use_first_minimum = False, 
                            use_first_point = False, reference_energy=None):
        if use_first_minimum :
            en_min = self._find_first_energy_minimum(energies)
        elif use_first_point:
            en_min = energies[0]
        elif reference_energy:
            en_min = reference_energy
        else:
            en_min = np.min(energies)
        return ( (energies - en_min ) * conv_factor )

    def compute_veff(self,forces=[0.5,1.0,2.0,3.0,4.0,5.0],use_first_minimum=False):
        """ compute the effective potential according to Brügner et al. as
            U_F(d) = U(d) - Fd
        
        input:  energies = np.array [Hartree]
                distance = np.array [A]
                force = float    [nN]

        """
        self.v_eff = []
        self.v_eff_label = []
        distances = self.rel_distances[:,0] * constants.angstrom_to_m   # convert from angstrom to m

        for force in forces: 
            self.v_eff_label.append(f"Veff {force:3.2f} [kJ mol-1]")
            work = force * distances * constants.nNm_to_kJ_mol  # in kJ / mol
            self.v_eff.append ( self._compute_rel_energy(self.rel_en_kj_mol - work,conv_factor=1.0,use_first_minimum=use_first_minimum) )
        for nn, veff in enumerate(self.v_eff):
            en_max = self._find_maxima(veff)
            for val in en_max:
                logger.info("Found maximum of {:3} point {:3} = {:7.3f} kJ/mol".format(self.v_eff_label[nn], val, veff[val]))
        self.v_eff = np.array(self.v_eff).T # we transpose the matrix, to facilitate printing later on.


    @staticmethod
    def _compute_distances(molecules, distance_tuples):
        """ 

        input: distance_tuples = list( tuple )

        return: distances = list ( list )
        """
        distances = []
        distances_label = []
        dist_tupl = distance_tuples
        for i,j in dist_tupl:
            distances_label.append(f"dist({i+1} - {j+1}) [Ang]")
        for mol in molecules:
            temp_dist = []
            for i,j in dist_tupl:
                temp_dist.append(mol.distance(i,j))
            distances.append(temp_dist)
        distances = np.array(distances)
        return(distances,distances_label)
    
    def save_distances(self,distance_tuples):
        self.distances = []
        self.rel_distances = []
        self.distances_label = []
        self.dist_tupl = distance_tuples
        self.distances, self.distances_label = self._compute_distances(self.molecules, self.dist_tupl)
        self._compute_rel_distance_min_energy_structure()
        self._compute_strain()

    def save_eps(self,eps_tuples):
        """ compute eps = difference of distances """
        temp_eps = []
        self.eps_label=[]
        # eps tuples are of form
        # [(1,2 ,3 ,4),(....),...]
        for tup in eps_tuples:
            temp_dist_tuples = [(tup[0],tup[1]),(tup[2],tup[3])]
            dist, label = self._compute_distances(self.molecules,temp_dist_tuples)
            temp_eps.append(np.squeeze(np.diff(dist)))
            self.eps_label.append("eps({}-{})".format(label[0],label[1]))
        # transpose array to get right shape (len(molecules),len(eps tuples))
        self.eps = np.array(temp_eps).T

    def _compute_rel_distance_min_energy_structure(self):
        en_minimum = np.argmin(self.energies)
        self.rel_distances  = self.distances - self.distances[en_minimum]
        self.rel_distances_label = ["rel " + str(x) for x in self.distances_label]

    def _compute_strain(self):
        """ engineering strain
        e = ( l - L ) / L 
            e : engineering strain, l current distance, L initial distance
        """
        self.strain = []
        l_zero = self.distances[0,:]
        self.strain = (self.distances - l_zero ) / l_zero
        self.strain_label = ["strain " + str(x)[:-6] for x in self.distances_label]
        
    def _filter_data_by_energy(self,cut_off_en=400.):
        """ filter the data and remove high energy points
            this function does not remove computed distances or V_eff data.
            those computed data must be recomputed after filtering!

        input:  cut_off_en = float # kj mol-1

        """
        for nn in range( len( self.rel_en_kj_mol) -1 , -1, -1):
            if self.rel_en_kj_mol[nn] >= cut_off_en:
                self.energies = np.delete(self.energies,[nn])
                self.rel_en_kj_mol = np.delete(self.rel_en_kj_mol ,[nn])
                self.spin.pop(nn)
                self.molecules.pop(nn)
                self.num_struc.pop(nn)

    @staticmethod
    def _find_maxima(energies):
        last = energies[0]
        rising = False
        maxima = []
        for nn, val in enumerate(energies[1:]):
            if val > last:
                rising = True
            elif val < last and rising:
                maxima.append(nn)
                rising = False
            last = val
        return(maxima)



    @staticmethod
    def _find_first_energy_minimum(rel_energies):
        min_en = rel_energies[0]
        for val in rel_energies[1:]:
            if val >= min_en:
                break
            min_en = val
        return(min_en)

    def compute_forces(self, dist_index=0):
        de = ( self.energies[1:] - self.energies[:-1] ) * constants.hartree_to_J
        dx = ( self.distances[1:,dist_index] - self.distances[:-1,dist_index] ) * constants.angstrom_to_m
        self.forces = np.append( [0.0], np.divide(de,dx) * 1.0E9 ) # from N to nN
        self.max_force = np.max(self.forces)
        logger.info("Max Force = {:7.3f} nN".format(self.max_force))

    def write_force_csv(self,of, print_strain=True, dist_index=0):
        header = ["Num Struc"]
        # set distances
        header.append(self.distances_label[dist_index])
        header.append(self.rel_distances_label[dist_index])
        if print_strain:
            header.append(self.strain_label[dist_index])
        header.append("Force nN")
        # start writing rows
        writer = csv.writer(of)
        writer.writerow(header)
        for nn in range(len(self.forces)):
            row = [self.num_struc[nn]]
            row.append(self.distances[nn,dist_index])
            row.append(self.rel_distances[nn,dist_index])
            if print_strain:
                row.append(self.strain[nn,dist_index])
            row.append(self.forces[nn])
            writer.writerow(row)
        

    def write_en_csv(self,of,print_veff=False,print_strain=False,print_spin=True):
        header = ["Num Struc"]
        # set distances
        if len(self.distances) == len(self.energies):
            header.extend(self.distances_label)
            header.extend(self.rel_distances_label)
            if print_strain:
                header.extend(self.strain_label)
        if len(self.eps) == len(self.energies):
            header.extend(self.eps_label)
        header.append("Energy [au]")
        header.append("Rel Energy [kJ mol-1]")
        if print_veff == True:
            header.extend(self.v_eff_label)
        if print_spin == True:
            header.append("<S**2>")
        writer = csv.writer(of)
        writer.writerow(header)
        for nn in range(len(self.energies)):
            row = [self.num_struc[nn]]
            if len(self.distances) == len(self.energies):
                row.extend(self.distances[nn])
                row.extend(self.rel_distances[nn])
                if print_strain:
                    row.extend(self.strain[nn])
            if len(self.eps) == len(self.energies):
                row.extend(self.eps[nn])
            row.append(self.energies[nn])
            row.append(self.rel_en_kj_mol[nn])
            if print_veff == True:
                row.extend(self.v_eff[nn])
            if print_spin == True:
                row.append(self.spin[nn])
            writer.writerow(row)




if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("xyzfiles",help="list of xyz files to analyse", nargs="+", type=argparse.FileType("r"))
    parser.add_argument("-csv",help="output csv file", type=argparse.FileType("w"), default="summary.csv")
    parser.add_argument("-fcsv",help="output force csv file", type=argparse.FileType("w"), default="summary_f.csv")
    parser.add_argument("-fil",help="filter out energies higher than nn", type=float, default=None)
    parser.add_argument("-d",help="distances to analyse as pairs, e.g. '0 1; 2 3'. the first distance pair"+\
                                  "will be used to compute veff.", default=None)
    parser.add_argument("-eps",help="compute difference of distances. needs 2 atom pairs "+\
                                "'0 1  2 3 = (d0-1 - d2-3)'.", default=None)                              
    parser.add_argument("-strain",help="compute strain", action="store_true", default=False)
    parser.add_argument("-compforce",help="compute forces", action="store_true", default=False)
    parser.add_argument("-veff",help="compute veff", action="store_true", default=False)
    parser.add_argument("-spin",help="print <S**2>", action="store_true", default=False)
    parser.add_argument("-movie",help="sample a very large number of forces for movies.", action="store_true", default=False)
    parser.add_argument("-forces",help="forces to use in veff calculations, e.g. '0.5; 1.0; 1.5'", default=False)
    parser.add_argument("-first_min",help="use first minimum instead of global minimum to compute veff", action="store_true", default=False)
    parser.add_argument("-first_point",help="use first point instead of global minimum to compute relative energies", action="store_true", default=False)
    parser.add_argument("-ref_en",help="use reference energy given in Hartree", type=float, default=None)


    group_logging = parser.add_argument_group("logging")
    group_logging.add_argument("-log_level", help="set the log level", choices=["DEBUG","INFO"], default="INFO")
    group_logging.add_argument("-stream_level", help="set the streaming level", choices=["DEBUG","INFO","CRITICAL"], default="CRITICAL")
    
    args = parser.parse_args()

    # START LOGGING HERE
    # this way no logs are created for e.g. -h or --version 
    cogef_logging.logging_init(log_file="analyse.log",log_level = args.log_level, stream_level = args.stream_level)
    logger = logging.getLogger("analyse")
    logger.info("Starting analyse.py .... " )
    logger.info("Version: {}".format(__version__))
    for fil in args.xyzfiles:
        logger.info("XYZfiles: {}".format(fil.name))
    logger.info("Command line Arguments: ")
    for arg,value in (vars(args)).items():
        logger.info("    -{:10s} : {}".format(arg,value))


    # extract the molecules as list of trajectories from the
    # xyzfiles. The elements are of type Molecule()
    trajectories = get_trajectories(args.xyzfiles)

    # convert our trajectory to a dictionary with lists as elements.
    data = analyse_cogef()
    data.fill_with_molecules(trajectories,filter_en=args.fil,use_first_minimum = args.first_min, 
                            use_first_point = args.first_point,ref_en=args.ref_en)

    # compute distances
    if args.d:
        dist_tuples = get_dist_tuple(args.d)
        data.save_distances(dist_tuples)

        if args.veff:
            if args.movie:
                data.compute_veff(forces=np.arange(0.01,5.01,0.01),use_first_minimum=args.first_min)
            elif args.forces:
                forces = [ float(x) for x in args.forces.split(";") ]
                data.compute_veff(forces=forces,use_first_minimum=args.first_min)
            else:
                data.compute_veff(use_first_minimum=args.first_min)
        if args.compforce:
            data.compute_forces()
            data.write_force_csv(args.fcsv)
    # eps is internally used for differences of distances. we will compute the difference between 
    # two pairs of atoms a b c d -> d(a-b) - d(c-d)
    if args.eps:
        eps_tuples = get_dist_tuple(args.eps)
        data.save_eps(eps_tuples)
    data.write_en_csv(args.csv,print_veff=args.veff,print_strain=args.strain)
    logger.info("Finished analyse.py.")
    logger.info("")
    logger.info("")
    


