import numpy as np
import sys

import logging

from io import TextIOWrapper

from cogef.constants import PSE

logger = logging.getLogger("molecule")

class Molecule():
    """ the Molecule base class """
    def __init__(self, elements=[], 
                    coordinates=[],
                    fragments=[],
                    charge=[0],
                    multiplicity=[1]):
        logger.debug("Initialise new molecule")
        self.elements = elements
        self.fragments = fragments
        self.coordinates = np.array(coordinates)
        self.charge = charge
        self.multiplicity = multiplicity

    def __str__(self):
        string=""
        for i,j in zip(self.charge,self.multiplicity):
            string += "{} {} ".format(i,j)
        string += "\n"
        for el, coord in zip(self.elements,self.coordinates):
            string += "{:3s} {:14.8f} {:14.8f} {:14.8f} \n".format(el,*coord)
        return(string)

    def read_xyz(self,inpstr=None):
        logger.debug("Reading coordinates...")
        if isinstance(inpstr,str):
            lines = inpstr.split("\n")
        elif isinstance(inpstr,TextIOWrapper):
            lines = inpstr.readlines()
        else:
            lines=0
        num_atoms = int(lines[0])
        coords=[]
        elements=[]
        for ii in range(num_atoms):
            temp = lines[ii+2].split()
            coords.append([float(temp[1]), float(temp[2]), float(temp[3]) ])
            elements.append(str(temp[0]))
        self.elements = elements
        self.coordinates = np.array(coords)
        logger.debug("Finished reading coordinates.")
    
    def write_xyz(self,filename,comment="", write_mode="w"):
        logger.debug(f"Writing xyz coordinates to file {str(filename)}...")
        if isinstance(filename,str):
            of = open(filename, write_mode)
            need_to_close = True
        else:
            of = filename 
            need_to_close = False
        num_atoms=len(self.coordinates)
        of.write("{:} \n".format(num_atoms))
        of.write("{:80s} \n".format(comment.strip()[:80]))
        for ii in range(num_atoms):
            of.write("{:10s} {:14.8f} {:14.8f} {:14.8f} \n".format(
                str(self.elements[ii]), *self.coordinates[ii]))
        if need_to_close:
            of.close()
        logger.debug("Finished writing xyz coordinates...")
    
    def read_gaussian_raw_coords(self,raw_coords):
        """
        convert raw coordinate string from gaussian log files to the molcule format
        since gaussian uses atom numbers instead of element names we convert them 
        with our implemented PSE.

        input:  raw_coords = list
        """
        elements = []
        coords = []
        for line in raw_coords:
            temp = line.split()
            elements.append(PSE[int(temp[1])])
            coords.append([float(x) for x in temp[3:]])
        self.elements = elements
        self.coordinates = np.array(coords)

    @staticmethod
    def read_xyz_trj(inpstr=None):
        logger.debug("Reading xyz trajectory ...")
        if isinstance(inpstr,str):
            lines = inpstr.split("\n")
        elif isinstance(inpstr,TextIOWrapper):
            lines = inpstr.readlines()
        else:
            lines=0
        molecules = []
        linecount = 0 
        while linecount < len(lines):
            num_atoms = int(lines[linecount])
            molecules.append( "".join(lines[linecount:linecount + num_atoms + 2]))
            linecount += num_atoms + 2   
        logger.debug("Finished reading xyz trajectory ...")        
        return(molecules)
        


class OniomMolecule(Molecule):
    def write_oniom(self):
        pass

    def write_connectivity(self):
        pass

    def write_parm(self):
        pass
        

              

 