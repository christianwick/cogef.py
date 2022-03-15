""" molecule.py contains all Molecule classes """

import numpy as np
import sys

import logging

from io import TextIOWrapper

from cogef.constants import PSE

logger = logging.getLogger("molecule")

class Molecule():
    """ the Molecule base class """
    def __init__(self, elements=None, 
                    coordinates=None,
                    fragments=None,
                    charge_multi="0 1",
                    comment=""):
        logger.debug("Initialise new molecule")
        self.elements = elements
        self.fragments = fragments
        self.coordinates = np.array(coordinates)
        self.charge_multi = charge_multi.strip()
        self.comment = comment
        self.energy = None
        self.spin = None
        self.point = None
        self.lowest_freq = None
        self.nimag = None 
        self.d_mat = None # distance matrix

    def __str__(self):
        string=""
        string += self.charge_multi
        string += "\n"
        for el, coord in zip(self.elements,self.coordinates):
            string += "{:3s} {:14.8f} {:14.8f} {:14.8f} \n".format(el,*coord)
        return(string)

    def read_xyz(self,inpstr=None, process_comment = False):
        logger.debug("Reading coordinates...")
        if isinstance(inpstr,str):
            lines = inpstr.split("\n")
        elif isinstance(inpstr,TextIOWrapper):
            lines = inpstr.readlines()
        else:
            lines=0
        num_atoms = int(lines[0])
        self.comment = lines[1]
        if process_comment:
            self._process_comment_line(lines[1])
        coords=[]
        elements=[]
        for ii in range(num_atoms):
            temp = lines[ii+2].split()
            coords.append([float(temp[1]), float(temp[2]), float(temp[3]) ])
            elements.append(str(temp[0]))
        self.elements = elements
        self.coordinates = np.array(coords)
        logger.debug("Finished reading coordinates.")
    
    def write_xyz(self,filename, comment="", write_mode="w"):
        logger.debug(f"Writing xyz coordinates to file {str(filename)}...")
        if isinstance(filename,str):
            of = open(filename, write_mode)
            need_to_close = True
        else:
            of = filename 
            need_to_close = False
        num_atoms=len(self.coordinates)
        of.write("{:} \n".format(num_atoms))
        of.write("{:80s} \n".format(comment.strip()))
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

    def distance(self,atom1,atom2):
        vec = self.coordinates[atom1] - self.coordinates[atom2] 
        return(np.linalg.norm(vec))

    def distance_matrix(self,broadcast=False):
        """ 
        we compute the distance matrix utilising
        the broadcasting feature in numpy or utilising a simpler and slower code
                
        see: https://stackoverflow.com/questions/22720864/efficiently-calculating-a-euclidean-distance-matrix-using-numpy
        (see also: scipy.spatial.distance_matrix for an alternative using scipy)

        explanation of the broadcasting feature:
        # None is here identical to np.newaxis
        >>> a=np.array([[1,0,0],[1,0,1]])
        >>> a.shape
        (2, 3)
        >>> a[None,:,:].shape
        (1, 2, 3)
        >>> a[:,None,:].shape
        (2, 1, 3)
        >>> d=a[None,:,:]-a[:,None,:]
        >>> d
        array([[[ 0,  0,  0],
                [ 0,  0,  1]],

               [[ 0,  0, -1],
                [ 0,  0,  0]]])
        >>> d.shape
        (2, 2, 3)
        """
        if broadcast:
            b = self.coordinates[None,:,:] - self.coordinates[:,None,:]
            self.d_mat = np.linalg.norm(b, axis=2)
        else:
            N = len(self.coordinates)
            b=[]
            for na in range(N):
                b.append(self.coordinates - self.coordinates[na]) # this uses also numpy broadcasting
            b = np.array(b)
            self.d_mat = np.linalg.norm(b, axis=2)

    def center_coordinates(self):
        """ 
        Center coordinates using the geometric center of all atoms
        """
        centroid = self.centroid(self.coordinates)
        self.coordinates -= centroid

    @staticmethod
    def rmsd(coords1,coords2):
        """
        Compute the rmsd between two coordinate sets
        
        Input: 
            coords1,coords2 = np.array() (shape (m,3))

        Returns:
            rmsd = float()
        """
        pair_distances = np.linalg.norm(coords1 - coords2,axis=1)
        rmsd = np.sqrt( np.sum(pair_distances **2 )/ len(pair_distances) )
        return(rmsd)

    @staticmethod
    def centroid(coords):
        """
        Compute the centroid of coordinates sets
        
        Input: 
            coords = np.array() (shape (m,3))

        Returns:
            centroid = np.array (shape (1,3))
        """
        c = np.sum(coords, axis=0)
        l = len(coords)
        return(c/l)

    @staticmethod
    def kabsch(A,B):
        """
        Superimpose A on B using the kabsch algorithm
        
        see https://en.wikipedia.org/wiki/Kabsch_algorithm
        and https://en.wikipedia.org/wiki/Wahba%27s_problem 

        the kabsch algorithm has 3 steps:
        1. center both molecules A B
        2. compute the covariance Matrix H = A.T * B
        3. compute the rotation matrix R with single-value-decomposition:
            H = U S V.T
            d = sign(det(U V.T) Note: U and V are orthogonal
            R = U [[1,0,0],[0,1,0][0,0,d]] V.T
        Input: 
            A,B = np.array() (shape (m,3))

        Returns:
            A np.array() (shape (m,3))  : A rotated by R
        """
        # Center A and B
        A -= Molecule.centroid(A)
        B -= Molecule.centroid(B)
        
        # Compute covariance matrix h
        h =  A.T @ B 
        
        #single-value-decomposition of h
        u,s,vt = np.linalg.svd( h )
        
        # check for proper rotation
        d = np.sign( np.linalg.det( u ) * np.linalg.det(vt) )

        # compute rotation matrix r
        e = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,d]])
        r = u @ e @ vt
        A = A @ r
        return(A)

    def _process_comment_line(self,line):
        """
        Try to extract information from xyz comment lines.
        assuming the following format:
        E(UB3LYP) = -118.996452168 | S**2 = 1.008 | point 069 | freq -400.000 | Nimag = 2
        """
        temp = line.split("|")
        self.energy = float(temp[0].split()[2])
        self.spin = float(temp[1].split()[2])
        self.point = int(temp[2].split()[1])
        if len(temp) == 4: self.lowest_freq = float(temp[3].split()[1])
        if len(temp) == 5: self.nimag = float(temp[4].split()[2])
    
    def comment_line(self, *args):
        comment = f"{self.scf_energy} | S**2 = {self.spin:.3f}"
        if isinstance(self.point,int):
            comment += f" | point {self.point:03d}"
        return(comment)

    @staticmethod
    def read_xyz_trj(inpstr=None):
        """ read in an xyz trj file and return a list of the single xyz files """
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
    """
    ONIOM Molecule. 

    Additional Parameters:
        self.connectivity : list
        self.coordinates : list 
        self.parm : list 
        self.opt_flags : list 
        self.layer : list 
        self.link_atom : list 
        self.mm_elements : list 
        self.high_layer_idx : list 
        self.medium_layer_idx : list 
        self.low_layer_idx : list 
            idx - internal atom numbers which belong to respective layer.

    Additional Functions:
        read_oniom_template()
        write_xyz()
        write_layer_xyz()
        write
    
    """
    def read_oniom_template(self,inpstr):
        self.connectivity = []
        self.coordinates = []
        self.parm = []
        self.opt_flags = []
        self.layer = []
        self.link_atom = []
        self.mm_elements = []
    
        logger.debug("Reading ONIOM template...")
        if isinstance(inpstr,str):
            lines = inpstr.split("\n")
        elif isinstance(inpstr,TextIOWrapper):
            lines = inpstr.readlines()
        else:
            lines=0
        block=0
        for line in lines:
            temp = line.strip()
            if temp == "":
                block += 1
            elif block == 0:
                data = self._read_molecule_spec(temp)
                self.mm_elements.append(data[0])
                self.opt_flags.append(data[1])
                self.coordinates.append(data[2])
                self.layer.append(data[3])
                self.link_atom.append(data[4])

            elif block == 1:
                self.connectivity.append(line.strip())
            elif block == 2:
                self.parm.append(line.strip())
        self.coordinates = np.array(self.coordinates)    
        self.high_layer_idx = [ x for x,l in enumerate(self.layer) if l == "H" ]
        self.medium_layer_idx = [ x for x,l in enumerate(self.layer) if l == "M" ]
        self.low_layer_idx = [ x for x,l in enumerate(self.layer) if l == "L" ]
        self.active_atoms_idx = [ x for x,l in enumerate(self.opt_flags) if l == 0 ]
        logger.debug("Finished reading ONIOM Template.")
        logger.debug("ONIOM template: \n" + str(self))

    def _read_molecule_spec(self,line):
        element = ""
        opt_flag = ""
        coordinate = []
        layer = ""
        link_atom = ""
        temp = line.split()
        num_entries = len(temp)
        if num_entries == 4:
            element = temp[0]
            coordinate = [float(temp[1]),float(temp[2]),float(temp[3])]
        if num_entries >= 6:
            element = temp[0]
            opt_flag = int(temp[1])
            coordinate = [float(temp[2]),float(temp[3]),float(temp[4])]
            layer = temp[5]
        if num_entries >= 7:
            link_atom = temp[6]
        return(element, opt_flag, coordinate, layer, link_atom)

    def write_xyz(self, filename, comment="", write_mode="w"):
        super().write_xyz(filename, comment, write_mode)
        if isinstance(filename,TextIOWrapper):
            filename = filename.name
        self.write_layer_xyz(filename + "_high.xyz", comment=comment, write_mode=write_mode, symbol="H")
        self.write_layer_xyz(filename + "_medium.xyz", comment=comment, write_mode=write_mode, symbol="M")
        #self.write_layer_xyz(filename + "_low.xyz", comment=comment, write_mode=write_mode, symbol="L")

    def write_layer_xyz(self,filename, comment="", write_mode="w", symbol="H"):
        mol = Molecule()
        if symbol == "H": layer_idx = self.high_layer_idx
        elif symbol == "M": layer_idx = self.medium_layer_idx
        else: layer_idx = self.low_layer_idx
        mol.elements = [ self.elements[i] for i in layer_idx ]
        mol.coordinates  = [ self.coordinates[i] for i in layer_idx ]
        mol.comment = self.comment
        mol.write_xyz(filename,comment,write_mode)

    def update_opt_flags(self,atoms,flag=-1):
        """
        read in a list of atoms and set the opt flags to flag
        
        Parameters: atoms : list of int
                    flag : int
        """
        logger.debug("Updating opt flags..")
        for idx in range(len(self.opt_flags)):
            if idx in atoms:
                self.opt_flags[idx] = int(flag)
        
    
    def __str__(self):
        string=""
        string += self.charge_multi
        string += "\n"
        for el, opt_flags, coord, layer, link_atom in zip(self.mm_elements, self.opt_flags, self.coordinates, self.layer, self.link_atom):
            string += "{:20s} {:3d} {:14.8f} {:14.8f} {:14.8f} {:3s} {:3s} \n".format(el,opt_flags,*coord,layer,link_atom)
        string += "\n"
        for line in self.connectivity:
            string += "{} \n".format(line)
        string += "\n"
        for line in self.parm:
            string += "{} \n".format(line)
        return(string)


        

              

 