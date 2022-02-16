
import logging
import re
import sys

from collections import OrderedDict
from cogef.molecule import Molecule
from cogef.molecule import OniomMolecule
from cogef.mulliken import Mulliken

logger = logging.getLogger("gaussian")


class classLink0():
    def __init__(self,chk=None,mem="30GB",nproc="12"):
        self.chk = chk
        self.mem = mem
        self.nproc = nproc
    
    def __str__(self):
        lines = ""
        if self.chk:
            lines += "%chk=" + self.chk.strip() + "\n"
        lines += "%mem=" + self.mem.strip() + "\n"
        lines += "%nproc=" + self.nproc.strip() + "\n"
        lines += "\n"
        return(lines)

class classGaussianRoute():
    def __init__(self):
        self.route = ["#P"]
        self.opt = []
        self.guess = []
        self.geom = []
        self.scf = []
        self.iop = []

    def __str__(self):
        line = ""
        line += self._keywords_to_string(options = self.route,separator=" ") + " "
        line += self._keywords_to_string("opt", self.opt) + " "
        line += self._keywords_to_string("guess", self.guess) + " "
        line += self._keywords_to_string("geom", self.geom) + " "
        line += self._keywords_to_string("IOp", self.iop) + " "
        line += self._keywords_to_string("scf", self.scf) + " "
        return (line)

    @staticmethod
    def _keywords_to_string(key=None,options=None,separator=","):
        line = ""
        if isinstance(options,list) and options != []:
            for idx in range(len(options)-1):
                line += options[idx].strip() + separator
            line += options[-1] 
            if key:
                line = key + "(" + line + ")"
        return(line)
    
class GaussianInput():
    """ the gaussian base class 
    
    Parameters:
        of : open file handle
        mem : str   
        nproc : str
        charge_multi : str
    
    """
    def __init__(self, of=sys.stdout, mem="12GB", nproc="12", charge_multi = "0 1"):
        logger.debug("Initialise gaussian base class")
        self.of = of
        self.molecule = Molecule(charge_multi=charge_multi)
        self.comment = "gen by cogef.py"
        self.link0 = classLink0(mem=mem,nproc=nproc)   

    def write_inputfile(self,link1=False, route="#P",geom=True, modredundant=None):
        """
        general routine to write gaussian input files.

        Parameters: route : str
                        gaussian route section
                    geom : boolean
                        True : write molecule specification (incl charge and multiplicity)
                    modredundant : list
                        if not None, print modredundant section at the end.
        """
        # link1
        if link1: self.of.write("--link1--\n")
        # link0
        self.of.write(str(self.link0))
        # route
        self.of.write(str(route))
        self.of.write("\n")
        self.of.write("\n")
        # comment and
        # molecule specification
        if geom:
            self.of.write(self.comment.strip())
            self.of.write("\n")
            self.of.write("\n")
            self.of.write(str(self.molecule))
            self.of.write("\n")
        # modredundant
        if modredundant:
            self._write_modredundant(modredundant)
        self.of.write("\n")
    
    def _write_modredundant(self, args=["1 11 F"]):
        for val in args:
            self.of.write("{} \n".format(val))
        self.of.write("\n")

class OniomInput(GaussianInput):
    def __init__(self,charge_multi = "0 1 0 1 0 1", **kwargs):
        super().__init__(**kwargs)
        self.molecule = OniomMolecule(charge_multi = charge_multi)
        self.high_layer = Molecule(charge_multi = charge_multi)
    
    def _write_molecule(self):
        self.of.write(self.molecule.charge_multi)
        self.of.write("\n")
        for el, opt_flag, coord, layer, link_atom in zip(self.molecule.mm_elements, self.molecule.opt_flags, self.molecule.coordinates, self.molecule.layer, self.molecule.link_atom):
            self.of.write("{:20s} {:3d} {:14.8f} {:14.8f} {:14.8f} {:3s} {:3s} \n".format(el,opt_flag,*coord,layer,link_atom))
        self.of.write("\n")

    def _write_connectivity(self):
        for line in self.molecule.connectivity:
            self.of.write("{} \n".format(line))
        self.of.write("\n")
    
    def _write_parm(self):
        for line in self.molecule.parm:
            self.of.write("{} \n".format(line))
        self.of.write("\n")
    
    def write_inputfile(self,link1=False, route="#P",geom=True,connectivity=True, modredundant=None, parm=True):
        """
        general routine to write ONIOM input files.

        Parameters: route : str
                        gaussian route section
                    geom : boolean
                    connectivity : boolean
                    parm : boolean
                        True : write molecule specification (incl charge and multiplicity)
                            / connectivity / parameters
                    modredundant : list
                        if not None, print modredundant section at the end.
        """
        # link1
        if link1: self.of.write("--link1--\n")
        # link0
        self.of.write(str(self.link0))
        # route
        self.of.write(str(route))
        self.of.write("\n")
        self.of.write("\n")
        # comment and
        # molecule specification
        if geom:
            self.of.write(self.comment.strip())
            self.of.write("\n")
            self.of.write("\n")
            self._write_molecule()
        if connectivity:
            self._write_connectivity()
        # modredundant
        if modredundant:
            self._write_modredundant(modredundant)
        if parm:
            self._write_parm()


class AmberInput(OniomInput):
    """ the gaussian base class for pure AMBER calculations """

    def write_input(self,of,modredundant, initial_stab_opt = False, instability=False, route_args={}, oniom_opt = False):
        logger.debug("Writing gaussian input file ...")
        """ write input files for oniom calculations

        Input:
            we ignore most of the input arguments. mainly we use the route_args.
        """
        # NOTE: ONIOM uses RFO as standard. we use the gaussian standard optimizer
        # we setup the opt args assuming that we start a simple optimisation. 
        # if we perform other steps prior to optimisation we update them in each section accordingly
        opt_args = {
            "guess" : ["Mix","always"],
            "geom" : ["connectivity"],
            "opt" : ["NoMicro"],
            "nosymm" : None,
            "scf" : ["XQC","MaxConven=75"] }
        opt_args.update(route_args)
        # SECTION 2: OPTIMISATION
        self._write_link0(of)
        self._write_route(of, args = opt_args)
        self._write_title(of)
        self._write_molecule(of)
        #self._write_modredundant(of,modredundant)
        self._write_parm(of)
        of.write("\n\n")
        logger.debug("Finished writing gaussian input file.")
        
class CheckGaussianLogfile():
    """ Process gaussian log files and check for errors

        self.instability -> boolean (True if an instability has been found)
        self.stationary -> boolean (True if stationary point has been found)
        self.error -> boolean (True if last termination was not a Normal termination.)
    """
    def __init__(self, filename=None):
        logger.debug("Initialise CheckGaussianLogfile class")
        self._find_normal_termination = re.compile("Normal termination")
        self._find_instab = re.compile("The wavefunction has an internal instability.")
        self._find_stable = re.compile("The wavefunction is stable under the perturbations considered.")
        self._find_stationary = re.compile("-- Stationary point found.")

        self._find_spin = re.compile("S\*\*2 before annihilation ")
        self._find_scf_energy = re.compile("SCF Done:")
        self._find_struct = re.compile("Input orientation")
        self._find_mulliken = re.compile("Mulliken charges and spin densities:")
        self._find_mulliken_h = re.compile("Mulliken charges and spin densities with hydrogens summed into heavy atoms:")

        self.filename = filename
        self.instability = False
        self.stationary = False
        self.error = False

        self.spin = 0.0
        self.scf_energy = 0.0
        self.molecule = Molecule()
        self.mulliken = Mulliken()
    
    def read_log(self):
        self.instability = False
        self.stationary = False
        self.error = False
        self.spin = 0.0
        self.scf_energy = 0.0
        
        last_line = ""
        logger.debug(f"Reading gaussian log file {self.filename}")
        with open(self.filename,"r") as of:
            while True:
                line = of.readline()
                if not line: break
                # check for errors
                self._read_stationary(line)
                self._read_instability(line)

                # check for important data
                self._read_spin(line)
                self._read_scf_energy(line)
                self._read_struct(line, of)
                self._read_mulliken_h(line, of)
                last_line = line
            # we check only the last line for error terminations:
            self._check_termination(last_line)
        logger.info(f"WFN is stable: {not self.instability}")
        logger.info(f"Found stationary point: {self.stationary}")
        logger.info(f"Normal termination: {not self.error}")
        logger.info(f"S**2 = {self.spin}")
        logger.info(f"SCF Energy: {self.scf_energy}")

    def check_bonding(self):
        logger.debug(f"Check bonding...")
        mul_spin_dens = self.mulliken.find_spin_larger_than(0.5)
        if len(mul_spin_dens.atom_number) >= 2:
            logger.info(f"Mulliken Spin densities larger than 0.5: \n{mul_spin_dens}")
            # get atom numbers of min and max spin dens
            mul_spin_dens.sort_by("spin")
            min_mulliken = mul_spin_dens.atom_number[0] - 1
            max_mulliken = mul_spin_dens.atom_number[-1] - 1 
            dist = self.molecule.distance(min_mulliken,max_mulliken)
            logger.info(f"Distance between atoms {min_mulliken +1 } and {max_mulliken + 1} = {dist:10.3f} A")
        logger.debug(f"Finished check bonding.")
        
    def _read_instability(self,line):
        if self._find_instab.search(line):
            logger.warning(f"Instability detected in file {self.filename}")
            self.instability = True
        #elif self._find_stable.search(line):
        #    logger.debug(f"Detected a stable Wavefunction in file {self.filename}")
        #    self.instability = False

    def _read_stationary(self,line):
        if self._find_stationary.search(line):
            logger.debug(f"Stationary point found in file {self.filename}")
            self.stationary = True

    def _check_termination(self,line):
        if self._find_normal_termination.search(line):
            logger.debug(f"Normal termination detected in file {self.filename}")
            self.error = False
        else: 
            logger.warning(f"Error termination detected in file {self.filename}")
            self.error = True 
    
    def _read_spin(self, line):
        if self._find_spin.search(line):
            self.spin = float(line.split()[3].strip(","))
    
    def _read_scf_energy(self, line):
        if self._find_scf_energy.search(line):
            temp = line.split()
            self.scf_energy = temp[2] + " = " + temp[4]

    def _read_struct(self, line, fin):
        """ 
        Read structure from gaussian log files.
        
        input:  fin = open file read
                line = string
        """
        if self._find_struct.search(line):
            self.molecule.read_gaussian_raw_coords(self._read_geom(fin))

    def _read_mulliken_h(self, line, of):
        """ 
        Read MULLIKEN  data from gaussian log files with hydrogens summed into
        heavy atoms
        
        input:  fin = open file read
                line = string
        """
        if self._find_mulliken_h.search(line):
            of.readline() # this is the line with 1 2
            self.mulliken._read_gaussian_raw_data(
                    self._read_mulliken_blocks(of,"Electronic spatial extent"))

    def comment_line(self, point=None, *args):
        comment = f"{self.scf_energy} | S**2 = {self.spin:.3f}"
        if isinstance(point,int):
            comment += f" | point {point:03d}"
        return(comment)
    
    @staticmethod
    def _read_mulliken_blocks(of, term = ""):
        """ 
        Read Mulliken data from gaussian log files
        Input:
            of = open file read
            term = str 
        The file "of" will be read until "term" is found.
        """
        temp = []
        while True:
            line = of.readline()
            if term in line or line == "":
                break
            temp.append(line)
        return(temp)

    @staticmethod
    def _read_geom(fin,reg=re.compile("----")):
        """
        Read geometry blocks in gaussian log files

        input:  fin = open file read
                reg = regex

        returns: raw_coords = list
        """
        raw_coords = []
        fin.readline()
        one = fin.readline().split(maxsplit=3)
        two = fin.readline().split(maxsplit=3)
        #title = [ one[x] + " " + two[x].strip() for x in range(4)  ]
        fin.readline()
        while True:
            line = fin.readline()
            if ( not line ) or reg.search(line): break 
            raw_coords.append(line)
        return(raw_coords)


class CheckOniomLogfile(CheckGaussianLogfile):
    """ Process gaussian log files and check for errors in ONIOM calculations

        self.instability -> boolean (True if an instability has been found)
        self.stationary -> boolean (True if stationary point has been found)
        self.error -> boolean (True if last termination was not a Normal termination.)
    """
    def __init__(self, *args):
        super().__init__(*args)
        self._find_scf_energy = re.compile("ONIOM: extrapolated energy =")
    
    def _read_scf_energy(self, line):
        if self._find_scf_energy.search(line):
            temp = line.split()
            self.scf_energy = temp[0] + " = " + temp[4]

class CheckAmberLogfile(CheckGaussianLogfile):
    """ Process gaussian log files and check for errors in pure AMBER calculations

        self.instability -> boolean (True if an instability has been found)
        self.stationary -> boolean (True if stationary point has been found)
        self.error -> boolean (True if last termination was not a Normal termination.)
    """
    def __init__(self, *args):
        super().__init__(*args)
        self._find_scf_energy = re.compile("^ Energy=")
    
    def _read_scf_energy(self, line):
        if self._find_scf_energy.search(line):
            temp = line.split()
            self.scf_energy = temp[0] + " = " + temp[1]

#the following functions are kept for backwards compatibility and will
#be removed in future versions
def read_instability_from_log(filename):
    instab=False
    with open(filename,"r") as of:
        while True:
            line = of.readline()
            if not line: break 
            if "The wavefunction has an internal instability." in line:
                instab = True
    return(instab)

def read_stationary_from_log(filename):
    stationary=False
    with open(filename,"r") as of:
        while True:
            line = of.readline()
            if not line: break 
            if "-- Stationary point found." in line:
                stationary = True
    return(stationary)

def check_normal_termination_from_log(filename):
    last_line = ""
    error = True
    with open(filename,"r") as of:
        while True:
            line = of.readline()
            if not line: break
            else: last_line = line
    if "Normal termination" in last_line:
        error = False
    return error