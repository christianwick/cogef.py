
import logging
import re

from collections import OrderedDict
from cogef.molecule import Molecule

logger = logging.getLogger("gaussian")

class GaussianInput():
    """ the gaussian base class """
    def __init__(self, filename="", 
                    title="Gen by cogef",
                    chkfile="cogef.chk",
                    mem="12GB",
                    nproc="12",
                    elements=[], 
                    coordinates=[],
                    fragments=[],
                    charge=[0],
                    multiplicity=[1],
                    route = "#P",
                    method="UB3LYP/D95(d,p)"):
        logger.debug("Initialise gaussian base class")
        self.filename = filename
        self.molecule = Molecule(
                    elements=elements, 
                    coordinates=coordinates,
                    fragments=fragments,
                    charge=charge,
                    multiplicity=multiplicity)
        self.title = title
        self.link0 = OrderedDict( {
                "%chk" : chkfile,
                "%mem" : mem,
                "%nproc" : nproc })
        self.route = OrderedDict({
                "route" : route ,
                "level_of_theory" : method })                      

    def _write_link0(self, of, args={}):
        temp = self.link0.copy()
        temp.update(args)
        for arg,val in temp.items():
            of.write("{}={}\n".format(arg,val))
        of.write("\n")

    def _write_route(self, of, args={}, method=None):
        temp = self.route.copy()
        temp.update(args) 
        if isinstance(method,str):
            temp.update( {"level_of_theory" : method })
        route = temp.pop("route")
        lot = temp.pop("level_of_theory")
        of.write("{} {} ".format(route,lot))
        for arg, val in temp.items():
            if val == None:
                of.write("{} ".format(arg))
            elif isinstance(val,list):
                of.write("{}=(".format(arg))
                for el in val[:-1]:
                    of.write("{},".format(el))
                of.write("{}) ".format(val[-1]))
            else:
                of.write("{}={} ".format(arg,val))
        of.write("\n\n")

    def _write_title(self, of, title=None):
        if isinstance(title,str):
            of.write("{}\n".format(title.strip()))
        else:
            of.write("{}\n".format(self.title.strip()))
        of.write("\n")

    def _write_molecule(self,of,molecule=None):
        if not isinstance(molecule, Molecule):
            molecule=self.molecule
        of.write(str(molecule))
        of.write("\n")
    
    def _write_modredundant(self, of, args=["1 11 F"]):
        for val in args:
            of.write("{} \n".format(val))
        of.write("\n")
    
    def _write_link1(self,of):
        of.write("--link1--\n")

    def write_input(self,of,modredundant, initial_stab_opt = False, instability=False, route_args={}):
        logger.debug("Writing gaussian input file ...")
        opt_args = {
            "guess" : ["Mix","always"],
            "geom" : "Modredundant",
            "opt" : ["CalcFc","RFO"],
            "nosymm" : None,
            "scf" : ["XQC","MaxConven=75"] }
        opt_args.update(route_args)
        # should we perform a stability analysis prior to optimisation?
        if initial_stab_opt:
            opt_args.update({"guess" : "read", "geom" : ["Modredundant","allcheck"]})
            self._write_link0(of)
            self._write_route(of, args = {
                "guess" : "mix",
                "stable" : "opt",
                "scf" :  ["XQC","MaxConven=75"] ,
                "nosymm" : None })
            self._write_title(of)
            self._write_molecule(of)
            self._write_link1(of)
        # in case we found an instability we use the stable wfn directly.        
        if instability:
            opt_args.update({"guess" : "read", "geom" : ["Modredundant","allcheck"]})
        self._write_link0(of)
        self._write_route(of, args = opt_args)
        if not "allcheck" in opt_args["geom"]:
            self._write_title(of)
            self._write_molecule(of)
        self._write_modredundant(of,modredundant)
        self._write_link1(of)
        self._write_link0(of)
        self._write_route(of, args = {
            "guess" : "read",
            "geom" : "allcheck",
            "stable" : "opt",
            "scf" :  ["XQC","MaxConven=75"] ,
            "nosymm" : None })
        of.write("\n\n")
        logger.debug("Finished writing gaussian input file.")

    def write_sp_input(self, of, route_args={}):
        logger.debug("Writing gaussian sp input file ...")
        opt_args = {
            "guess" : ["Mix"],
            "stable" : "opt",
            "scf" : "XQC",
            "nosymm" : None }
        opt_args.update(route_args)
        self._write_link0(of)
        self._write_route(of, args = opt_args)
        self._write_title(of)
        self._write_molecule(of)
        of.write("\n\n")
        logger.debug("Finished writing gaussian sp input file.")


    

class GaussianInputWithFragments(GaussianInput):  
    def _write_molecule_with_fragments(self,of,molecule=None,fragment=None):
        if not isinstance(molecule, Molecule):
            molecule=self.molecule
        if fragment == None:
            of.write(str(molecule))
        else:
            for ii in range(len(molecule.elements)):
                if ii in fragment:
                    molecule.elements[ii] += "(Fragment=1)"
                else:
                    molecule.elements[ii] += "(Fragment=2)"
            for i,j in zip(molecule.charge,molecule.multiplicity):
                of.write("{} {} ".format(i,j))
            of.write("\n")
            for el, coord in zip(molecule.elements,molecule.coordinates):
                of.write("{:3s} {:14.8f} {:14.8f} {:14.8f} \n".format(el,*coord))
        of.write("\n")

    def write_input_fragment_guess(self,of,modredundant,fragment):
        logger.debug("Writing gaussian input file with fragments...")
        self._write_link0(of)
        self._write_route(of, args = {"guess": ["fragment=2","only"]})
        self._write_title(of)
        self._write_molecule_with_fragments(of,fragment=fragment)
        self._write_link1(of)
        self._write_link0(of)
        self._write_route(of, args = {
            "guess" : ["read"],
            "geom" : ["Modredundant","allcheck"],
            "opt" : ["CalcFc","RFO","maxcyc=100"],
            "nosymm" : None,
            "scf" : "XQC" })
        self._write_modredundant(of,modredundant)
        self._write_link1(of)
        self._write_link0(of)
        self._write_route(of, args = {
            "guess" : "read",
            "geom" : "allcheck",
            "stable" : "opt",
            "scf" : "XQC",
            "nosymm" : None })
        of.write("\n\n")
        logger.debug("Finished writing gaussian input file with fragments.")


class OniomInput(GaussianInput):
    def __init__(self, template, **kwargs):
        super().__init__(**kwargs)
        pass
        #self.molecule = OniomMolecule(): 
        
class CheckGaussianLogfile():
    """ Process gaussian log files and check for errors

        self.instability -> boolean (True if an instability has been found)
        self.stationary -> boolean (True if stationary point has been found)
        self.error -> boolean (True if last termination was not a Normal termination.)
    """
    def __init__(self, filename):
        logger.debug("Initialise CheckGaussianLogfile class")
        self._find_normal_termination = re.compile("Normal termination")
        self._find_instab = re.compile("The wavefunction has an internal instability.")
        self._find_stable = re.compile("The wavefunction is stable under the perturbations considered.")
        self._find_stationary = re.compile("-- Stationary point found.")

        self._find_spin = re.compile("S\*\*2 before annihilation ")
        self._find_scf_energy = re.compile("SCF Done:")
        self._find_struct = re.compile("Input orientation")

        self.filename = filename
        self.instability = False
        self.stationary = False
        self.error = False

        self.spin = 0.0
        self.scf_energy = 0.0
        self.molecule = Molecule()
    
    def read_log(self):
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
                last_line = line
            # we check only the last line for error terminations:
            self._check_termination(last_line)
        logger.info(f"WFN is stable: {not self.instability}")
        logger.info(f"Found stationary point: {self.stationary}")
        logger.info(f"Normal termination: {not self.error}")
        logger.info(f"S**2 = {self.spin}")
        logger.info(f"SCF Energy: {self.scf_energy}")

    def _read_instability(self,line):
        if self._find_instab.search(line):
            logger.warning(f"Instability detected in file {self.filename}")
            self.instability = True
        elif self._find_stable.search(line):
            logger.debug(f"Detected a stable Wavefunction in file {self.filename}")
            self.instability = False

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

    def comment_line(self, point=None, *args):
        comment = f"{self.scf_energy} | S**2 = {self.spin:.3f}"
        if isinstance(point,int):
            comment += f" | point {point:03d}"
        return(comment)

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