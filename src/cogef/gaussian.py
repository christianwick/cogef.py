
import logging

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