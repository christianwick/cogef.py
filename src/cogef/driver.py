""" molecule.py contains all Molecule classes """

from re import S
import numpy as np
import sys
import os

import logging

from io import TextIOWrapper
import subprocess

from cogef.constants import PSE
from cogef import gaussian
from cogef import modstruct
from cogef import cogef_logging

logger = logging.getLogger("driver")

class cogef_loop():
    """
    this is the base class for cogef calculations.
    we set all important options on init.

    functions:  run(): run the cogef calculation.

    helper_functions:
                check_gaussian_logfile() 
                run_gauss()

    """
    def __init__(self, atom1, atom2, dx = 0.1, runtype=None,  xyz = None,
                    level_of_theory = "B3LYP/6-31G*" , mem = "4GB", nproc = "12", print_level="#P",
                    maxcyc = None, maxconv = 75, startchk = False, cm = "0 1",
                    cycles = 1, reverse = None, restart = 0, modredundant=None, symm_stretch=True, 
                    dp=1.0, fragment=[], max_error_cycles = 5, mulliken_h = False,
                    trajectory = None, checkpoint = None , no_mix = False, **kwargs):
        """
        initalise all important parameters

        Parameters: 
                    runtype : str 
                        choices: "restricted/unrestricted/None"
                        this will set the formalism and dependend checks / input keywords
                    xyz : TextIOwrapper,None
                        xyz file to read. will be read if specified.
                    atom1 : int 
                    atom2 : int 
                        atoms to stretch starting at 0
                    dx : float
                        stretch increment
                    level_of_theory : str
                    mem = str
                    nproc = str
                    print_level = str
                        gaussian print level, e.g. #P, #N
                    MaxCyc : int, None
                        gaussian opt option maxcyc
                    MaxConv : int, None
                        gaussian scf=xqc option maxconventional
                    cm = str
                        charge multiplicity line
                    startchk : boolean
                        if True, read guess from guess.chk file at first iteration.
                    cycles : int 
                        cycles to run 
                    reverse : Boolean
                        if True, run in other direction usind abs(dx)
                    restart : int
                        restart a previuous calculation at cycle
                    modredunant : list
                        list of additional modredundant options. elements e.g. "1 2 F"
                    symm_stretch : boolean 
                        stretch symmetric (True) or asymmetric (False)
                    dp : float
                        dp used to scale the stretching of the fragment atoms
                    fragment : list 
                        contains all the atoms of the first fragment to be stretched.
                    max_error_cycles : int
                        set the maximum number of errors to restart each cycle.
                    mulliken_h : boolean
                        if True, use mulliken charges with Hydrogen atoms summed into
                        heavy atoms to check bonging. use standard Mulliken charges otherwise.
                    trajectory : of,None
                    checkpoint : of,None
                        write trajectory / checkpoint file if points to an open file
                    no_mix : Bool
                        If False, allow mixing of initial HOMO,LUMO orbitals, 
                        if True, turn of mixing completely.

        internal Parameters:
                    self.check_stability : boolean   
                    self.check_error : boolean
                    self.check_stationary : boolean
                        check for stability / error termination / stationary 
                            point if True
                    self.allow_mixing : Bool
                        is the opposite of no_mix. Eg. if no_mix = True, allow_mixing = False
        
        """
        self.level_of_theory = level_of_theory.strip()
        self.startchk = startchk
        self.dx = dx 
        self.atom1 = atom1
        self.atom2 = atom2
        self.modredundant = ["{} {} F".format(atom1 + 1, atom2 + 1 )]
        self.symm_stretch= symm_stretch
        self.dp = dp
        self.fragment = fragment
        self.max_error_cycles = max_error_cycles
        self.glog = gaussian.CheckGaussianLogfile()
        self.ginp = gaussian.GaussianInput(mem=mem, nproc = nproc, charge_multi = cm)
        self.print_level = print_level
        self.maxcyc = maxcyc
        self.maxconv = maxconv
        self.mulliken_h = mulliken_h
        self.allow_mixing = not no_mix
        if xyz:
            self.ginp.molecule.read_xyz(xyz)
        self.trajectory = trajectory
        self.checkpoint = checkpoint
        if runtype == "restricted":
            self.write_ginp = self._write_ginp_restricted # here we are not calling the function!
            self.check_error = True
            self.check_stationary = True
            self.check_stability = False
            self.mix_guess = False
        elif runtype == "unrestricted":
            self.write_ginp = self._write_ginp_unrestricted # here we are not calling the function!
            self.check_error = True
            self.check_stationary = True
            self.check_stability = True
            self.mix_guess = True
        if modredundant:
            self.modredundant.extend(modredundant)
        # loop_range for cogef main loop
        if reverse:
            # we loop from reverse n cycles backwards. smallest point is 000
            # and set dx to a negative value
            # and finally tell the glog class that we already found a broken
            # structure to avoid printing of another structure with a broken bond
            # by setting self.glog.found_broken_bond = True
            self.loop_range = range(reverse, max(-1, reverse - cycles -1 ), -1)
            self.dx = -abs(self.dx)
            self.glog.found_broken_bond = True
        else:
            self.loop_range = range(restart, cycles +1 )

    def run(self):
        """
        run the cogef calculation.

        We use two nested loops, the outer loop is increasing the actual cogef cycles. 
        At each cycle of the outer loop, the distance between atom1 and atom2 will be 
        increased (decreased)
        At each cycle of the inner loop, we check the convergence / regular termination of
        the cogef calculation.

        Internal Parameters:
            read_guess : boolean
                if True, we will set reading of the guess from the wfn. 
            mix_guess : boolean
                if True, we will mix HOMO LUMO orbitals of the guess orbs. This will only affect 
                the runtype == unrestricted calculations.
        """
        # start the outer cogef loop
        read_guess = False
        mix_guess = False
        for cycle in self.loop_range:
            logger.info("Starting cycle {} out of a maximum of {}".format(cycle,max(self.loop_range)))
            basename = "cogef_{:03d}".format(int(cycle))
            glog_filename = basename + ".log"
            ginp_filename = basename + ".com"   
            if self.allow_mixing: mix_guess = True   
            error_cycle = 0
            # start the inner loop 
            while error_cycle < self.max_error_cycles:
                logger.info("Starting iteration {} of {} on cycle {}".format(error_cycle,self.max_error_cycles,cycle))
                # write the gaussian input file. the actual function will be set at init
                # read_guess = True will read the guess from the chk file. 
                # mix_guess = True will mix the input HOMO and LUMO orbitals.  
                if cycle > 1 or self.startchk: read_guess = True 
                self.write_ginp(ginp_filename, read_guess=read_guess,mix_guess=mix_guess)
                self.run_gauss(ginp_filename)
                # checking the logfile for errors and updated coordinates:
                new_coordinates, glog_is_ok = self.check_gaussian_logfile(glog_filename,cycle,error_cycle)
                self.ginp.molecule.coordinates = new_coordinates
                if glog_is_ok:
                    if self.trajectory:
                        logger.info(f"Adding structure to trajectory {self.trajectory}")
                        self.ginp.molecule.write_xyz(self.trajectory, comment=self.glog.comment_line(point=cycle), write_mode="a")
                    break
                else:
                    os.rename(glog_filename, "{}_error_{}.log".format(basename, error_cycle))
                    error_cycle += 1
                    mix_guess = False # turn off mixing and read 
            # Modify coords for next cycle
            self.ginp.molecule.coordinates = modstruct.mod_fragments(self.ginp.molecule.coordinates, 
                               self.atom1, self.atom2, self.dx, symmetric=self.symm_stretch, dp=self.dp, fragment=self.fragment)
            # write checkpoint structure to disk.
            if self.checkpoint:
                logger.info(f"Writing checkpoint xyz file {self.checkpoint}")
                self.ginp.molecule.write_xyz(self.checkpoint,comment=f"input structure for cycle {cycle+1}")

    def _write_ginp_restricted(self, ginp_filename="", read_guess=False, **kwargs):
        """
        generate input file for restricted cogef calculations. This will NOT CHECK the stability of the 
        wfn at each cycle.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, use the previous wfn as guess.
        """
        logger.info(f"setting chk point file to {'guess.chk'}")
        self.ginp.link0.chk="guess.chk"    
        with open(ginp_filename, "w") as of:
            self.ginp.of = of
            # JOB 1 geometry optimisation
            route = gaussian.classGaussianRoute()
            route.route = [self.print_level, self.level_of_theory , "nosymm", "test"] 
            route.opt = ["modredundant", "RFO"]
            if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
            if read_guess: route.guess = ["read"]
            self.ginp.write_inputfile(link1=False,route=route,geom=True,modredundant=self.modredundant)

    def _write_ginp_unrestricted(self, ginp_filename="", read_guess=False, mix_guess=False):
        """
        generate input file for unrestricted cogef calculations. This will also check the stability of the 
        wfn at each cycle and restart the optimisation if the wfn is not stable.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, use the previous wfn as guess.
                mix_guess : boolean
                        if True, mix the inital orbitals.
        """
        logger.info(f"setting chk file: guess.chk")   
        self.ginp.link0.chk="guess.chk"    
        with open(ginp_filename, "w") as of:
            self.ginp.of = of
            # JOB 1 geometry optimisation
            route = gaussian.classGaussianRoute()
            route.route = [self.print_level, self.level_of_theory , "nosymm", "test"]
            route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
            route.opt = ["modredundant","RFO"]
            if mix_guess: route.guess.append("mix")
            if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
            if read_guess: route.guess.append("read")
            self.ginp.write_inputfile(link1=False,route=route,geom=True,modredundant=self.modredundant)
            # JOB 2 stability analysis
            route = gaussian.classGaussianRoute()
            route.route = [self.print_level, self.level_of_theory , "nosymm", "stable=opt", "test"]
            route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
            route.guess = ["read"]
            route.geom = ["allcheck"]
            self.ginp.write_inputfile(link1=True,route=route,geom=False)

    def check_gaussian_logfile(self,filename,cycle,error_cycle):
        """
        Check gaussian log files for error terminations.

        Parameters: filename : filename of the glog file to parse
                    cycle : int   current cycle
                    error_cycle : int current error cycle


        Returns:    coordinates : last coordinates from glog file
                    glog_is_ok : boolean
                            True if no error was found, False otherwise
        """
        self.glog.filename = filename
        self.glog.read_log(self.mulliken_h)
        if self.glog.check_bonding():
            self.glog.molecule.write_xyz(f"broken_bond_at_{cycle}.xyz", 
                comment=self.glog.comment_line(point=cycle))
        # check for instability or break out 
        # we also check for stationary points and move the old log files
        glog_is_ok = True
        if self.check_stability and self.glog.instability:
            logger.warning("Instability detected at cycle {} {}".format(cycle,error_cycle))
            glog_is_ok = False
        if self.check_error and self.glog.error:
            logger.warning("Error termination detected at cycle {} {}".format(cycle,error_cycle ))
            glog_is_ok = False
        if self.check_stationary and not self.glog.stationary:
            logger.warning("No stationary point found at cycle {} {}".format(cycle,error_cycle ))
            glog_is_ok = False 
        return ( self.glog.molecule.coordinates, glog_is_ok )
                         


    @staticmethod
    def run_gauss(ginp_filename):
        rungauss = subprocess.run(["cogef_rung16",ginp_filename], check=True, capture_output=True, text=True)
        



class oniom_cogef_loop(cogef_loop):
    """
    this modifies the cogef_loop class to run ONIOM type cogef calculations.

    additional Parameters: (to all parameters of the base class)
        oniomtemplate : of
        oniomopt : string 
            select oniom opt method. choices "hybrid", "me", "ee_sp"
        constraint : string
            choices: "modredundant", "opt_flag"
            select method used to constrain atoms in the cogef calculation. 

    """
    def __init__(self, oniomtemplate, oniomopt="hybrid", constraint="modredundant", xyz=None, **kwargs):
        super().__init__(**kwargs)
        self.oniomopt = oniomopt 
        self.constraint = constraint
        self.glog = gaussian.CheckOniomLogfile()
        self.ginp = gaussian.OniomInput(mem=self.ginp.link0.mem, nproc = self.ginp.link0.nproc,
             charge_multi = self.ginp.molecule.charge_multi)
        self.ginp.molecule.read_oniom_template(oniomtemplate)

        if xyz:
            self.ginp.molecule.read_xyz(xyz)
        if self.constraint == "opt_flag":
            logger.debug("setting opt flags..")
            self.convert_modredundant_to_opt_flags()
            self.modredundant = None

    def convert_modredundant_to_opt_flags(self):
        """
        convert modredundant input of fixed bonds e.g. "1 2 F" to 
        opt flags with -1 to fix atoms in an oniom calculation. 
        this is particularly necessary if you want to fix atoms of the real system
        which are not part of the high layer.
        """
        atoms = []
        for val in self.modredundant:
            temp = val.split()
            atoms.append(int(temp[0])-1)
            atoms.append(int(temp[1])-1)
        self.ginp.molecule.update_opt_flags(atoms)

    def _write_ginp_restricted(self, ginp_filename="", read_guess=False, **kwargs):
        """
        generate input file for restricted cogef calculations. This will NOT CHECK the stability of the 
        wfn at each cycle. Additionally, in ONIOM, we select the optimisation method:

            me : run an me only optimisation
            ee : run an ee only optimisation
            ee_sp : run an me optimisation followed by an ee single point calculation.
            hybrid : run an me optimisation followed by an ee optimisation

        The constraint used, can be added as a modredundant section or using the optimisation flag.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, read wfn from guess.chk
        """
        logger.info(f"setting chk point file to {'guess.chk'}")
        self.ginp.link0.chk="guess.chk"
        # select optimisation method and prepare input accordingly
        # we keep here some redundancy to increase readability and 
        # to facilitate further modifications.
        # -------------
        # HYBRID method
        if self.oniomopt == "hybrid":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom opt me
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.opt = ["Mic120","Quadmac"]
                if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess =  ["read"]
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)
                # JOB 2 oniom opt ee
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory+"=EMBEDCHARGE", "nosymm", "test"] 
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.guess = ["read"]
                route.opt = ["Quadmac"]
                if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
        # ------------
        # EE SP method
        elif self.oniomopt == "ee_sp":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom me opt
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.opt = ["Mic120","Quadmac"]
                if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess =  ["read"]
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)
                # JOB 2  oniom ee sp
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory+"=EMBEDCHARGE", "nosymm", "test"] 
                route.iop = ["2/15=3"]
                route.geom = ["allcheck"]
                route.guess = ["read"]
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
        # ---------
        # EE method
        elif self.oniomopt == "ee":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory+"=EMBEDCHARGE", "nosymm", "test"] 
                route.iop = ["2/15=3"]
                route.geom = ["connect"]
                route.opt = ["Quadmac"]
                if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess =  ["read"]
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)
        # ---------
        # ME method
        else:
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm", "test"] 
                route.iop = ["2/15=3"]
                route.geom = ["connect"]
                route.opt = ["Mic120", "Quadmac"]
                if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess =  ["read"]
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)

    def _write_ginp_unrestricted(self, ginp_filename="", read_guess=False, mix_guess = False):
        """
        generate input file for restricted cogef calculations. This will also check the stability of the 
        wfn at each cycle at the highest level applied in the optimisation method:

            me : run an me only optimisation
            ee : run an ee only optimisation
            ee_sp : run an me optimisation followed by an ee single point calculation.
            hybrid : run an me optimisation followed by an ee optimisation

        The constraint used, can be added as a modredundant section or using the optimisation flag.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, read guess from guess.chk
        """
        logger.info(f"setting chk point file to {'guess.chk'}")
        self.ginp.link0.chk="guess.chk"
        # select optimisation method and prepare input accordingly
        # we keep here some redundancy to increase readability and 
        # to facilitate further modifications.
        # -------------
        # HYBRID method
        if self.oniomopt == "hybrid":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom opt me
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.opt = ["Mic120","Quadmac"]
                if mix_guess: route.guess.append("mix")
                if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess.append("read")
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)
                # JOB 2 oniom opt ee
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory+"=EMBEDCHARGE", "nosymm","test"] 
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.opt = ["Quadmac"]
                route.guess.append("read")
                if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
                # JOB 3 oniom ee stable=opt
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory + "=EMBEDCHARGE", "nosymm", "stable=opt", "test"]
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.guess = ["read"]
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
        # ------------
        # EE SP method
        elif self.oniomopt == "ee_sp":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom opt me
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.opt = ["Mic120","Quadmac"]
                if mix_guess: route.guess.append("mix")
                if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess.append("read")
                self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                        modredundant=self.modredundant,parm=True)
                # JOB 2 oniom me stable=opt
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory + "=EMBEDCHARGE", "nosymm", "stable=opt", "test"]
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.guess = ["read"]
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
        # ---------
        # EE method
        elif self.oniomopt == "ee":
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom opt ee
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory+ "=EMBEDCHARGE", "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.opt = ["Quadmac"]
                if mix_guess: route.guess.append("mix")
                if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess.append("read")
                # JOB 2 oniom ee stable=opt
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory + "=EMBEDCHARGE", "nosymm", "stable=opt", "test"]
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.guess = ["read"]
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
        # ---------
        # ME method
        else:
            with open(ginp_filename, "w") as of:
                self.ginp.of = of
                # JOB 1 oniom opt me
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm","test"] 
                route.geom = ["connect"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.opt = ["Mic120","Quadmac"]
                if mix_guess: route.guess.append("mix")
                if self.maxcyc: route.opt.append("MaxCyc="+str(self.maxcyc))
                if self.constraint == "modredundant": route.opt.append("modredundant")
                if read_guess: route.guess.append("read")
                # JOB 2 oniom ee stable=opt
                route = gaussian.classGaussianRoute()
                route.route = [self.print_level, self.level_of_theory, "nosymm", "stable=opt", "test"]
                route.geom = ["allcheck"]
                route.iop = ["2/15=3"]
                route.scf = ["XQC", "Maxconv="+str(self.maxconv)]
                route.guess = ["read"]
                self.ginp.write_inputfile(link1=True,route=route, geom=False, connectivity=False, 
                        modredundant=self.modredundant,parm=False)
    

class amber_cogef_loop(oniom_cogef_loop):
    """
    this modifies the cogef_loop class to run AMBER type cogef calculations with gaussian.

    additional Parameters: (to all parameters of the base class)
        ambertemplate : of

        constraint : string
            choices: "modredundant", "opt_flag"
            select method used to constrain atoms in the cogef calculation. 

    """
    def __init__(self, ambertemplate, constraint="modredundant", xyz=None, **kwargs):
        """ 
        we reuse most of the code for ONIOM type calculations. 

        since FF calculations do not distinguish between unrestricted and restricted,
        the runtype does not affect this calculation.

        """
        cogef_loop.__init__(self,**kwargs)
        self.constraint = constraint
        self.glog = gaussian.CheckAmberLogfile()
        self.ginp = gaussian.AmberInput(mem=self.ginp.link0.mem, nproc = self.ginp.link0.nproc,
             charge_multi = self.ginp.molecule.charge_multi)
        self.ginp.molecule.read_oniom_template(ambertemplate)
        if xyz:
            self.ginp.molecule.read_xyz(xyz)
        if self.constraint == "opt_flag":
            logger.debug("setting opt flags..")
            self.convert_modredundant_to_opt_flags()
            self.modredundant = None
        self.write_ginp = self._write_ginp
        # independent of runtype:
        self.check_error = True
        self.check_stationary = True
        self.check_stability = False
        self.mix_guess = False


    def _write_ginp(self,ginp_filename="",**kwargs):
        """
        generate input file for amber cogef calculations. T
        The constraint used, can be added as a modredundant section or using the optimisation flag.

        Paramter: 
                ginp_filename : str 
                        filename of the ginp file that should be generated.

        """
        with open(ginp_filename, "w") as of:
            self.ginp.of = of
            # JOB 1
            route = gaussian.classGaussianRoute()
            route.route = [self.print_level, self.level_of_theory, "nosymm", "test"] 
            route.iop = ["2/15=3"]
            route.geom = ["connect"]
            route.opt = ["NoMicro"]
            if self.maxcyc: route.opt.append("maxcyc="+str(self.maxcyc))
            if self.constraint == "modredundant": route.opt.append("modredundant")
            self.ginp.write_inputfile(link1=False,route=route, geom=True, connectivity=True, 
                    modredundant=self.modredundant,parm=True)