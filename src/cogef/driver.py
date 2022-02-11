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
    def __init__(self,runtype, atom1, atom2, dx = 0.1, 
                    level_of_theory = "B3LYP/6-31G*" , read_guess = False, cycles = 1, 
                    reverse = None, restart = 0, modredundant=None, symmetric=True, 
                    dp=1.0, fragment=[], max_error_cycles = 5, 
                    trajectory = None, checkpoint = None ):
        """
        initalise all important parameters

        Parameters: 
                    runtype : str 
                        choices: "restricted"
                        this will set the formalism and dependend checks / input keywords
                    atom1 : int 
                    atom2 : int 
                        atoms to stretch starting at 0
                    dx : float
                        stretch increment
                    level_of_theory : str
                    read_guess : boolean
                        if True, use previous wfn as guess from chk file.
                    cycles : int 
                        cycles to run 
                    reverse : Boolean
                        if True, run in other direction usind abs(dx)
                    restart : int
                        restart a previuous calculation at cycle
                    modredunant : list
                        list of additional modredundant options. elements e.g. "1 2 F"
                    symmetric : boolean 
                        stretch symmetric (True) or asymmetric (False)
                    dp : float
                        dp used to scale the stretching of the fragment atoms
                    fragment : list 
                        contains all the atoms of the first fragment to be stretched.
                    max_error_cycles : int
                        set the maximum number of errors to restart each cycle.
                    trajectory : of,None
                    checkpoint : of,None
                        write trajectory / checkpoint file if points to an open file

        internal Parameters:
                    check_stability : boolean   
                    check_error : boolean
                    check_stationary : boolean
                    check for stability / error termination / stationary 
                            point if True
        
        """
        self.level_of_theory = level_of_theory.strip()
        self.read_guess = read_guess
        self.dx = dx 
        self.atom1 = atom1
        self.atom2 = atom2
        self.modredundant = ["{} {} F".format(atom1, atom2)]
        self.symmetric = symmetric
        self.dp = dp
        self.fragment = fragment
        self.max_error_cycles = max_error_cycles
        self.glog = gaussian.CheckGaussianLogfile()
        self.ginp = gaussian.GaussianInput()
        self.trajectory = trajectory
        self.checkpoint = checkpoint
        if runtype == "restricted":
            self.write_ginp = self._write_ginp_restricted # here we are not calling the function!
            self.check_error = True
            self.check_stationary = True
            self.check_stability = False
        elif runtype == "unrestricted":
            self.write_ginp = self._write_ginp_unrestricted # here we are not calling the function!
            self.check_error = True
            self.check_stationary = True
            self.check_stability = True
        if modredundant:
            self.modredundant.extend(modredundant)
        # loop_range for cogef main loop
        if reverse:
            # we loop from reverse n cycles backwards. smallest point is 000
            self.loop_range = range(reverse, max(-1, reverse - cycles -1 ), -1)
        else:
            self.loop_range = range(restart, cycles +1 )

    def run(self):
        """

        """
        # start the outer cogef loop
        for cycle in self.loop_range:
            logger.info("Starting cycle {}".format(cycle))
            basename = "cogef_{:03d}".format(int(cycle))
            glog_filename = basename + ".log"
            ginp_filename = basename + ".com"        
            error_cycle = 0
            # start the inner loop 
            while error_cycle < self.max_error_cycles:
                # write the gaussian input file. the actual function will be set at init
                if self.read_guess and cycle >= 1: 
                    self.write_ginp(ginp_filename, read_guess=True)
                else:
                    self.write_ginp(ginp_filename, read_guess=False)
                self.run_gauss(ginp_filename)
                # checking the logfile for errors and updated coordinates:
                new_coordinates, glog_is_ok = self.check_gaussian_logfile(glog_filename,cycle,error_cycle)
                self.ginp.molecule.coordinates = new_coordinates
                if glog_is_ok:
                    if self.trajectory:
                        logger.info(f"Adding structure to trajectory {self.trajectory}")
                        self.ginp.molecule.write_xyz(self.trajectory, comment=self.glog.comment_line(point=cycle))
                    break
                else:
                    os.rename(glog_filename, "{}_error_{}.log".format(basename, error_cycle))
                    error_cycle += 1
            # Modify coords for next cycle
            self.ginp.molecule.coordinates = modstruct.mod_fragments(self.ginp.molecule.coordinates, 
                               self.atom1, self.atom2, self.dx, symmetric=self.symmetric, dp=self.dp, fragment=self.fragment)
            # write checkpoint structure to disk.
            if self.checkpoint:
                logger.info(f"Writing checkpoint xyz file {self.checkpoint}")
                self.ginp.molecule.write_xyz(self.checkpoint,comment=f"input structure for cycle {cycle+1}")

    def _write_ginp_restricted(self, ginp_filename="", read_guess=False):
        """
        generate input file for restricted cogef calculations. This will NOT CHECK the stability of the 
        wfn at each cycle.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, use the previous wfn as guess.
        """
        if self.read_guess:
            self.ginp.link0.chk="guess.chk"
            logger.info(f"Adding structure to trajectory {self.trajectory}")
        with open(ginp_filename, "w") as of:
            self.ginp.of = of
            route = "#P " + self.level_of_theory + " nosymm OPT=(modredundant,RFO) "
            if read_guess:
                route += "guess=read"
            self.ginp.write_inputfile(link1=False,route=route,geom=True,modredundant=self.modredundant)

    def _write_ginp_unrestricted(self, ginp_filename="", read_guess=False):
        """
        generate input file for unrestricted cogef calculations. This will also check the stability of the 
        wfn at each cycle and restart the optimisation if the wfn is not stable.

        Paramter: 
                ginp_filename : str
                        filename of the ginp file that should be generated.
                read_guess : boolean
                        if True, use the previous wfn as guess.
        """
        logger.info(f"setting chk file: guess.chk")   
        self.ginp.link0.chk="guess.chk"    
        with open(ginp_filename, "w") as of:
            self.ginp.of = of
            # JOB 1 geometry optimisation
            route = "#P " + self.level_of_theory + " nosymm OPT=(modredundant,RFO) "
            if read_guess:
                route += "guess=read"
            self.ginp.write_inputfile(link1=False,route=route,geom=True,modredundant=self.modredundant)
            # JOB 2 stability analysis
            route = "#P " + self.level_of_theory + " nosymm stable=opt guess=read geom=allcheck"
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
        self.glog.read_log()
        self.glog.check_bonding()
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
        

