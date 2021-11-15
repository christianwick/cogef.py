#!/usr/bin/env python3
""" analyse_cogef.py - read xyz files with energy line as comment,
     compute distance between atoms
     plot energies
     write energy files and relative energy files for processing.
"""

import sys
import numpy as np
import csv


# read command line args
# atom1 atom2: use distance between those atoms to compute 
# forces, distances and plots. gaussian convention: first atom = 1 
filename = sys.argv[1] 
atom1 = int(sys.argv[2])-1
atom2 = int(sys.argv[3])-1

# define constants
hartree_to_kcal_mol = 627.5095
hartree_to_kJ_mol = 2625.498 # https://www.cup.uni-muenchen.de/oc/zipse/teaching/computational-chemistry-1/constants-and-conversion-units/
hartree_to_J = 4.3597447222E-18  # wikipedia.org
angstrom_to_m = 1.0E-10 
#kcal_molA_to_nN = 1.66053892103219E-11 * (4.184**-1) * 1E9

def write_csv(filename,fieldnames,data):
    with open (filename, "w") as csvfile:
        writer = csv.DictWriter(csvfile,fieldnames = fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in data:
            writer.writerow(row)

# read in xyz file
comments = []
molecules = []
with open(filename,"r") as fil:
    while True:
        try:
            natom = int(fil.readline())
        except:
            break
        else:
            comments.append(fil.readline())
            itern=0
            atoms=[]
            while itern < natom:
                atoms.append(fil.readline())
                itern+=1
            molecules.append(atoms)


# compute distances
distances = []
for atoms in molecules:
    temp1 = atoms[atom1].split()
    temp2 = atoms[atom2].split()
    r1 = np.array([float(temp1[1]),float(temp1[2]),float(temp1[3])])
    r2 = np.array([float(temp2[1]),float(temp2[2]),float(temp2[3])])
    dist = np.linalg.norm(r1-r2)
    distances.append(dist)

# try to read energy
energies = []
for line in comments:
    temp = line.split()
    energies.append(float(temp[2]))

# compute relative energies
min_en = np.min(energies)
equi_dist = distances[np.argmin(energies)]
data = []
for ii in range(len(distances)):
    temp = {}
    temp["Structure"] = ii
    temp["Distance [angstrom]"] = distances[ii]
    temp["Relative Distance [angstrom]"] = distances[ii] - equi_dist
    temp["Absolute Energy [Hartree]"] = energies[ii]
    temp["Relative Energy [Hartree]"] = energies[ii] - min_en
    temp["Relative Energy [kJ/mol]"] = ( energies[ii] - min_en ) * hartree_to_kJ_mol
    temp["Relative Energy [kcal/mol]"] = ( energies[ii] - min_en ) * hartree_to_kcal_mol
    data.append(temp)
rel_en_kj_mol = [ x["Relative Energy [kJ/mol]"] for x in data ]
# compute forces in nN
forces = []
de_array = []
dx_array = []
for ii in range(0,len(distances)-1):
    de = (energies[ii+1] - energies[ii]) * hartree_to_J
    dx = (distances[ii+1] - distances[ii]) * angstrom_to_m
    forces.append(de/dx * 1.0E9 )  # in nN
    de_array.append(de)  # in J
    dx_array.append(dx * 1.0E9 )  # in nm 
forces = np.array(forces) 
# print some info
print("File: {}".format(filename))
print("Atoms {} {}".format(atom1+1,atom2+1))
print("re (using inital distance) [A]: {:10.4f}  ".format(distances[0]))
print("De [kJ/mol] (maximum)         : {:10.4f}  ".format(np.max(rel_en_kj_mol)))
print("De [kJ/mol] (last point)      : {:10.4f}  ".format(rel_en_kj_mol[-1]))
print("Max F [nN]                    : {:10.4f} ".format(np.max(forces)))
print("corresponding distance [A]    : {:10.3f} ".format(distances[1+np.argmax(forces)]))
print("corresponding structure       : {:3d} ".format(1+np.argmax(forces)))
# write force data into array of dictionaries.
force_data = []
for ii in range(len(forces)):
    temp = {}
    temp["Structure"] = ii + 1
    temp["Force [nN]"] = forces[ii] 
    temp["de [J]"] = de_array[ii]
    temp["dx [nm]"] = dx_array[ii] 
    force_data.append(temp)

# compute Veff
rel_dist_angstrom = np.array([ x["Relative Distance [angstrom]"] for x in data ])
rel_en_kj_mol = np.array(rel_en_kj_mol)
work_kJ_mol = rel_dist_angstrom * 1.0E-10 *1.0E-9 * 1.0E-3 * 6.022E23 
Veff_1_kj_mol = rel_en_kj_mol - work_kJ_mol 
Veff_2_kj_mol = rel_en_kj_mol - work_kJ_mol * 2.0
Veff_3_kj_mol = rel_en_kj_mol - work_kJ_mol * 3.0
Veff_4_kj_mol = rel_en_kj_mol - work_kJ_mol * 4.0
Veff_5_kj_mol = rel_en_kj_mol - work_kJ_mol * 5.0
Veff_1_kj_mol -= np.min(Veff_1_kj_mol[:int(len(Veff_1_kj_mol)/2)]) 
Veff_2_kj_mol -= np.min(Veff_2_kj_mol[:int(len(Veff_1_kj_mol)/2)]) 
Veff_3_kj_mol -= np.min(Veff_3_kj_mol[:int(len(Veff_1_kj_mol)/2)]) 
Veff_4_kj_mol -= np.min(Veff_4_kj_mol[:int(len(Veff_1_kj_mol)/2)]) 
Veff_5_kj_mol -= np.min(Veff_5_kj_mol[:int(len(Veff_1_kj_mol)/2)]) 
for ii in range(len(distances)):
    data[ii]["Veff 1.0 [kJ/mol]"] = Veff_1_kj_mol[ii]   
    data[ii]["Veff 2.0 [kJ/mol]"] = Veff_2_kj_mol[ii]   
    data[ii]["Veff 3.0 [kJ/mol]"] = Veff_3_kj_mol[ii]   
    data[ii]["Veff 4.0 [kJ/mol]"] = Veff_4_kj_mol[ii]   
    data[ii]["Veff 5.0 [kJ/mol]"] = Veff_5_kj_mol[ii]   


fieldnames = ["Structure","Distance [angstrom]","Relative Distance [angstrom]", "Absolute Energy [Hartree]", "Relative Energy [Hartree]", "Relative Energy [kJ/mol]",
    "Relative Energy [kcal/mol]"]
write_csv(filename+"_{}_{}".format(atom1+1,atom2+1)+"_en.csv",fieldnames,data)
fieldnames =["Structure","dx [nm]","de [J]","Force [nN]"]
write_csv(filename+"_{}_{}".format(atom1+1,atom2+1)+"_f.csv",fieldnames,force_data)
fieldnames = ["Structure","Distance [angstrom]","Relative Distance [angstrom]","Relative Energy [kJ/mol]",
    "Veff 1.0 [kJ/mol]","Veff 2.0 [kJ/mol]","Veff 3.0 [kJ/mol]","Veff 4.0 [kJ/mol]","Veff 5.0 [kJ/mol]"]
write_csv(filename+"_{}_{}".format(atom1+1,atom2+1)+"_veff.csv",fieldnames,data)