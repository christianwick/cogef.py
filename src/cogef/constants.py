""" constants used to convert units """
__all__ = ["hartree_to_kcal_mol","hartree_to_kJ_mol","hartree_to_J","angstrom_to_m"]

# define constants
hartree_to_kcal_mol = 627.5095
hartree_to_kJ_mol = 2625.498 # https://www.cup.uni-muenchen.de/oc/zipse/teaching/computational-chemistry-1/constants-and-conversion-units/
hartree_to_J = 4.3597447222E-18  # wikipedia.org
kcal_mol_to_kJ_mol = 4.184 # https://www.cup.uni-muenchen.de/oc/zipse/teaching/computational-chemistry-1/constants-and-conversion-units/
angstrom_to_m = 1.0E-10
avogadro = 6.022E23 
# 1nNm = 1.0E-9 Nm = 1.0E-09 J = 1.0E-09 * 1.0E-03 kJ = 1.0E-12 kJ
nNm_to_kJ_mol = 1.0E-12 * avogadro 
# 1nN = 1.0E-9 N = 1.0E-09 J/m = 1.0E-09 * 1.0E-03 kJ/m = 1.0E-12 kJ / m 
nN_to_kJ_mol_angstrom = 1.0E-12 * avogadro * angstrom_to_m

# PSE
PSE=["","H" ,                                                                                                                                                      "He",
        "Li","Be",                                                                                                                        "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
        "Na","Mg",                                                                                                                        "Al","Si","P" ,"S" ,"Cl","Ar",
        "K" ,"Ca",                                                                      "Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr",                                                                      "Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"]

# source: www.webelements.com
PSE_MASS = ["",1.008,4.0026,
               6.94,9.0122,10.81,12.011,14.007,15.999,18.998,20.180,
               22.990,24.305,26.982,28.085,30.974,32.06,35.45,39.948,
               39.098,40.078,44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.38,69.723,72.63,74.922,78.96,79.904,83.798,
               85.468,87.62,88.906,91.224,92.906,95.96,97.91,101.07,102.91,106.42,107.87,112.41,114.82,118.71,121.76,127.60,126.90,131.29,
               132.91,137.33,138.91,140.12,140.91,144.24,144.91,150.36,151.96,157.25,158.93,162.50,164.93,167.26,168.93,173.05,174.97,178.49,180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,204.38,207.2,208.98,208.98,209.99,222.02,
               223.02,226.03,227.03,232.04,231.04,238.03,237.05,244.06,243.06,247.07,247.07,251.08,252.08,257.10,258.10,259.10,262.11,265.12,268.13,271.13,270.,277.15,276.15,281.16,280.16,285.17]
               
# Cordero et. al., DOI: 10.1039/b801115j and CCDB recommendation of 1.5 as default
PSE_COVALENT_RADII = ["", 0.31, 0.28,
                      1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
                      1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
                      2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.61, 1.52, 1.50, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
                      2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
                      2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50,
                      2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69 ]
PSE_COVALENT_RADII.extend([1.5 for x in range(len(PSE)-len(PSE_COVALENT_RADII))])