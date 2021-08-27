""" constants used to convert units """
__all__ = ["hartree_to_kcal_mol","hartree_to_kJ_mol","hartree_to_J","angstrom_to_m"]

# define constants
hartree_to_kcal_mol = 627.5095
hartree_to_kJ_mol = 2625.498 # https://www.cup.uni-muenchen.de/oc/zipse/teaching/computational-chemistry-1/constants-and-conversion-units/
hartree_to_J = 4.3597447222E-18  # wikipedia.org
angstrom_to_m = 1.0E-10


# PSE
PSE=["","H" ,                                                                                                                                                      "He",
        "Li","Be",                                                                                                                        "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
        "Na","Mg",                                                                                                                        "Al","Si","P" ,"S" ,"Cl","Ar",
        "K" ,"Ca",                                                                      "Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr",                                                                      "Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn"]
