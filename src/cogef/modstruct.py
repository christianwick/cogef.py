""" modify structures
 
"""


import numpy as np
import logging

logger = logging.getLogger("modstruct")

def mod_single_atom(coords,atom1,atom2,dx=0.02,symmetric=False):
    logger.debug("moving single atom: {} {} by {} A".format(atom1,atom2,dx))
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    r1 = np.array([float(coords[atom1][0]),float(coords[atom1][1]),float(coords[atom1][2])])
    r2 = np.array([float(coords[atom2][0]),float(coords[atom2][1]),float(coords[atom2][2])])
    vec = ( r2 - r1 ) / np.linalg.norm( r2 - r1 )
    new_coords = coords 
    if symmetric:
        new_r1 = r1 - vec * dx * 0.5 
        new_r2 = r2 + vec * dx * 0.5
        new_coords[atom1] = [new_r1[0], new_r1[1], new_r1[2]]
        new_coords[atom2] = [new_r2[0], new_r2[1], new_r2[2]]
    else:
        new_r1 = r1 - vec * dx
        new_coords[atom1] = new_r1
    if len(coords[0]) == 4:
        return(zip(elements,new_coords))
    else:
        return(new_coords)


def mod_fragments(coords, atom1, atom2, dx=0.02, dp=1.0, fragment=[0,1,2,3],
        symmetric=False, current_strain = 0.0, ds = 1.0):
    """ modify the coordinates of molecules using fragments

    compute the vector between atom1 and atom2
        vec = atom2 - atom1

    Parameter:
        coords: list 
            lsit of coordinates [x,y,z],...]
        atom1: int 
            first atom to stretch by dx
        atom2: int
            second atom to stretch by dx
        dx: float
        fragment: list
            litst of atom numbers starting at 0
            atom1 must be part of fragments!
        dp: float
            add an additional perturbation for atoms exept atom1 and atom2
            no perturbation for dp=1.0.
        current_strain : float
            add an additive perturbation to separate fragments
        ds: float
            adds an additional multiplicative damping to the strain perturbation

    
    returns
        mat : np.array(n,3)
    """
    logger.debug("moving atoms: {} {} by {} A, symmetric = {}".format(atom1,atom2,dx,symmetric))
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    frag_a = []
    frag_b = []
    if fragment != []: 
        frag_a = fragment[:]
        frag_b = [ x for x in range(len(coords)) if x not in frag_a ]
        frag_a.remove(atom1)
        frag_b.remove(atom2)
    logger.debug("fragment A: {} ; dp: {:4.3f} A; strain: {:4.3f}; ds: {:4.3f}".format(frag_a, dp, current_strain, ds))
    logger.debug("fragment B: {} ; dp: {:4.3f} A; strain: {:4.3f}; ds: {:4.3f}".format(frag_b, dp, current_strain, ds))
    mat = np.array(coords)
    vec = ( mat[atom2,:] - mat[atom1,:] ) / np.linalg.norm( mat[atom1,:] - mat[atom2,:] )
    vec_dx = vec * dx
    vec_dp = vec * (dx * dp + current_strain * ds )

    u = np.zeros_like(mat)
    if symmetric:
        u[frag_a] -= vec_dp * 0.5
        u[frag_b] += vec_dp * 0.5
        u[atom1] -= vec_dx * 0.5
        u[atom2] += vec_dx * 0.5
    else: 
        u[frag_a] -= vec_dp 
        u[atom1] -= vec_dx 

    mat += u

    if len(coords[0]) == 4:
        return(zip(elements,mat))
    else:
        return(mat)


