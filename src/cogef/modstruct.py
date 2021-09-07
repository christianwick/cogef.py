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

def mod_atoms(coords,atom1,atom2,dx=0.02,symmetric=False):
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    mat = np.array(coords)
    r1 = mat[atom1,:]
    r2 = mat[atom2,:]
    vec = (r2 - r1) / np.linalg.norm(r2 - r1)
    dist_mat_r1 = np.linalg.norm(mat - r1, axis=1)
    dist_mat_r2 = np.linalg.norm(mat - r2, axis=1)
    print(dist_mat_r1,dist_mat_r2)

    new_coords = coords 
    if symmetric:
        for ii in range(len(dist_mat_r1)):
            if dist_mat_r2[ii] < dist_mat_r1[ii] :
                r = mat[ii,:] + vec * dx * 0.5
                print("move atom {} in + direction.".format(ii) )
            else:
                r = mat[ii,:] - vec * dx * 0.5
                print("move atom {} in - direction.".format(ii) )
            new_coords[ii] = [coords[ii][0],r[0],r[1],r[2]]
    else:    
        for ii in range(len(dist_mat_r1)):
            if dist_mat_r2[ii] < dist_mat_r1[ii] :
                r = mat[ii,:] + vec * dx
                print("move atom {}".format(ii) )
                new_coords[ii] = [coords[ii][0],r[0],r[1],r[2]]
    if len(coords[0]) == 4:
        return(zip(elements,new_coords))
    else:
        return(new_coords)

def mod_fragments(coords, atom1, atom2, dx=0.02, dp=1.0, fragment=[0,1,2,3], symmetric=False):
    """ 
        coords: list of coordinates [x,y,z],...]
        atom1: first atom to stretch
        atom2: second atom to stretch
        dx: amount of stretching
        fragment: needs to contain atom1
        dp: add an additional perturbation for atoms exept atom1 and atom2
    """
    logger.debug("moving fragments: {} {} by {} A".format(atom1,atom2,dx))
    logger.debug("fragment atoms: {} ; dp: {} A".format(fragment, dp))
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    mat = np.array(coords)
    r1 = mat[atom1,:]
    r2 = mat[atom2,:]
    vec = (r2 - r1) / np.linalg.norm(r2 - r1)
    vec_dx = vec * dx
    vec_dp = vec_dx * dp
    u = np.zeros_like(mat)
    if symmetric:
        vec_dx = vec_dx * 0.5
        vec_dp = vec_dp * 0.5 
        for ii in range(len(coords)):
            if ii == atom1:
                logger.debug("move anchor {} in - direction.".format(ii) )
                u[ii,:] -=  vec_dx
            elif ii == atom2:
                logger.debug("move anchor {} in + direction.".format(ii) )
                u[ii,:] +=  vec_dx
            elif ii in fragment:
                logger.debug("move atom {} in - direction.".format(ii) )
                u[ii,:] -=  vec_dp
            else:
                logger.debug("move atom {} in + direction.".format(ii) )
                u[ii,:] +=  vec_dp
    else:    
        for ii in range(len(coords)):
            if ii == atom1:
                logger.debug("move anchor {} in - direction.".format(ii) )
                u[ii,:] -=  vec_dx
            elif ii in fragment:
                logger.debug("move atom {} in - direction.".format(ii) )
                u[ii,:] -=  vec_dp
    mat += u
    if len(coords[0]) == 4:
        return(zip(elements,mat))
    else:
        return(mat)