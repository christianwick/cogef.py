""" modify structures
 
"""


import numpy as np
import logging

logger = logging.getLogger("modsruct")

def mod_single_atom(coords,atom1,atom2,dx=0.02,symmetric=False):
    logger.debug("moving single atom: {} {} by {} A".format(atom1,atom2,dx))
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    r1 = np.array([float(coords[atom1][0]),float(coords[atom1][1]),float(coords[atom1][2])])
    r2 = np.array([float(coords[atom2][0]),float(coords[atom2][1]),float(coords[atom2][2])])
    vec = ( r2 - r1 ) / np.linalg.norm( r2 - r1 )
    new_coorbeta = coords 
    if symmetric:
        new_r1 = r1 - vec * dx * 0.5 
        new_r2 = r2 + vec * dx * 0.5
        new_coorbeta[atom1] = [new_r1[0], new_r1[1], new_r1[2]]
        new_coorbeta[atom2] = [new_r2[0], new_r2[1], new_r2[2]]
    else:
        new_r1 = r1 - vec * dx
        new_coorbeta[atom1] = new_r1
    if len(coords[0]) == 4:
        return(zip(elements,new_coorbeta))
    else:
        return(new_coorbeta)


def mod_fragments(coords, atom1, atom2, dx=0.02, alpha=1.0, frag1=None, frag2=None, exclude = None,
        symmetric=False, current_strain = 0.0, beta = 1.0, nbonds=1.0, cguess=None):
    """ modify the coordinates of molecules using fragments

    compute the vector between atom1 and atom2 and stretch all atoms by a certain amount
        vec = atom2 - atom1

        Symmetric = True:
        atom1 stretched by -dx ; atom2 stretched by +dx
        frag a atoms: stretched by -dx * alpha + strain * beta

    Parameter:
        coords: list 
            lsit of coordinates [x,y,z],...]
        atom1: int 
            first atom to stretch by dx
        atom2: int
            second atom to stretch by dx
        dx: float
        frag1: list
            litst of atom numbers starting at 0
            atom1 must be part of fragments!
        frag2: list
            list of atom numbers starting at 0
            atom2 must be part of fragments!
        exclude: list
            list of atom numbers starting at 0 to be excluded from movement.
        alpha: float
            add an additional perturbation for atoms exept atom1 and atom2
            no perturbation for alpha=1.0.
        current_strain : float
            add an additive perturbation to separate fragments
        beta: float
            add an additional multiplicative damping to the strain perturbation
        nbonds: float
            number of bonds along the chain. The lenght of the translation vector 
            for the fragment atoms will be devided by this value if the strain method 
            is used 

    returns
        mat : np.array(n,3)
    """
    logger.debug("moving atoms: {} {} by {} A, symmetric = {}".format(atom1,atom2,dx,symmetric))
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    frag_2 = frag2
    frag_1 = frag1
    #if fragment != []: 
    #    frag_2 = fragment[:]
    #    frag_1 = [ x for x in range(len(coords)) if x not in frag_2 ]
    if exclude:
        frag_2 = [ x for x in frag_2 if x not in exclude ]
        frag_1 = [ x for x in frag_1 if x not in exclude ]
    if atom2 in frag2: frag_2.remove(atom2)
    if atom1 in frag1: frag_1.remove(atom1)
    if atom1 in frag2 or atom2 in frag1:
        logger.warning("Wrong atom definitions. Found stretching atom in wrong fragment. CHECK YOUR INPUT.")
    logger.debug("fragment 2: {} ; alpha: {:4.3f} A; strain: {:4.3f}; beta: {:4.3f}".format(frag_2, alpha, current_strain, beta))
    logger.debug("fragment 1: {} ; alpha: {:4.3f} A; strain: {:4.3f}; beta: {:4.3f}".format(frag_1, alpha, current_strain, beta))
    mat = np.array(coords)
    vec = ( mat[atom2,:] - mat[atom1,:] ) / np.linalg.norm( mat[atom1,:] - mat[atom2,:] )
    vec_dx = vec * dx
    logger.debug(f"effective dx: {np.linalg.norm(vec_dx):4.3f}")

    if cguess == "dx":
        u = _guess_dx(mat=mat,vec=vec,vec_dx=vec_dx,dx=dx,frag_1=frag_1,frag_2=frag_2,atom1=atom1,
            atom2=atom2,alpha=alpha,symmetric=symmetric)
    elif cguess == "strain":
        u = _guess_using_strain(mat=mat,vec=vec,vec_dx=vec_dx,frag_1=frag_1,frag_2=frag_2,atom1=atom1,atom2=atom2,
            current_strain=current_strain,beta=beta,nbonds=nbonds,symmetric=symmetric)    
    elif cguess == "new_strain":
        u = _guess_using_new_strain(mat=mat,vec=vec,vec_dx=vec_dx,dx=dx,frag_1=frag_1,frag_2=frag_2,
            atom1=atom1,atom2=atom2,beta=beta,symmetric=symmetric)
    else:
        u = _guess_dx(mat=mat,vec=vec,vec_dx=vec_dx,dx=dx,frag_1=frag_1,frag_2=frag_2,atom1=atom1,
            atom2=atom2,alpha=alpha,symmetric=symmetric)

    mat += u

    if len(coords[0]) == 4:
        return(zip(elements,mat))
    else:
        return(mat)



def _guess_using_strain(mat,vec,vec_dx,frag_1,frag_2,atom1,atom2,current_strain,beta,nbonds,symmetric):
    u = np.zeros_like(mat)
    if symmetric:
        vec_2 = -1 * ( mat[frag_2] - mat[atom2] ) * current_strain * 0.5 / nbonds
        vec_1 = -1 * ( mat[frag_1] - mat[atom1] ) * current_strain * 0.5 / nbonds
        logger.debug(f"norm vec_2: {np.linalg.norm(vec_2,axis=1)}")
        logger.debug(f"norm vec_1: {np.linalg.norm(vec_1,axis=1)}")
        dot_2 = np.abs(np.matmul(vec_2,vec))
        dot_1 = np.abs(np.matmul(vec_1,vec))
        u[frag_2] += vec  * ( beta * dot_2[:,None] )
        u[frag_1] -= vec  * ( beta * dot_1[:,None] )
        u[atom1] -= vec_dx * 0.5
        u[atom2] += vec_dx * 0.5
        
    else: 
        vec_2 = -1 * ( mat[frag_2] - mat[atom2] ) * current_strain / nbonds
        logger.debug(f"norm vec_2: {np.linalg.norm(vec_2,axis=1)}")
        dot_2 = np.abs(np.matmul(vec_2,vec))
        u[frag_2] += vec * ( beta * dot_2[:,None] )
        u[atom2] += vec_dx 
    logger.debug(f"norm u: {np.linalg.norm(u,axis=1)}")
    return(u)

def _guess_dx(mat,vec,vec_dx,dx,frag_1,frag_2,atom1,atom2,alpha,symmetric):
    u = np.zeros_like(mat)
    if symmetric:
        u[frag_2] += vec  * 0.5 * ( dx * alpha ) 
        u[frag_1] -= vec  * 0.5 * ( dx * alpha ) 
        u[atom1] -= vec_dx * 0.5
        u[atom2] += vec_dx * 0.5
        
    else: 
        u[frag_2] += vec * ( dx * alpha )
        u[atom2] += vec_dx 
    logger.debug(f"norm u: {np.linalg.norm(u,axis=1)}")
    return(u)

def _guess_using_new_strain(mat,vec,vec_dx,dx,frag_1,frag_2,atom1,atom2,beta,symmetric):
    u = np.zeros_like(mat)
    ac = np.linalg.norm( mat[atom1,:] - mat[atom2,:] )
    if symmetric:
        factor = 0.5 * dx / (ac + 0.5 * dx)
        vec_2 = -1 * ( mat[frag_2] - ( mat[atom2] + vec_dx * 0.5 )) * factor
        vec_1 = -1 * ( mat[frag_1] - ( mat[atom1] - vec_dx * 0.5 )) * factor
        logger.debug(f"norm vec_2: {np.linalg.norm(vec_2,axis=1)}")
        logger.debug(f"norm vec_1: {np.linalg.norm(vec_1,axis=1)}")
        dot_2 = np.abs(np.matmul(vec_2,vec))
        dot_1 = np.abs(np.matmul(vec_1,vec))
        u[frag_2] += vec  * ( 0.5 * dx + beta * dot_2[:,None] )
        u[frag_1] -= vec  * ( 0.5 * dx + beta * dot_1[:,None] )
        u[atom1] -= vec_dx * 0.5
        u[atom2] += vec_dx * 0.5
        
    else: 
        factor = 0.5 * dx / (ac + 0.5 * dx)
        vec_2 = -1 * ( mat[frag_2] - ( mat[atom2] + vec_dx * 0.5 )) * factor
        logger.debug(f"norm vec_2: {np.linalg.norm(vec_2,axis=1)}")
        dot_2 = np.abs(np.matmul(vec_2,vec))
        u[frag_2] += vec * ( dx + beta * dot_2[:,None] )
        u[atom2] += vec_dx 
    logger.debug(f"norm u: {np.linalg.norm(u,axis=1)}")
    return(u)

