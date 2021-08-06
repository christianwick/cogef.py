""" modify structures
 
"""


import numpy as np


def mod_single_atom(coords,atom1,atom2,dx=0.02,symmetric=False):
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
        new_r2 = r2 + vec * dx
        new_coords[atom2] = [new_r2[0], new_r2[1], new_r2[2]]
    if len(coords[0]) == 4:
        return(zip(elements,new_coords))
    else:
        return(new_coords)

def mod_atoms(coords,atom1,atom2,dx=0.02,symmetric=False):
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    mat = np.zeros([len(coords),3])
    for ii,atom in enumerate(coords):
        mat[ii,:] = [atom[1],atom[2],atom[3]]
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
                print("move atom {} in + direction.".format(ii+1) )
            else:
                r = mat[ii,:] - vec * dx * 0.5
                print("move atom {} in - direction.".format(ii+1) )
            new_coords[ii] = [coords[ii][0],r[0],r[1],r[2]]
    else:    
        for ii in range(len(dist_mat_r1)):
            if dist_mat_r2[ii] < dist_mat_r1[ii] :
                r = mat[ii,:] + vec * dx
                print("move atom {}".format(ii+1) )
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
    if len(coords[0]) == 4:
        elements = [x[0] for x in coords]
        coords = [[x[1],x[2],x[3]] for x in coords ]
    mat = np.zeros([len(coords),3])
    for ii,atom in enumerate(coords):
        mat[ii,:] = [atom[0],atom[1],atom[2]]
    r1 = mat[atom1,:]
    r2 = mat[atom2,:]
    vec = (r2 - r1) / np.linalg.norm(r2 - r1)

    new_coords = coords 
    if symmetric:
        for ii in range(len(coords)):
            if ii in fragment:
                print("move atom {} in - direction.".format(ii+1) )
                if ii != atom1 and ii != atom2:
                    r = mat[ii,:] - vec * dx * dp * 0.5
                else:
                    r = mat[ii,:] - vec * dx * 0.5
            else:
                print("move atom {} in + direction.".format(ii+1) )
                if ii != atom1 and ii != atom2:
                    r = mat[ii,:] + vec * dx * dp * 0.5
                else:
                    r = mat[ii,:] + vec * dx * 0.5 
            new_coords[ii] = [r[0],r[1],r[2]]
    else:    
        for ii in range(len(coords)):
            if ii in fragment:
                print("move atom {}".format(ii+1) )
                if ii != atom1 and ii != atom2:
                    r = mat[ii,:] - vec * dx * dp
                else:
                    r = mat[ii,:] - vec * dx
                new_coords[ii] = [r[0],r[1],r[2]]
    if len(coords[0]) == 4:
        return(zip(elements,new_coords))
    else:
        return(new_coords)