""" fragments.py contains logic to parse fragment atom specifications 

fragments can be specified as strings: e.g. "1,2" or "1-4" or "1-4,5-6,7"
this needs to be converted to a list of atoms.

"""

def convert_inputstr(inpstr=None):
    """ parse input string and produce fragment list 
    
    Parameters:
        inpstr : str
            inputstring e.g. "1,2,4-5,6-7,12"
    
    Returns: 
        fragment : list ( int )

    """
    fragment = []
    if inpstr != "":
        for val in inpstr.split(","):
            temp = val.split("-",maxsplit=1)
            if len(temp) == 2:
                fragment.extend([ x - 1  for x in range( int(temp[0]), int(temp[1])+1) ])
            elif len(temp) == 1:
                fragment.append(int(temp[0])-1)
    return( sorted(set(fragment)) )
                
if __name__ == "__main__":
    print("--------------------")
    print("test fragments.py")
    frag = "1,2,3,3,4,4,4,5-6,5-8,12-15"
    print(convert_inputstr(frag))
    frag = ""    
    print(convert_inputstr(frag))

    