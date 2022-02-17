""" mulliken.py - read mulliken spin density and charges from qm log files """


import logging

logger = logging.getLogger("mulliken")

class Mulliken():
    def __init__(self):
        """ The Mullikan base class.

        Internal Parameters:
            self.atom_name = list (atom name associated with)
            self.atom_number = list (atom number)
            self.charge = list (mulliken charges)
            self.spin = list (mulliken spin density)

        """
        logger.debug("Initialise new Mulliken class")
        self.atom_name = []
        self.atom_number = []
        self.charge = []
        self.spin = [] 

    def __str__(self):
        mystr = ""
        for nn, atom_num in enumerate(self.atom_number):
            mystr += f"{self.atom_number[nn]} {self.atom_name[nn]} {self.charge[nn]} {self.spin[nn]} \n"
        return(mystr)

    def _read_gaussian_raw_data(self,raw_data):
        """
        convert lines of type :
        2  C   -0.014288  -0.000000
        atom number, atom name, charge, spin

        Input: 
            raw_data = list (str)
        """
        atom_number, atom_name, charge, spin = [], [], [], []
        for line in raw_data:
            temp = line.split()
            atom_number.append(int(temp[0]))
            atom_name.append(str(temp[1]))
            charge.append(float(temp[2]))
            spin.append(float(temp[3]))
        self.atom_number = atom_number
        self.atom_name = atom_name 
        self.charge = charge 
        self.spin = spin
        
    def sort_by(self, key="spin"):
        """
        sort data in class Mulliken by axis

        Input: key = str  [spin,charge,atom_name,atom_number]
        """
        logger.debug(f"sorting Mulliken data by {key}")
        all_data = zip(self.atom_number,self.atom_name,self.charge,self.spin)
        values = {
            "atom_number" : 0,
            "atom_name" : 1,
            "charge" : 2,
            "spin" : 3 
            }
        sorted_data = sorted(all_data, key= lambda val: val[values[key]])
        self.atom_number, self.atom_name, self.charge, self.spin = [list(l) for l in zip(*sorted_data) ]

    def find_spin_larger_than(self, thresh = 0.5 ):
        temp = []
        for nn, sp in enumerate(self.spin):
            if abs(sp) >= thresh:
                temp.append(nn)
        temp_mul = Mulliken()
        temp_mul.atom_number = [ self.atom_number[i] for i in temp ]
        temp_mul.atom_name = [ self.atom_name[i] for i in temp ]
        temp_mul.charge = [ self.charge[i] for i in temp ]
        temp_mul.spin = [ self.spin[i] for i in temp ]
        return ( temp_mul )




            

    
