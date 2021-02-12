import numpy as np

import ase
import ase.units
import ase.calculators.calculator as calc

import cellconstructor as CC 
import cellconstructor.Structure



class ToyModelCalculator(calc.Calculator):

    def __init__(self,  *args, **kwargs):
        """
        
        V(r) = shift + De *(1 - exp(-a*(r - R0)))^2 + E * x

        
        """
        calc.Calculator.__init__(self, *args, **kwargs)

        # The equilibrium bond lenght 
        self.R0 = 1

        # The rigid energy shift (to have 0 at R -> oo)
        self.shift = 0

        # Depth of the potential
        self.De = 1

        # The force constant of the potential
        self.a = 1

        # Crystal field
        self.E = 0


    
    def calculate(self, atoms = None,  *args, **kwargs):
        calc.Calculator.calculate(self, atoms, *args, **kwargs)

        
        coords = atoms.get_positions()

        r_vec =  coords[0,:] - coords[1,:]
        r_mod = np.sqrt(r_vec.dot(r_vec))


        energy = 0
        force = np.zeros((2, 3), dtype = np.double)

        #TODO: get energy and force
        
        self.results = {"energy": energy, "forces": force}
