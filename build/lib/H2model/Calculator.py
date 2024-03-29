import numpy as np

import ase
import ase.units
import ase.calculators.calculator as calc

import cellconstructor as CC 
import cellconstructor.Structure
import cellconstructor.Units as units


class ToyModelCalculator(calc.Calculator):

    def __init__(self,  *args, **kwargs):
        """
        Computes the potential and forces for an H2 molecule:
        
            V(r) = H2_shift + H2_D * {1 - exp[- H2_a * (r - H2_re)]}^2 + E * (x_1 - x_2),
        
            where r = |r_1 - r_2|.
        """
        calc.Calculator.__init__(self, *args, **kwargs)
        
        # Setup what properties the calculator can load
        self.implemented_properties = ["energy", "forces"]
        
        # Chose the model
        self.model = 'rotating'
        
        # Chose if we want a com fixed or not
        self.fix_com = False
        
        # Parameters of the Morse potential in ATOMIC UNITS
        
        # The rigid energy shift in HARTREE
        self.H2_shift = -1.17225263
        
        # Depth of the potential in HARTREE
        self.H2_D     =  0.13681332
        
        # The force constant of the potential in  1/BOHR
        self.H2_a     =  1.21606669
        
        # The equilibrium bond lenght  in BOHR
        self.H2_re    =  1.21606669
        
        # The harmnonic elastic constant in HARTREE /BOHR^2
        self.k_harm = 2 * (self.H2_a**2) * self.H2_D
        
        # Crystal field in HARTREE /BOHR
        self.E = 0.001
        
        # The spring constant for the com
        self.k_com = 1. * self.k_harm
        
    def minimum(self):
        """
        Get the minimum of the full potential since this can be done analitycally

        Returns:
        -------
            -Vmin: double, the minimum of the potential in HARTREE
        """

        expmin = 0.5 * (1. + np.sqrt(1 + 2. * self.E * np.cos(np.pi)/(self.H2_D * self.H2_a)))

        rmin = self.H2_re - np.log(expmin) / self.H2_a

        Vmin = self.E * rmin * np.cos(np.pi) + self.H2_shift + self.H2_D * (1 - expmin)**2
    
        return Vmin
    
    
    def calculate(self, atoms = None,  *args, **kwargs):
        """
        COMPUTES ENERGY AND FORCES IN eV and eV/ ANGSTROM
        =================================================
        
        Returns:
        -------
            -self.results: a dict with energy and forces.
        """
        calc.Calculator.calculate(self, atoms, *args, **kwargs)
        
        # Energy and force in HARTREE and HARTREE/BOHR
        energy, force = 0., np.zeros((2, 3), dtype = np.double)

        # Position in ANGSTROM converted in BOHR, np.array with shape = (2, 3)
        coords = atoms.get_positions() * units.A_TO_BOHR
        
        # Get the relative coordinate
        rel_coord =  (coords[0,:] - coords[1,:])
        
        # Get the com cooridnate (same masses)
        com_coord = 0.5 *  (coords[0,:] + coords[1,:])
        
        if self.model == 'rotating':
            # Get the radial distance
            r = np.sqrt(rel_coord.dot(rel_coord))

            # Get the energy in HARTREE subtrating the minimum of the Morse + crystal field potential
            energy = self.H2_shift + self.H2_D * (1. - np.exp(-self.H2_a * (r - self.H2_re)))**2 + self.E * rel_coord[0] - self.minimum()
            
            # Derivative with respect the radial distance
            diff_V_r = 2. * self.H2_a * self.H2_D * (1. - np.exp(-self.H2_a * (r - self.H2_re))) * np.exp(-self.H2_a * (r - self.H2_re))

            # Get the forces for the first particle in HARTREE /BOHR
            force[0,:]  = - diff_V_r * rel_coord[:] /r
            force[0,0] += - self.E

            # Get the forces for the second particle in HARTREE /BOHR
            force[1,:] = - force[0,:]
            
        else:
            # Get the energy in HARTREE 
            energy = 0.5 * self.k_harm * ((rel_coord[0] + self.H2_re)**2 + rel_coord[1]**2 + rel_coord[2]**2)

            # Get the forces for the first particle in HARTREE /BOHR
            force[0,0] = - self.k_harm * (rel_coord[0] + self.H2_re)
            force[0,1] = - self.k_harm * rel_coord[1]
            force[0,2] = - self.k_harm * rel_coord[2]
            
            # Get the forces for the second particle in HARTREE /BOHR
            force[1,:] = - force[0,:]
            
                    
        if self.fix_com:
            energy     += 0.5 * self.k_com * com_coord.dot(com_coord)
            force[0,:] += - 0.5 * self.k_com * com_coord[:]
            force[1,:] += - 0.5 * self.k_com * com_coord[:]
            

        # CONVERT from HARTREE, HARTREE /BOHR in -> eV, eV /ANGSTROM
        self.results = {"energy": energy * 2. * units.RY_TO_EV, "forces": force * 2. * units.RY_TO_EV /units.BOHR_TO_ANGSTROM}
            
    
        return self.results
