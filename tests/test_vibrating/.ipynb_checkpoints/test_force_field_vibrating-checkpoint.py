import numpy as np
import time
import os
import shutil
import cellconstructor as CC

import H2model 
import H2model.Calculator

import matplotlib.pyplot as plt

import nonlinear_sscha
import nonlinear_sscha.NonLinearStructure as NLS
import nonlinear_sscha.NonLinearEnsemble as NLE
import nonlinear_sscha.Conversion as conv
import nonlinear_sscha.Extra as extra

import shutil

import pytest

def test_forces():
    # Get the current directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Set the force field
    ff_calculator = H2model.Calculator.ToyModelCalculator()
    # ff_calculator.E = 0.01
    ff_calculator.fix_com = False
    ff_calculator.model = 'vibrating'
    re = ff_calculator.H2_re * conv.AU_TO_ANGSTROM
    k_harm =  (2. * ff_calculator.H2_a**2 * ff_calculator.H2_D) * conv.HA_TO_RY

    # Set the initial guess for the FC constant
    Cart_dyn = extra.get_dyn(ff_calculator.model)

    struct = CC.Structure.Structure(2)

    # We must setup the masses (in Ha)
    struct.masses = {"H" : 918.58996499058958}

    struct.coords = np.array([[-re/2., 0., 0.], [+re/2., 0., 0.]])

    struct.unit_cell = 3. * np.eye(3)

    struct.has_unit_cell = True

    Cart_dyn = CC.Phonons.Phonons(struct)

    Cart_dyn.dynmats[0] = k_harm * np.eye(6)
    
    w, pols = Cart_dyn.DiagonalizeSupercell()
    if np.any(w * conv.RY_TO_mEV < 1.):
        raise ValueError('The FC matrix is not positive definite!')

    initial_pos = np.copy(struct.coords)
    # In angstrom
    x_range = np.linspace(-re/5., re/5., 200) 
    delta_x = x_range[1] - x_range[0]

    directory = 'Results'
    if os.path.isdir(directory):
        shutil.rmtree(directory)

    os.mkdir(directory)
    os.chdir(directory)

    for atom in range(2):
        for coord in range(3):
            energy = []
            force  = []

            for i in range(len(x_range)):
                Cart_dyn.structure.coords = np.copy(initial_pos)
                Cart_dyn.structure.coords[atom,coord] = np.copy(initial_pos[atom,coord] + x_range[i])
                
                results = ff_calculator.calculate(atoms = Cart_dyn.structure.get_ase_atoms())

                energy.append(results['energy'] )
                force.append(results['forces'][atom,coord] )

            fig, ax = plt.subplots(2, 1, figsize = (10,10))

            force_num = -np.gradient(np.asarray(energy), delta_x)

            y = initial_pos[atom,coord] + x_range
            ax[0].plot(y, np.asarray(energy) * 1000, 'ro', label='Energy (meV)')
            ax[1].set_ylabel('Energy (meV)')
            ax[0].legend()

            ax[1].plot(y, force_num, 'ro', label = 'Numerical')
            ax[1].plot(y, force, 'bx', label = 'Exact')
            ax[1].set_ylabel('Forces (eV/Angstrom)')
            ax[1].set_xlabel('position at = {} coord = {} (Angstrom)'.format(atom, coord))
            ax[1].legend()

            fig.tight_layout()
            plt.savefig('forces_at={}_comp={}.eps'.format(atom,coord))
            
            if np.any(np.abs(force_num - force)>1e-1):
                raise ValueError('The force field is not working!')
    return 

if __name__ == "__main__":
    test_forces()

