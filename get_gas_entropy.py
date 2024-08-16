#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 18:38:35 2022
https://janaf.nist.gov/tables/H-083.html #Database
Symmetry Number: Table 10.1 and Appendix B of C. Cramer “Essentials of Computational Chemistry”, 2nd Ed.
spin:  0 for molecules in which all electrons are paired, 0.5 for a free radical with a single unpaired electron, 1.0 for a triplet with two unpaired electrons, such as O_2.
@author: qli
"""
import numpy as np
from ase.io import read, write
from ase.thermochemistry import IdealGasThermo
from scipy import constants as con

atoms = read('./freq/POSCAR')
print(atoms)
sym = 3 # symmetry number of NH3 
spin = 0 # spin of NH3. 
tem = 300 # Temperature
out = read('./OUTCAR',  format='vasp-out')
potentialenergy = out.get_potential_energy()

model = read('./freq/POSCAR')
# model_positions = model.get_positions()

vib_energies = []
with open('./freq/OUTCAR') as f_in:
    lines = f_in.readlines()
    for num, line in enumerate(lines):
        if 'cm-1' in line:
            vib_e = float(line.rstrip().split()[-2])
            vib_energies.append(vib_e)

vib_energies = np.array(vib_energies[:-5])/1000 # For Gas, the last six are translation and rotation
# zpe = sum(vib_energies)/2

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=atoms,
                        geometry='linear',
                        symmetrynumber=3, spin=1.0)

zpe = thermo.get_ZPE_correction()
entropy = thermo.get_entropy(temperature=tem, pressure=101325,verbose=True)  # Unit eV/K
TS = tem * entropy   
G = potentialenergy + zpe - TS

print('G', G)
print('S', con.Avogadro * con.electron_volt * entropy, 'J/K/mol')  
