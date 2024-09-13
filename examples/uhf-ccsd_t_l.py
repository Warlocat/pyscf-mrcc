from pyscf import gto, scf
from mrcc import mrcc_interface
import os

mole = gto.Mole()
mole.atom = '''
F    0.00000000     0.00000000    -0.00000000
'''
mole.basis = "ccpvdz"
mole.unit = "Bohr"
mole.spin = 1

mf = scf.UHF(mole)
mf.kernel()

# CCSD(T)_Lambda calculation
# see Chem. Phys. Lett. 1997, 281, 130.
#     Int. J. Quantum Chem. 1998, 70, 601.
#     J. Chem. Phys. 2008, 128, 44110.
mrcc_interface(mf, "CCSD(T)_L", frozen_core=1)
os.system("dmrcc")
# reference energy
# Total CCSD[T] energy [au]:                   -99.527625642991
# Total CCSD(T) energy [au]:                   -99.527574087159
# Total CCSD(T)_L energy [au]:                 -99.527600128429
