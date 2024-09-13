from pyscf import gto, scf
from mrcc import mrcc_interface
import os

mole = gto.Mole()
mole.atom = '''
H   -0.00000000     0.00000000     1.61507743
F    0.00000000     0.00000000    -0.08567644
'''
mole.basis = "ccpvdz"
mole.unit = "Bohr"

mf = scf.RHF(mole)
mf.kernel()

mrcc_interface(mf, "CCSDTQ", frozen_core=1)
os.system("dmrcc")
# reference energy: -100.22816276
