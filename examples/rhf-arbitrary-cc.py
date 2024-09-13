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

# MRCC supports many arbitrary-order CC methods.

# CCn: CC2, CC3, ...
mrcc_interface(mf, "CC2", frozen_core=1)
os.system("dmrcc")
# reference energy: -100.22128627

# CC(n): CCSD, CCSDT, ...
mrcc_interface(mf, "CC(3)", frozen_core=1)
os.system("dmrcc")
# reference energy: -100.22778038

# CC(n-1)(n): CCSD(T), CCSDT(Q), ...
mrcc_interface(mf, "CC(3)(4)", frozen_core=1)
os.system("dmrcc")
# reference energy: -100.22818183