import numpy
from functools import reduce

def mrcc_calc(calc_level):
    calc_level = calc_level.strip()
    calc_level = calc_level.upper()
    icalc = 1
    if calc_level == "CCSD":
        ex_level = 2
    elif calc_level == "CCSDT":
        ex_level = 3
    elif calc_level == "CCSDTQ":
        ex_level = 4
    elif calc_level == "CCSDT-1":
        ex_level = 3
        icalc = 5
    elif calc_level == "CCSDT-1B":
        ex_level = 3
        icalc = 6
    elif calc_level == "CCSDT-3":
        ex_level = 3
        icalc = 8
    elif calc_level == "CCSD(T)":
        ex_level = 3
        icalc = 3
    elif calc_level == "CCSDT(Q)":
        ex_level = 4
        icalc = 3
    elif calc_level == "CCSD(T)_L":
        ex_level = 3
        icalc = 4
    elif calc_level == "CCSDT(Q)_L":
        ex_level = 4
        icalc = 4
    elif calc_level == "CISD":
        ex_level = 2
        icalc = 0
    elif calc_level == "CCSDT[Q]":
        ex_level = 4
        icalc = 2
    elif calc_level[:3] == "CC(":
        if len(calc_level.split("(")) == 2 and calc_level[-1] == ")":
            # CC(n) arbitrary-order CC
            ex_level = int(calc_level[3:-1])
        elif len(calc_level.split("(")) == 2 and calc_level[-1] == "]":
            # CC(n-1)[n] arbitrary-order non-iterative CC[]
            tmp = calc_level.split("[")[1]
            ex_level = int(tmp[:-1])
            icalc = 2
        elif len(calc_level.split("(")) == 3 and calc_level[-1] == ")":
            # CC(n-1)(n) arbitrary-order non-iterative CC()
            tmp = calc_level.split("(")[2]
            ex_level = int(tmp[:-1])
            icalc = 3
        elif len(calc_level.split("(")) == 3 and calc_level[-3:] == ")_L":
            # CC(n-1)(n)_L arbitrary-order non-iterative CC()_Lambda
            tmp = calc_level.split("(")[2]
            ex_level = int(tmp[:-3])
            icalc = 4
        else:
            raise NotImplementedError("Input arbitrary order calculations are not supported: %s" % calc_level)
    elif calc_level[:3] == "CI(":
        if len(calc_level.split("(")) == 2:
            # CI(n) arbitrary-order CI
            ex_level = int(calc_level[3:-1])
            icalc = 0
        else:
            raise NotImplementedError("Only CC(n), CC(n-1)(n), and CC(n-1)(n)_L are supported")
    # CCn arbitrary-order CCn
    elif calc_level[:2] == "CC" and calc_level[2].isdigit():
        ex_level = int(calc_level[2:])
        icalc = 7
    else:
        raise NotImplementedError("Input method not supported")
    return ex_level, icalc

def mrcc_config(method, nsing = 1, ntrip = 0, ndoub = 0, irestart = 0, iden = 0, isym = 1, closed_shell = 1, icanonical = 1, tol = 7, freq = 0.0, idboc = 0, mem = 4000, file = "fort.56"):
    """
    Write the input file (fort.56) for MRCC program
    Args:
        method: str
            The method for MRCC calculation. It will be turned into excitation level and calculation flag by mrcc_calc.
        nsing: int
            Number of singlet roots.
        ntrip: int
            Number of triplet roots.
        ndoub: int
            Number of doublet roots.
        irestart: int
            Whether to restart the calculation.
        iden: int
            Flag to choose which density to calculate.
            0: no density calculation
            1: unrelaxed density
            2: relaxed density (correct one to use for properties)
            >2: for higher-order properties
        isym: int
            Symmetry flag, determine the symmetry for target state.
        closed_shell: int
            Flag for RHF calculation. RHF -> 1, ROHF/UHF -> 0.
        icanonical: int
            Flag for canonical orbitals. RHF/UHF -> 1, ROHF -> 0.
        tol: int
            Set the CC convergence threshold in MRCC to 10^(-tol).
        freq: float
            Frequency for frequency-dependent polarizabilities.
        idboc: int
            Flag for diagonal Born-Oppenheimer correction.
        mem: int
            Memory in MB.
    """
    iex, icalc = mrcc_calc(method) # excitation level and calculation flag
    nacto = nactv = idiag = inconver = imaxex = isacc = 0 # not used
    freq = 0.00000 # for calculation of freq-dependent polarizabilities
    ispatial = closed_shell # whether to use spatial orbital
    # if use ROHF and semicanonical orbitals, set ispatial to 1
    
    with open(file, "w") as ofs:
        ofs.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%d\t%d\n" % (iex, nsing, ntrip, irestart, icalc, iden, inconver, isym, idiag, closed_shell, ispatial, icanonical, ndoub, nacto, nactv, tol, imaxex, isacc, freq, idboc, mem))
        # Below are comments. MRCC won't read them.
        ofs.write("ex.lev, nsing, ntrip, restart, calc, dens, conver, symm, diag, closed-shell, spatial, canonical, ndoub, nacto, nactv, tol,maxex, sacc, freq, dboc, mem\n")


def mrcc_moints(h1e, eri, nelec, nmo, nuc_attraction, rhf = True, symm = False, file = "fort.55"):
    if symm:
        raise NotImplementedError("Molecular symmetry is not supported")
    else:
        symm_list = numpy.ones(nmo, dtype = int)

    output = "\t%d\t%d\n" % (nmo, nelec)
    for i in range(nmo):
        output += "\t%d" % (symm_list[i])
    output += "\n 150000\n" # seems not used
    if not rhf:
        raise NotImplementedError("Only RHF is supported")
    else:
        assert h1e.shape == (nmo, nmo) and eri.shape == (nmo*(nmo+1)//2, nmo*(nmo+1)//2)
        assert nelec % 2 == 0
    from pyscf.tools import fcidump
    with open(file, "w") as fout:
        fout.write(output)
        fcidump.write_eri(fout, eri, nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
        fcidump.write_hcore(fout, h1e, nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
        output_format = fcidump.DEFAULT_FLOAT_FORMAT + '  0  0  0  0\n'
        fout.write(output_format % nuc_attraction)

    # orbital mapping
    with open(file, "a") as fout:
        for i in range(nmo):
            fout.write("\t%d" % (i+1))
        fout.write("\n")
        for i in range(nmo):
            fout.write("\t%d" % (i+1))

    



def mrcc_interface(mf, method, frozen_core = 0, mol_sym = False, tol = 7, mem = 4000):
    e_nuc = mf.energy_nuc()
    if type(mf) == scf.hf.RHF:
        rhf = True
        nelec = mf.mol.nelec[0] + mf.mol.nelec[1] - frozen_core*2
        mo_coeff = mf.mo_coeff[:,frozen_core:]
        h1e = reduce(numpy.dot, (mo_coeff.T.conj(), mf.get_hcore(), mo_coeff))
        from pyscf import ao2mo
        if mf._eri is None:
            if getattr(mf, 'exxdiv', None):  # PBC system
                raise NotImplementedError("PBC system is not supported")
            else:
                eri = ao2mo.full(mf.mol, mo_coeff)
        else:  # Handle cached integrals or customized systems
            eri = ao2mo.full(mf._eri, mo_coeff)
        
        # If there is frozen core, we need to add an energy correction from core contribution to the reference energy
        if frozen_core > 0:
            occ_tmp = numpy.zeros_like(mf.mo_occ)
            occ_tmp[:frozen_core] = 2
            den_tmp = mf.make_rdm1(mf.mo_coeff, occ_tmp)
            print(mf.energy_elec(dm = den_tmp)[0])
            e_nuc += mf.energy_elec(dm = den_tmp)[0]
    else:
        raise NotImplementedError("Only RHF is supported")
    
    mrcc_config(method, tol = tol, mem = mem)
    mrcc_moints(h1e, eri, nelec, mo_coeff.shape[1], e_nuc, rhf, mol_sym)

    return

from pyscf import gto, scf
import os

mole = gto.Mole()
mole.atom = '''
H     -0.00000000     0.00000000     1.61507743
F      0.00000000     0.00000000    -0.08567644
'''
mole.basis = "ccpvdz"
mole.unit = "Bohr"

mf = scf.RHF(mole)
mf.kernel()
mrcc_interface(mf, "CCSDT", frozen_core=1)
os.system("dmrcc")

