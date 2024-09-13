'''
    This script is used to interface PySCF (https://pyscf.org/) with MRCC (https://mrcc.hu/) 
    program. The interface is to utilize the MRCC program performing high-order and various
    approximated coupled-cluster calculations.

    The interface is compatible with most versions of MRCC. It processes the PySCF integrals
    and rewrite them in MRCC format (fort.55, similar to FCIDUMP) as well as write the input 
    file (fort.56) for MRCC calculation. One can run os.system("dmrcc") to start the MRCC 
    calculation.
    
    The MRCC 2023 will raise an error. This can be solved by commenting out the 664-667 lines 
    in the MRCC source code driver.f. It might influence some of the F12 calculation in MRCC. 

    Written by Chaoqun Zhang <cq_zhang@outlook.com>.
'''
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
    '''
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
    '''
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
        # assume UHF; ROHF is not supported yet
        assert h1e.shape == (2, nmo, nmo) and eri.shape == (3, nmo*(nmo+1)//2, nmo*(nmo+1)//2)
        from pyscf.tools import fcidump
        output_format = fcidump.DEFAULT_FLOAT_FORMAT + ' %4d %4d %4d %4d\n'
        with open(file, "w") as fout:
            fout.write(output)
            fcidump.write_eri(fout, eri[0], nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (0.0, 0, 0, 0, 0))
            fcidump.write_eri(fout, eri[1], nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (0.0, 0, 0, 0, 0))
            fcidump.write_eri(fout, eri[2], nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (0.0, 0, 0, 0, 0))
            fcidump.write_hcore(fout, h1e[0], nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (0.0, 0, 0, 0, 0))
            fcidump.write_hcore(fout, h1e[1], nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (0.0, 0, 0, 0, 0))
            fout.write(output_format % (nuc_attraction, 0, 0, 0, 0))
    else:
        assert h1e.shape == (nmo, nmo) and eri.shape == (nmo*(nmo+1)//2, nmo*(nmo+1)//2)
        assert nelec % 2 == 0
        from pyscf.tools import fcidump
        output_format = fcidump.DEFAULT_FLOAT_FORMAT + ' %4d %4d %4d %4d\n'
        with open(file, "w") as fout:
            fout.write(output)
            fcidump.write_eri(fout, eri, nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fcidump.write_hcore(fout, h1e, nmo, tol=fcidump.TOL, float_format=fcidump.DEFAULT_FLOAT_FORMAT)
            fout.write(output_format % (nuc_attraction, 0, 0, 0, 0))

    # orbital mapping
    with open(file, "a") as fout:
        for i in range(nmo):
            fout.write("\t%d" % (i+1))
        fout.write("\n")
        for i in range(nmo):
            fout.write("\t%d" % (i+1))

def mrcc_interface(mf, method, frozen_core = 0, mol_sym = False, tol = 7, mem = 4000):
    '''
    The current version only supports RHF and UHF for ground state calculations.
    In principle, the extension to excited states is possible by changing nsing, ntrip, and ndoub.
    Args:
        mf: pyscf.scf.hf.RHF or pyscf.scf.uhf.UHF
            The PySCF SCF object.
        method: str
            The method for MRCC calculation. It will be turned into excitation level and calculation flag by mrcc_calc.
        frozen_core: int
            Number of frozen core orbitals. Frozen virtual orbitals are not supported yet.
        mol_sym: bool
            Whether to use molecular symmetry. Not supported yet.
        tol: int
            Set the CC convergence threshold in MRCC to 10^(-tol).
        mem: int
            Memory in MB.
    '''
    e_nuc = mf.energy_nuc()
    nelec = mf.mol.nelec[0] + mf.mol.nelec[1] - frozen_core*2
    nmo = mf.mol.nao_nr() - frozen_core
    spin = mf.mol.spin
    nsing = ntrip = ndoub = 0
    if spin == 0:
        nsing = 1
    elif spin == 1:
        ndoub = 1
    elif spin == 2:
        ntrip = 1
    else:
        raise NotImplementedError("Only singlet, doublet, and triplet are supported")

    from pyscf import scf
    if type(mf) == scf.hf.RHF:
        rhf = True
        mo_coeff = mf.mo_coeff[:,frozen_core:]
        h1e_ao = mf.get_hcore()

        # If there are frozen core orbitals, the e_nuc and h1e need to be corrected
        # to ensure the correct energy and fock matrices in the active space.
        if frozen_core > 0:
            occ_core = numpy.zeros_like(mf.mo_occ)
            occ_core[:frozen_core] = 2
            den_core = mf.make_rdm1(mf.mo_coeff, occ_core)
            e_nuc += mf.energy_elec(dm = den_core)[0]
            h1e_ao = mf.get_fock(dm = den_core)

        h1e = reduce(numpy.dot, (mo_coeff.T, h1e_ao, mo_coeff))
        from pyscf import ao2mo
        if mf._eri is None:
            if getattr(mf, 'exxdiv', None):  # PBC system
                raise NotImplementedError("PBC system is not supported")
            else:
                eri = ao2mo.full(mf.mol, mo_coeff)
        else:  # Handle cached integrals or customized systems
            eri = ao2mo.full(mf._eri, mo_coeff)
    elif type(mf) == scf.uhf.UHF:
        rhf = False
        mo_coeff_a = mf.mo_coeff[0][:,frozen_core:]
        mo_coeff_b = mf.mo_coeff[1][:,frozen_core:]
        h1e_ao_a = mf.get_hcore()
        h1e_ao_b = mf.get_hcore()

        # If there are frozen core orbitals, the e_nuc and h1e need to be corrected
        # to ensure the correct energy and fock matrices in the active space.
        if frozen_core > 0:
            occ_core = numpy.zeros_like(mf.mo_occ)
            occ_core[0][:frozen_core] = 1
            occ_core[1][:frozen_core] = 1
            den_core = mf.make_rdm1(mf.mo_coeff, occ_core)
            e_nuc += mf.energy_elec(dm = den_core)[0]
            h1e_ao = mf.get_fock(dm = den_core)
            h1e_ao_a = h1e_ao[0]
            h1e_ao_b = h1e_ao[1]

        h1e_a = reduce(numpy.dot, (mo_coeff_a.T, h1e_ao_a, mo_coeff_a))
        h1e_b = reduce(numpy.dot, (mo_coeff_b.T, h1e_ao_b, mo_coeff_b))
        from pyscf import ao2mo
        if mf._eri is None:
            if getattr(mf, 'exxdiv', None):  # PBC system
                raise NotImplementedError("PBC system is not supported")
            else:
                eri_aa = ao2mo.full(mf.mol, mo_coeff_a)
                eri_bb = ao2mo.full(mf.mol, mo_coeff_b)
                eri_ab = ao2mo.general(mf.mol, (mo_coeff_a,mo_coeff_a,mo_coeff_b,mo_coeff_b))
        else:  # Handle cached integrals or customized systems
            eri_aa = ao2mo.full(mf._eri, mo_coeff_a)
            eri_bb = ao2mo.full(mf._eri, mo_coeff_b)
            eri_ab = ao2mo.general(mf._eri, (mo_coeff_a,mo_coeff_a,mo_coeff_b,mo_coeff_b))
        h1e = numpy.array([h1e_a, h1e_b])
        eri = numpy.array([eri_aa, eri_bb, eri_ab])
    else:
        raise NotImplementedError("Only RHF and UHF are supported")
    
    mrcc_config(method, nsing, ntrip, ndoub, closed_shell = rhf, tol = tol, mem = mem)
    mrcc_moints(h1e, eri, nelec, nmo, e_nuc, rhf, mol_sym)

    return

