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
