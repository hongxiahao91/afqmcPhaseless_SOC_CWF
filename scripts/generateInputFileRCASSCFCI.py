import os
import numpy as np
import sys; sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
from molecule import *
from rhf import *
from model import *
from rcas import *

atoms             = [['O',(0.000,0.000,0.000)]]
chrg              = 1
spn               = 1
basis             = 'sto-3g'
psp               = None
sym               = 'c1'
ncas              = 5
neleca            = 4
nelecb            = 3

mol = setupMolecule(atoms, chrg, spn, basis, psp, sym)
rhf  = doROHFCalculation(mol, 1e-12, 1e-8, 0.4, 16000, 5000)
canonic = CanonicalBais(mol, rhf, 1e-8)

writeInputForModel(mol, rhf, canonic, 1e-12)

writeROHFSD2is(mol, rhf, canonic, "phi.dat", noise=1e-5)

rcas = doRCASSCFCalculation(mol, rhf, ncas, (neleca, nelecb), True, 16000, 16000, 1e-15, 200)
rcas = doRCASCICalculation(mol, rhf, ncas, (neleca, nelecb), 3, rcas.mo_coeff, True, 16000, 16000, 1e-12, 200)
rcasInfo = RCasInfo(canonic, rcas, rcas.ci[0])
rcasInfo.write()
writeRCasMDCas2s(mol, rcasInfo, 1e-15, "phiT.dat")

# rcasInfo = RCasInfo()
# writeRCasMDCas2s(mol, rcasInfo, 1e-3, "phiT.dat")

# from fci import *
# doFCICalculation(mol, rhf.mo_coeff)