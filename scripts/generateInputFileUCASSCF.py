import os
import numpy as np
import sys; sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
from molecule import *
from rhf import *
from model import *
from uhf import *
from ucas import *

atoms             = [['O',(0.000,0.000,0.000)]]
chrg              = 2
spn               = 2
basis             = 'sto-3g'
psp               = None
sym               = 'c1'
ncas              = 5
neleca            = 4
nelecb            = 2

mol = setupMolecule(atoms, chrg, spn, basis, psp, sym)
rhf  = doROHFCalculation(mol, 1e-12, 1e-8, 0.4, 16000, 5000)
canonic = CanonicalBais(mol, rhf, 1e-8)

writeInputForModel(mol, rhf, canonic, 1e-12)

uhf = doUHFCalculation(mol, 1e-12, 1e-8, 0.4, 16000, 5000)
writeUHFSD2s(mol, uhf, canonic, "phi.dat", 1e-5)

ucas = doUCASSCFCalculation(mol, uhf, ncas, (neleca, nelecb), 16000, 16000, 1e-12, 200)
ucasInfo = UCasInfo(canonic, ucas)
ucasInfo.write()
writeUCasMDCas2s(mol, ucasInfo, 1e-8, "phiT.dat")

# ucasInfo = UCasInfo()
# writeUCasMDCas2s(mol, ucasInfo, 1e-3, "phiT.dat")

from fci import *
doFCICalculation(mol, rhf.mo_coeff)
