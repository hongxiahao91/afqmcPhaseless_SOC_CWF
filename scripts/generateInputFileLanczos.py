import os
import subprocess
import numpy as np
import sys; sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
from molecule import *
from rhf import *
from model import *
from rcas import *
from lanczos import *

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

writeROHFSD2is(mol, rhf, canonic, "phi.dat", noise=1e-5)

rcas = doRCASSCFCalculation(mol, rhf, ncas, (neleca, nelecb), False, 16000, 16000, 1e-12, 200)

subprocess.call('rm -rf lanczos', shell=True)
subprocess.call('mkdir lanczos', shell=True)

orbital = np.dot( canonic.XInv, rcas.mo_coeff )
np.savetxt("lanczos/rcasOrbital.dat", orbital.flatten(), fmt='%26.18e')
writeInputForLanczos(orbital, "lanczos/model_param")

matrixSize       = 20           #Lanczos matrix size
accuracy         = 1E-10        #when new b smaller than accuracy, converge
convergeFlag     = 'E'          #'E' or 'W', converge by wave function or energy
maxLoop          = 100          #The max Lanczos matrix loop
litForProjection = 0.01         #When b is smaller, we need to project wave function.
lanwfsFlag       = 'F'          #'R' or 'F', 'R' use recurse wf, 'F' store full Lanczos wf

f = open("lanczos/lanczos_param", 'w')
f.write( '{:16d} \n'.format(matrixSize) )
f.write( '{:26.18e} \n'.format(accuracy)   )
f.write( '{:>16} \n'.format(convergeFlag))
f.write( '{:16d} \n'.format(maxLoop) )
f.write( '{:26.18e} \n'.format(litForProjection)   )
f.write( '{:>16} \n'.format(lanwfsFlag)   )
f.close()

os.chdir('lanczos')

subprocess.call('/Users/hshi/myproject/AFQMCLAB/tutorials/lanczosMolecueForMDCas2s/build/lanczosHubbard 3', shell=True)
writeRCasFromLanczos(0)

rcasInfo = RCasInfo()
writeRCasMDCas2s(mol, rcasInfo, 1e-15, "../phiT.dat")

os.chdir('..')

# from fci import *
# doFCICalculation(mol, rhf.mo_coeff)