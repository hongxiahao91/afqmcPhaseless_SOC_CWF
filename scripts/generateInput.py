import os
import numpy as np
import sys; ##sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
from molecule import *
from rhf import *
from model import *
from uhf import *

global NMO, N, Nup, Ndn, spin
NMO = 14
SMO = 2*NMO
N = 5
Nup = 3
Ndn = 2
spin = 1


def readFCIDUMPSOC():
	Tout = np.zeros(shape=(SMO,SMO),dtype=complex)
	Vout = np.zeros(shape=(SMO,SMO, SMO, SMO),dtype=float)
	P_upup = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	P_dndn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	P_updn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	EngShf = []
	SumEngShf = 0.0
	count = 0
	fin = open('./FCIDUMPSOC','r')
	while(True):
		line = fin.readline()
		if not line:
			break
		tmp = line.split()
		if len(tmp) == 6:
			real = float(tmp[0])
			imag = float(tmp[1])
			ind_i = int(tmp[2])
			ind_l = int(tmp[3])
			ind_j = int(tmp[4])
			ind_k = int(tmp[5])
			value = real+ 1j* imag
			if ind_i == 0 and ind_l == 0 and ind_j == 0 and ind_k == 0:
				EngShf[count] = real+ 1j* imag
				SumEngShf += EngShf[count]
				count = count+1
				continue
			if count == 0:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i][2*ind_l][2*ind_j][2*ind_k] = real
					P_upup[ind_i][ind_l][ind_j][ind_k] = real
			elif count == 1:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1] = real
					P_dndn[ind_i][ind_l][ind_j][ind_k] = real
			elif count == 2:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i][2*ind_l][2*ind_j+1][2*ind_k+1] = real
					Vout[2*ind_i+1][2*ind_l+1][2*ind_j][2*ind_k] = real
					P_updn[ind_i][ind_l][ind_j][ind_k] = real
			elif count ==3:
				if abs(value) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					Tout[ind_i][ind_l] = value
			elif count ==4:
				if abs(value) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					Tout[ind_i+NMO][ind_l+NMO] = value
			elif count ==5:
				if abs(value) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					Tout[ind_i][ind_l+NMO] = value
					Tout[ind_i+NMO][ind_l] = value

	return Tout, Vout

def readFCIDUMP():
	Tout = np.zeros(shape=(SMO,SMO),dtype=complex)
	Vout = np.zeros(shape=(SMO,SMO, SMO, SMO),dtype=float)
	P_upup = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	P_dndn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	P_updn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	EngShf = []
	SumEngShf = 0.0
	count = 0
	fin = open('./FCIDUMPSOC','r')
	while(True):
		line = fin.readline()
		if not line:
			break
		tmp = line.split()
		if len(tmp) == 5:
			real = float(tmp[0])
			ind_i = int(tmp[1])
			ind_l = int(tmp[2])
			ind_j = int(tmp[3])
			ind_k = int(tmp[4])
			if ind_i == 0 and ind_l == 0 and ind_j == 0 and ind_k == 0:
				EngShf[count] = real
				SumEngShf += EngShf[count]
				count = count+1
				continue
			if count == 0:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i][2*ind_l][2*ind_j][2*ind_k] = real
					P_upup[ind_i][ind_l][ind_j][ind_k] = real
			elif count == 1:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1] = real
					P_dndn[ind_i][ind_l][ind_j][ind_k] = real
			elif count == 2:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					ind_j = ind_j -1
					ind_k = ind_k -1
					Vout[2*ind_i][2*ind_l][2*ind_j+1][2*ind_k+1] = real
					Vout[2*ind_i+1][2*ind_l+1][2*ind_j][2*ind_k] = real
					P_updn[ind_i][ind_l][ind_j][ind_k] = real
			elif count ==3:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					Tout[ind_i][ind_l] = real
			elif count ==4:
				if abs(real) > 0.00000001:
					ind_i = ind_i -1
					ind_l = ind_l -1
					Tout[ind_i+NMO][ind_l+NMO] = real


	return Tout, Vout, P_upup, P_dndn, P_updn

if __name__ == '__main__':

	t,V,pUpUp, pDnDn, pUnDn = readFCIDUMP()
	V2 = np.transpose(V, (0,1,3,2))
	diff = V2-V
	print(diff)
	V = V.reshape( SMO*SMO, SMO*SMO )
	w, vec = np.linalg.eigh(V)

	eigenNum = NMO**2

	f = open("model_param", 'w')
	f.write( '{:16d} \n'.format(NMO) )
	f.write( '{:16d} \n'.format(N) )
	f.write( '{:16d} \n'.format(Nup) )
	f.write( '{:16d} \n'.format(Ndn) )
	f.write( '{:16d} \n'.format(eigenNum) )
	for i in range( 2*NMO ):
		for j in range( 2*NMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(t[j, i].real, t[j, i].imag))
	for i in range( 2*NMO ):
		for j in range( 2*NMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(K[j, i].real, K[j, i].imag))
	for i in range( 2*NMO ):
		for j in range( 2*NMO ):
			for l in range( eigenNum ):
			f.write( '{:26.18e} {:26.18e} \n'.format(V[l, j, i].real, V[l, j, i].imag))
	for i in range( eigenNum ):
		f.write( '{:26.18e} \n'.format(0.0) )
	f.close()

#	f = h5py.File("model_param", "w")
#	f.create_dataset("L",              (1,),                        data=[NMO],                 dtype='int')
#	f.create_dataset("N",              (1,),                        data=[N],                   dtype='int')
#	f.create_dataset("Nup",            (1,),                        data=[Nup],                 dtype='int')
#	f.create_dataset("Ndn",            (1,),                        data=[Ndn],                 dtype='int')
#	f.create_dataset("choleskyNumber", (1,),                        data=[choleskyNum],         dtype='int')
#	f.create_dataset("t",              (4*NMO**2,),                 data=t.ravel(),             dtype='complex')
#	f.create_dataset("K",              (4*NMO**2,),                 data=K.ravel(),             dtype='complex')
#	f.create_dataset("choleskyVecs",   (choleskyNum*NMO**2,),       data=V.ravel(),             dtype='float64')
#	f.create_dataset("choleskyBg",     (choleskyNum,),              data=np.zeros(choleskyNum), dtype='float64')
#	f.close()


#	V = readVMatrix()

