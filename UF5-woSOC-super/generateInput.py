import os
import numpy as np
import sys; ##sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
import h5py

global NMO, SMO, N, Nup, Ndn, spin
NMO = 14
SMO = 2*NMO
N = 5
Nup = 3
Ndn = 2
spin = 1

def readTUpMatrix(fout):
	fin=open('./TUpMatrix.txt','r')
	while(True):
		line= fin.readline()
		if(not line) :
			break
		tmp = line.split()
		real = float(tmp[0])
		imag = float(tmp[1])
		ind_i = int(tmp[2])
		ind_j = int(tmp[3])
		fout[ind_i-1][ind_j-1] = real+ 1j* imag
	fin.close()
	return fout

def readTDnMatrix(fout):
	fin=open('./TDnMatrix.txt','r')
	while(True):
		line= fin.readline()
		if(not line) :
			break
		tmp = line.split()
		real = float(tmp[0])
		imag = float(tmp[1])
		ind_i = int(tmp[2])
		ind_j = int(tmp[3])
		fout[ind_i+NMO-1][ind_j+NMO-1] = real+ 1j* imag
	fin.close()
	return fout

def readTSOCMatrix(fout):
	fin=open('./TSOCMatrix.txt','r')
	while(True):
		line= fin.readline()
		if(not line) :
			break
		tmp = line.split()
		real = float(tmp[0])
		imag = float(tmp[1])
		ind_i = int(tmp[2])
		ind_j = int(tmp[3])
		fout[ind_i-1][ind_j+NMO-1] = real+ 1j* imag
		fout[ind_i+NMO-1][ind_j-1] = real+ 1j* imag
	fin.close()
	return fout

def readVMatrix(fout):
	fin=open('./VMatrix.txt','r')
	while(True):
		line= fin.readline()
		if(not line) :
			break
		tmp = line.split()
		value = float(tmp[0])
		ind_i = int(tmp[1])-1
		ind_l = int(tmp[2])-1
		ind_j = int(tmp[3])-1
		ind_k = int(tmp[4])-1

		if abs(value) > 0.00000001:
			fout[2*ind_i][2*ind_l][2*ind_j][2*ind_k] = value
			fout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1] = value
			fout[2*ind_i][2*ind_l][2*ind_j+1][2*ind_k+1] = value
			fout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1] = value
	fin.close()
	return fout


if __name__ == '__main__':

	Tout = np.zeros(shape=(2*NMO,2*NMO),dtype=complex)
	tup = readTUpMatrix(Tout)
	t = readTDnMatrix(tup)
	#t = readTSOCMatrix(tdn)
	#print(t)

	K = t.copy()

	Vout = np.zeros(shape=(SMO,SMO, SMO, SMO),dtype=float)
	V = readVMatrix(Vout)
	V = V.reshape( SMO*SMO, SMO*SMO )  
	#Vdiag = V.diagonal().copy()
	w, vec = np.linalg.eigh(V)
	#vec2 = vec.sum(axis=1)
	#vec2 = vec2.reshape(SMO, SMO)
	#print (vec2 - vec2.T)
	final = np.dot(np.dot(vec, np.diag(w)),vec.T)
	print(V-final)		

	choleskyNum = NMO**2

#	V.resize(NMO,NMO,choleskyNum)
	#print(V)

	f = open("model_param", 'w')
	f.write( '{:16d} \n'.format(NMO) )
	f.write( '{:16d} \n'.format(N) )
	f.write( '{:16d} \n'.format(Nup) )
	f.write( '{:16d} \n'.format(Ndn) )
	f.write( '{:16d} \n'.format(choleskyNum) )
	for i in range( 2*NMO ):
		for j in range( 2*NMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(t[j, i].real, t[j, i].imag))
	for i in range( 2*NMO ):
		for j in range( 2*NMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(K[j, i].real, K[j, i].imag))
#	for i in range( choleskyNum ):
#		for j in range( NMO ):
#			for l in range( NMO ):
#				f.write( '{:26.18e} {:26.18e} \n'.format(V[l, j, i].real, V[l, j, i].imag))
	for i in range( choleskyNum ):
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

	w, v = np.linalg.eigh(K)
	print ("Eigen of Kmatrix")
	print(w)

	print ("FE total energy")
	print( "{:<26.18f}".format( np.sum( w[0:N] ) ) )
