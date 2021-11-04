import os
import numpy as np
import sys;  ##sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/pyscf" )
#from molecule import *
#from rhf import *
#from model import *
#from uhf import *

global NMO, SMO, N, Nup, Ndn, spin
NMO = 10
SMO = 2*NMO
N = 5
Nup = 3
Ndn = 2
spin = 1


def readFCIDUMP():
	Tout = np.zeros(shape=(SMO,SMO),dtype=complex)
	One = np.zeros(shape=(SMO,SMO),dtype=float)
	Vupup = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	Vdndn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	Vupdn = np.zeros(shape=(NMO,NMO, NMO, NMO),dtype=float)
	Vout = np.zeros(shape=(SMO,SMO, SMO, SMO),dtype=float)
	EngShf = []
	SumEngShf = 0.0
	count = 0
	fin = open('./FCIDUMP','r')
	print('open succ')
	line = fin.readline()
	line = fin.readline()
	line = fin.readline()
	line = fin.readline()
	while(True):
		line = fin.readline()
		#print(line)
		if (not line):
			break
		tmp = line.split()
		real = float(tmp[0])
		imag = float(tmp[1])
		ind_i = int(tmp[2])
		ind_l = int(tmp[3])
		ind_j = int(tmp[4])
		ind_k = int(tmp[5])
		value = real+ 1j* imag
		if (ind_i == 0 and ind_l == 0 and ind_j == 0 and ind_k == 0):
			#print(value)
			EngShf.append(value)
			SumEngShf += EngShf[count]
			count = count+1
			continue
		if count == 0:
			if abs(real) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				ind_j = ind_j -1
				ind_k = ind_k -1
				Vout[2*ind_i][2*ind_l][2*ind_j][2*ind_k] = real
				Vupup[ind_i][ind_l][ind_j][ind_k]= real
				if (ind_l == ind_j):
					One[ind_i][ind_k] += real
			#print (ind_i, ind_l, ind_j, ind_k, Vout[2*ind_i][2*ind_l][2*ind_j][2*ind_k])
		elif count == 1:
			if abs(real) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				ind_j = ind_j -1
				ind_k = ind_k -1
				Vout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1] = real
				Vdndn[ind_i][ind_l][ind_j][ind_k]= real
				if (ind_l == ind_j):
					One[ind_i+NMO][ind_k+NMO] += real
			#print (Vout[2*ind_i+1][2*ind_l+1][2*ind_j+1][2*ind_k+1])
		elif count == 2:
			if abs(real) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				ind_j = ind_j -1
				ind_k = ind_k -1
				Vout[2*ind_i][2*ind_l][2*ind_j+1][2*ind_k+1] = real
				Vout[2*ind_i+1][2*ind_l+1][2*ind_j][2*ind_k] = real
				Vupdn[ind_i][ind_l][ind_j][ind_k]= real
				#print (Vout[2*ind_i][2*ind_l][2*ind_j+1][2*ind_k+1])
		elif count ==3:
			if abs(value) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				Tout[ind_i][ind_l] = value
		elif count ==4:
			if abs(value) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				Tout[ind_i+NMO][ind_l+NMO] = value
		elif count ==5:
			if abs(value) > 1.0e-8:
				ind_i = ind_i -1
				ind_l = ind_l -1
				Tout[ind_i][ind_l+NMO] = value
				Tout[ind_i+NMO][ind_l] = value
		#print(Vout)
	#print ("Most of the elements are 0, so you didn't see the nonzero ones!")
	# for ii in range(2*NMO):
	# 	for jj in range(2*NMO):
	# 		for kk in range(2*NMO):
	# 			for ll in range(2*NMO):
	# 				if Vout[ii][jj][kk][ll] != 0:
	# 					print (ii,jj,kk,ll,Vout[ii][jj][kk][ll])
	# for ii in range(2*NMO):
	# 	for jj in range(2*NMO):
	# 		if Tout[ii][jj] !=0:
	# 			print (ii, jj, Tout[ii][jj])
	# print(Tout)
	fin.close()
	#print(Vout)
	return Tout, One, Vout, Vupup, Vdndn, Vupdn, SumEngShf



if __name__ == '__main__':

	t, One, V, Pupup, Pdndn, Pupdn, EngShf_total = readFCIDUMP()
	K = t.copy()
	K = K - 0.5*One
	#print(t)
	#exit(1)
	#print(V)
	V = np.transpose(V, (0,1,3,2))
	#diff = V2-V
	# for ii in range(2*NMO):
	# 	for jj in range(2*NMO):
	# 		for kk in range(2*NMO):
	# 			for ll in range(2*NMO):
	# 				if abs(V2[ii][jj][kk][ll]-V[ii][jj][kk][ll]) > 0.000001:
	# 					print (ii,jj,kk,ll,V[ii][jj][kk][ll],V2[ii][jj][kk][ll]-V[ii][jj][kk][ll])
	#print(diff)

	V = V.reshape( SMO*SMO, SMO*SMO )
	w, vec = np.linalg.eigh(V)

	# print(vec.shape)
	# print(w.shape)
	# print(vec[:,0])
	eigen = []
	eigenVec = np.zeros(shape=(SMO*SMO, 1),dtype=float)
	# print(eigenVec.shape)
	for ii in range(SMO*SMO):
		if abs(w[ii]) >= 1.0e-8:
			eigen.append(w[ii])
			eigenVec = np.concatenate((eigenVec,vec[:,ii:ii+1]),axis=1)
	eigenNum = len(eigen)
	#print(eigenVec.shape)
	eigenVec = np.delete(eigenVec, 0, axis=1)
	#print(eigenVec.shape)
	# print (eigenVec[:,0])
	#eigenVec = eigenVec.reshape(SMO,SMO,len(eigen))
	# print(eigenVec.shape)
	phoUp = np.zeros(shape=(1,eigenNum),dtype=float)
	phoDn = np.zeros(shape=(1,eigenNum),dtype=float)
	for ii in range(SMO):
		for jj in range(SMO):
			if(ii%2==0 and jj%2 ==0):
				phoUp = np.concatenate((phoUp, eigenVec[ii*SMO+jj:ii*SMO+jj+1,:]),axis=0)
			if(ii%2==1 and jj%2 ==1):
				phoDn = np.concatenate((phoDn, eigenVec[ii*SMO+jj:ii*SMO+jj+1,:]),axis=0)
	phoUp = np.delete(phoUp, 0, axis=0)
	phoDn = np.delete(phoDn, 0, axis=0)
	phoUp = np.transpose(phoUp)
	phoDn = np.transpose(phoDn)
	phoUp = phoUp.reshape(eigenNum, NMO, NMO)
	phoDn = phoDn.reshape(eigenNum, NMO, NMO)
	#print(phoUp.shape)
	#print(phoDn.shape)
	#print(phoUp[:,:,0])

	phoUpSum = phoUp + np.transpose(phoUp, (0,2,1))
	phoUpMinus = phoUp - np.transpose(phoUp, (0,2,1))
	phoDnSum = phoDn + np.transpose(phoDn, (0,2,1))
	phoDnMinus = phoDn = np.transpose(phoDn, (0,2,1))

	#print(phoUpMinus[:,:,0])
	phoSum = np.zeros(( eigenNum, 2 * NMO, 2 * NMO), dtype='float', order='F')
	phoMinus = np.zeros(( eigenNum, 2 * NMO, 2 * NMO), dtype='float', order='F')
	for i in range(eigenNum):
		phoSum[i, 0:NMO, 0:NMO] = phoUpSum[i]
		phoSum[i,NMO:2*NMO, NMO:2*NMO] = phoDnSum[i]
		phoMinus[i, 0:NMO, 0:NMO] = phoUpMinus[i]
		phoMinus[i,NMO:2*NMO, NMO:2*NMO] = phoDnMinus[i]

	#exit(1)
	f = open("model_param", 'w')
	f.write( '{:16d} \n'.format(NMO) )
	f.write( '{:16d} \n'.format(N) )
	f.write( '{:16d} \n'.format(Nup) )
	f.write( '{:16d} \n'.format(Ndn) )
	f.write( '{:16d} \n'.format(eigenNum) )
	f.write( '{:26.18e} \n'.format(EngShf_total.real) )
	for i in range( SMO ):
		for j in range( SMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(t[j, i].real, t[j, i].imag))
	for i in range( SMO ):
		for j in range( SMO ):
			f.write( '{:26.18e} {:26.18e} \n'.format(K[j, i].real, K[j, i].imag))
	for i in range(eigenNum):
		f.write( '{:26.18e} \n'.format(eigen[i]))
	for l in range( eigenNum ):
		for i in range( SMO ):
			for j in range( SMO ):
				f.write( '{:26.18e} \n'.format(phoSum[l, j, i].real))
	for l in range( eigenNum ):
		for i in range( 2*NMO ):
			for j in range( 2*NMO ):
				f.write( '{:26.18e} \n'.format(phoMinus[l, j, i].real))

	for i in range(NMO):
		for j in range(NMO):
			for k in range(NMO):
				for l in range(NMO):
					f.write( '{:26.18e} \n'.format(Pupup[j, i, l, k].real))

	for i in range(NMO):
		for j in range(NMO):
			for k in range(NMO):
				for l in range(NMO):
					f.write( '{:26.18e} \n'.format(Pdndn[j, i, l, k].real))

	for i in range(NMO):
		for j in range(NMO):
			for k in range(NMO):
				for l in range(NMO):
					f.write( '{:26.18e} \n'.format(Pupdn[j, i, l, k].real))

	# for l in range( eigenNum ):
	# 	for i in range( NMO ):
	# 		for j in range( NMO ):
	# 			f.write( '{:26.18e} {:26.18e} \n'.format(phoUpSum[i, j, l].real, 0.0))
	# for l in range( eigenNum ):
	# 	for i in range( NMO ):
	# 		for j in range( NMO ):
	# 			f.write( '{:26.18e} {:26.18e} \n'.format(phoUpMinus[i, j, l].real, 0.0))
	# for l in range( eigenNum ):
	# 	for i in range( NMO ):
	# 		for j in range( NMO ):
	# 			f.write( '{:26.18e} {:26.18e} \n'.format(phoDnSum[i, j, l].real, 0.0))
	# for l in range( eigenNum ):
	# 	for i in range( NMO ):
	# 		for j in range( NMO ):
	# 			f.write( '{:26.18e} {:26.18e} \n'.format(phoDnMinus[i, j, l].real, 0.0))
	for i in range( eigenNum ):
		f.write( '{:26.18e} \n'.format(0.0) )
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

