import matplotlib.pyplot as plt
import numpy as np

filename1 = "spectra_p.txt"
filename2 = "phenix_proton.txt"

filepion1 = "spectra_pi.txt"
filepion2 = "phenix_pion.txt"

filekaon1 = "spectra_K.txt"
filekaon2 = "phenix_kaon.txt"


mydata1 = np.genfromtxt(filename1, delimiter = None)[:,:]
mydata2 = np.genfromtxt(filename2, delimiter = None)[:,:]

#mypion1 = np.genfromtxt(filepion1, delimiter = None)[:,:]
#mypion2 = np.genfromtxt(filepion2, delimiter = None)[:,:]

#mykaon1 = np.genfromtxt(filekaon1, delimiter = None)[:,:]
#mykaon2 = np.genfromtxt(filekaon2, delimiter = None)[:,:]


plt.plot(mydata1[:,0],mydata1[:,1], 'o')
plt.plot(mydata2[:,0]*1000,mydata2[:,3], 'x')

#plt.plot(mypion1[:,0],mypion1[:,1], 'o')
#plt.plot(mypion2[:,0]*1000,mypion2[:,3], 'x')

#plt.plot(mykaon1[:,0],mykaon1[:,1], 'o')
#plt.plot(mykaon2[:,0]*1000,mykaon2[:,3], 'x')

plt.title('proton')
plt.xlabel('pt')
plt.ylabel('spectra')
plt.show()


