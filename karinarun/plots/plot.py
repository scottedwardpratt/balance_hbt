import matplotlib.pyplot as plt
import numpy as np

filename1 = "cf_exp.dat"
filename2 = "cf_outsidelong.dat"

mydata1 = np.genfromtxt(filename1, delimiter = None)[:,:]
mydata2 = np.genfromtxt(filename2, delimiter = None)[:,:]

trimdata = mydata2[:30,:]
print(trimdata)
#mypion1 = np.genfromtxt(filepion1, delimiter = None)[:,:]
#mypion2 = np.genfromtxt(filepion2, delimiter = None)[:,:]

#mykaon1 = np.genfromtxt(filekaon1, delimiter = None)[:,:]
#mykaon2 = np.genfromtxt(filekaon2, delimiter = None)[:,:]


plt.plot(mydata1[:,0],mydata1[:,1], 'o')
plt.plot(trimdata[:,0]/1000*2,1+trimdata[:,1], 'x')

#plt.plot(mypion1[:,0],mypion1[:,1], 'o')
#plt.plot(mypion2[:,0]*1000,mypion2[:,3], 'x')

#plt.plot(mykaon1[:,0],mykaon1[:,1], 'o')
#plt.plot(mykaon2[:,0]*1000,mykaon2[:,3], 'x')

plt.title('CF')
plt.xlabel('q')
plt.ylabel('CFout')
plt.show()


