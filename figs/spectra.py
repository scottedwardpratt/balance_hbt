import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *

plt.figure(figsize=(6,12))
fig = plt.figure(1)

fileproton1 = "../scottrun/spectra/spectra_p.txt"
fileproton2 = "../run/phenix_data/phenix_proton.txt"

filepion1 = "../scottrun/spectra/spectra_pi.txt"
filepion2 = "../run/phenix_data/phenix_pion.txt"

filekaon1 = "../scottrun/spectra/spectra_K.txt"
filekaon2 = "../run/phenix_data/phenix_kaon.txt"

height=0.28
width=0.88

jpanel=0
ax = fig.add_axes([0.1,0.15,width,height])
myproton1 = np.loadtxt(fileproton1,skiprows=0,unpack=True)
myproton2 = np.loadtxt(fileproton2,skiprows=0,unpack=True)
plt.plot(myproton1[0],myproton1[1], 'o',label='protons')
plt.plot(myproton2[0]*1000,myproton2[3], 's')
plt.legend()

jpanel=1
ax = fig.add_axes([0.1,0.15+jpanel*height,width,height])
mykaon1 = np.loadtxt(filekaon1,skiprows=0,unpack=True)
mykaon2 = np.loadtxt(filekaon2,skiprows=0,unpack=True)
plt.plot(mykaon1[0],mykaon1[1], 'o',label='kaons')
plt.plot(mykaon2[0]*1000,mykaon2[3], 's')
ax.set_xticklabels([])
plt.legend()

jpanel=2
ax = fig.add_axes([0.1,0.15+jpanel*height,width,height])
mypion1 = np.loadtxt(filepion1,skiprows=0,unpack=True)
mypion2 = np.loadtxt(filepion2,skiprows=0,unpack=True)
plt.plot(mypion1[0],mypion1[1], 'o',label='pions')
plt.plot(mypion2[0]*1000,mypion2[3], 's')
ax.set_xticklabels([])
plt.legend()

plt.savefig('cf_outlongside.pdf',format='pdf')
os.system('open -a Preview cf_outlongside.pdf')

