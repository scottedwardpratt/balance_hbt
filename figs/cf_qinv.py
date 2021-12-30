import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-4,3))
#plt.ticklabel_format(style='sci', axis='y')

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(10,15))
fig = plt.figure(1)

x0=0.15
width=(0.96-x0)
y0=0.07
height=0.45

root2=sqrt(2.0)
xmin=0.0
xmax=100

dNdY_pi=649 # both + and - multiplicity
dNdY_K=100.8
dNdY_p=61.4
dNdY_pi=0.5*dNdY_pi
dNdY_K=0.5*dNdY_K
dNdY_p=0.5*dNdY_p

#--------PIONS--------

results = np.loadtxt('../scottrun/results_bal/pipi/bf0_outsidelong.dat',skiprows=1,unpack=True)
xbal=results[0]
cfbal_pipi=results[7]*dNdY_pi

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiplus=results[10]
results = np.loadtxt('../scottrun/results_direct/pipluspiminus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiminus=results[10]
cfdirect_pipi=cfdirect_pipluspiminus-cfdirect_pipluspiplus

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiplus=results[10]
results = np.loadtxt('../scottrun/results_allwfs/pipluspiminus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiminus=results[10]
cfallwfs_pipi=cfallwfs_pipluspiminus-cfallwfs_pipluspiplus

#-------- Balance Distortions --------

results = np.loadtxt('../scottrun/results_bal/pipi_denom/bf0_outsidelong.dat',skiprows=1,unpack=True)
xbal=results[0]
cfbal_denom=results[7]*dNdY_pi
etawidth=1.8
cfbal_pipi=(cfbal_pipi/cfbal_denom)/(etawidth*sqrt(2.0*np.pi)*dNdY_pi)


ax = fig.add_axes([x0,y0,width,height])
#ax.tick_params(axis='both', which='major', labelsize=14)
plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')

type='$\pi\pi$'
#plt.plot(x,cfallwfs_pipi+cfbal_pipi,linestyle='-',linewidth2,color='r',markersize=10,marker='o',label=type)
plt.plot(2*x,cfallwfs_pipi,linestyle='-',linewidth=1,color='b',markersize=10,marker='s',label=type)
plt.plot(2*xbal,cfbal_pipi,linestyle='-',linewidth=3,color='r',label=type)
plt.plot(2*xbal,cfallwfs_pipi+cfbal_pipi,linestyle='-',linewidth=3,color='k',label=type)
ymin=-0.2
ymax=0.8
ax.set_yticks(np.arange(-1,1,0.2),minor=False)
ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans',size=24)
ax.set_yticks(np.arange(-1,1,0.1),minor=True)
plt.ylabel('$C_{\pi^+\pi^-}-C_{\pi^+\pi^+}$',fontsize=32)
plt.xlabel('$Q_{\\rm inv}$ [MeV/$c$]', fontsize=32,family='sans')
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
text(45,-0.07,'final-state interactions',color='blue',fontsize='28',family='sans')
text(12,0.02,'charge balance',color='red',fontsize='28',family='sans')
text(5,0.6,'sum',color='black',fontsize='28',family='sans')
ax.set_xticks(np.arange(xmin,xmax+1,50), minor=False)
ax.set_xticks(np.arange(xmin,xmax+1,10), minor=True)
ax.set_xticklabels(np.arange(xmin,xmax+1,50), minor=False, family='sans',size=24)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
text(93,0.72,'b)',color='black',fontsize=32,family='sans')


#upper panel
ax = fig.add_axes([x0,y0+height,width,height])

plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')

type='$\pi\pi$'
#plt.plot(x,cfallwfs_pipi+cfbal_pipi,linestyle='-',linewidth2,color='r',markersize=10,marker='o',label=type)
plt.plot(2*x,cfallwfs_pipi,linestyle='-',linewidth=1,color='b',markersize=10,marker='s',label=type)
plt.plot(2*xbal,cfbal_pipi,linestyle='-',linewidth=3,color='r',label=type)
plt.plot(2*xbal,cfallwfs_pipi+cfbal_pipi,linestyle='-',linewidth=3,color='k',label=type)
ymin=-0.05
ymax=0.02
ax.set_yticks(np.arange(-1,1,0.02),minor=False)
ax.set_yticklabels(np.arange(-1,1,0.02),minor=False,family='sans',size=24)
ax.set_yticks(np.arange(-1,1,0.01),minor=True)
plt.ylabel('$C_{\pi^+\pi^-}-C_{\pi^+\pi^+}$',fontsize=32)
plt.xlabel(None, fontsize=32)
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
text(45,-0.014,'final-state interactions',color='blue',fontsize='28',family='sans')
text(12,0.004,'charge balance',color='red',fontsize='28',family='sans')
text(65,0.007,'sum',color='black',fontsize='28',family='sans')
ax.set_xticks(np.arange(xmin,xmax+1,50), minor=False)
ax.set_xticks(np.arange(xmin,xmax+1,10), minor=True)
ax.set_xticklabels(np.arange(xmin,xmax+1,50), minor=False, family='sans',size=16)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

text(93,-0.047,'a)',color='black',fontsize=32,family='sans')

   
ax.set_xticks(np.arange(xmin,xmax+1,50), minor=False)
ax.set_xticks(np.arange(xmin,xmax+1,10), minor=True)
ax.set_xticklabels([], minor=False)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
  
plt.savefig('cf_qinv.pdf',format='pdf')
os.system('open -a Preview cf_qinv.pdf')
#plt.show()
quit()
