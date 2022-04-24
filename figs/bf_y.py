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
height=(1.0-y0-0.02)/3.0

root2=sqrt(2.0)
xmin=0.0
xmax=1.8001

dNdY_pi=649 # both + and - multiplicity
dNdY_K=100.8
dNdY_p=61.4
dNdY_pi=0.5*dNdY_pi
dNdY_K=0.5*dNdY_K
dNdY_p=0.5*dNdY_p

#--------PIONS--------

results = np.loadtxt('../scottrun/results_bal/pipi/bf0_y.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_pipi=results[1]

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiplus=results[1]
Npi=results[2]
errordirect_pipi=results[3]*dNdY_pi
results = np.loadtxt('../scottrun/results_direct/pipluspiminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiminus=results[1]
Npi=Npi+results[2]
errordirect_pipi+=results[3]*dNdY_pi
bfdirect_pipi=dNdY_pi*(cfdirect_pipluspiminus-cfdirect_pipluspiplus)

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiplus=results[1]
errorallwfs_pipi=results[3]*dNdY_pi
results = np.loadtxt('../scottrun/results_allwfs/pipluspiminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiminus=results[1]
errorallwfs_pipi+=results[3]*dNdY_pi
bfallwfs_pipi=dNdY_pi*(cfallwfs_pipluspiminus-cfallwfs_pipluspiplus)

#--------KAONS--------

results = np.loadtxt('../scottrun/results_bal/KK/bf0_y.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_KK=results[1]

results = np.loadtxt('../scottrun/results_direct/KplusKplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKplus=results[1]
NK=results[2]
errordirect_K=results[3]*dNdY_K
results = np.loadtxt('../scottrun/results_direct/KplusKminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKminus=results[1]
NK=NK+results[2]
errordirect_K+=results[3]*dNdY_K
bfdirect_KK=dNdY_K*(cfdirect_KplusKminus-cfdirect_KplusKplus)

results = np.loadtxt('../scottrun/results_allwfs/KplusKplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKplus=results[1]
errorallwfs_K=results[3]*dNdY_K
results = np.loadtxt('../scottrun/results_allwfs/KplusKminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKminus=results[1]
errordirect_K+=results[3]*dNdY_K
bfallwfs_KK=dNdY_K*(cfallwfs_KplusKminus-cfallwfs_KplusKplus)

#--------PROTONS--------

results = np.loadtxt('../scottrun/results_bal/pp/bf0_y.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_pp=results[1]

results = np.loadtxt('../scottrun/results_direct/pp/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspplus=results[1]
Np=results[2]
errordirect_pp=results[3]*dNdY_p
results = np.loadtxt('../scottrun/results_direct/ppbar/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspminus=results[1]
Np=Np+results[2]
errordirect_pp+=results[3]*dNdY_p
bfdirect_pp=dNdY_p*(cfdirect_ppluspminus-cfdirect_ppluspplus)


results = np.loadtxt('../scottrun/results_allwfs/pp/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pp=results[1]
errorallwfs_pp=results[3]*dNdY_p
results = np.loadtxt('../scottrun/results_allwfs/ppbar/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_ppbar=results[1]
errorallwfs_pp+=results[3]*dNdY_p
bfallwfs_pp=dNdY_p*(cfallwfs_ppbar-cfallwfs_pp)

for jpanel in range(0,3):
   ax = fig.add_axes([x0,y0+jpanel*height,width,height])
   #ax.tick_params(axis='both', which='major', labelsize=14)
   plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')

   if jpanel==0:
      type='$\pi\pi$'
      plt.errorbar(x,bfallwfs_pipi+bfbal_pipi,errordirect_pipi,linestyle='-',linewidth=2,color='k',markersize=10,marker='o',label=type)
      plt.errorbar(x,bfallwfs_pipi,errorallwfs_pipi,linestyle='-',linewidth=2,color='b',markersize=10,marker='s',label=type)
      plt.plot(xbal,bfbal_pipi,linestyle='-',linewidth=4,color='r',label=type)
      ymin=-0.1
      ymax=0.8
      ax.set_yticks(np.arange(-1,1,0.2),minor=False)
      ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans',fontsize=24)
      ax.set_yticks(np.arange(-1,1,0.1),minor=True)
      plt.xlabel('$\Delta y$', fontsize=32, weight='normal')
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
      text(0.05,0.08,'final-state interactions',color='blue',fontsize='28',family='sans')
      text(0.05,0.49,'charge balance',color='red',fontsize='28',family='sans')
      text(0.35,0.72,'sum',color='k',fontsize='28',family='sans')
   if jpanel==1:
      type='$KK$'
      ymin=-0.05
      ymax=0.5
      plt.errorbar(x,bfallwfs_KK+bfbal_KK,errordirect_K,linestyle='-',linewidth=2,color='k',markersize=10,marker='o',label=type)
      plt.errorbar(x,bfallwfs_KK,errorallwfs_K,linestyle='-',linewidth=2,color='b',markersize=10,marker='s',label=type)
      plt.plot(xbal,bfbal_KK,linestyle='-',color='r',linewidth=4,label=type)
      ax.set_yticks(np.arange(-1,1,0.2),minor=False)
      ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans',fontsize=24)
      ax.set_yticks(np.arange(-1,1,0.1),minor=True)
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      plt.ylabel('$B(\Delta y)$',fontsize='32')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
   if jpanel==2:
     type='$pp$'
     ymin=-0.05
     ymax=0.5
     plt.errorbar(x,bfallwfs_pp+bfbal_pp,errordirect_pp,linestyle='-',linewidth=2,color='k',markersize=10,marker='o',label=type)
     plt.errorbar(x,bfallwfs_pp,errorallwfs_pp,linestyle='-',linewidth=2,color='b',markersize=10,marker='s',label=type)
     plt.plot(xbal,bfbal_pp,linestyle='-',color='r',linewidth=4,label=type)
     ax.set_yticks(np.arange(-1,1,0.2),minor=False)
     ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans',fontsize=24)
     ax.set_yticks(np.arange(-1,1,0.1),minor=True)
     plt.xlim(xmin,xmax)
     plt.ylim(ymin,ymax)
     plt.ylabel(None)
     ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
     
   ax.set_xticks(np.arange(xmin,xmax,0.2), minor=False)
   ax.set_xticks(np.arange(xmin,xmax,0.1), minor=True)
   if jpanel==0:
      ax.set_xticklabels(np.arange(xmin,xmax,0.2), minor=False, family='sans',fontsize=24)
      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
   else:
      ax.set_xticklabels([])

   if jpanel==0:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),"c) "+type,size=32,color='black',ha='right')
   if jpanel==1:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),"b) "+type,size=32,color='black',ha='right')
   if jpanel==2:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),"a) "+type,size=32,color='black',ha='right')
   if jpanel==3:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=32,color='black',ha='right') 
  
plt.savefig('bf_y.pdf',format='pdf')
os.system('open -a Preview bf_y.pdf')
#plt.show()
quit()
