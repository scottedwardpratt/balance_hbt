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

x0=0.14
width=(0.97-x0)
y0=0.1
height=(1.0-y0-0.04)/3.0

root2=sqrt(2.0)
xmin=0.0
xmax=180.001

dNdY_pi=649 # both + and - multiplicity
dNdY_K=100.8
dNdY_p=61.4
dNdY_pi=0.5*dNdY_pi
dNdY_K=0.5*dNdY_K
dNdY_p=0.5*dNdY_p

dNdphi_pi=dNdY_pi*1.8/(2.0*np.pi)
dNdphi_K=dNdY_K*1.8/(2.0*np.pi)
dNdphi_p=dNdY_p*1.8/(2.0*np.pi)

#--------PIONS--------

results = np.loadtxt('../scottrun/results_bal/pipi/bf0_phi.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_pipi=results[1]

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiplus=results[1]
Npi=results[2]
errordirect_pipi=results[3]*dNdphi_pi
results = np.loadtxt('../scottrun/results_direct/pipluspiminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiminus=results[1]
Npi=Npi+results[2]
errordirect_pipi+=results[3]*dNdphi_pi
bfdirect_pipi=dNdphi_pi*(cfdirect_pipluspiminus-cfdirect_pipluspiplus)

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiplus=results[1]
errorallwfs_pipi=results[3]*dNdphi_pi
results = np.loadtxt('../scottrun/results_allwfs/pipluspiminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiminus=results[1]
errorallwfs_pipi+=results[3]*dNdphi_pi
bfallwfs_pipi=dNdphi_pi*(cfallwfs_pipluspiminus-cfallwfs_pipluspiplus)

#--------KAONS--------

results = np.loadtxt('../scottrun/results_bal/pipi/bf0_phi.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_KK=results[1]

results = np.loadtxt('../scottrun/results_direct/KplusKplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKplus=results[1]
NK=results[2]
errordirect_K=results[3]*dNdphi_K
results = np.loadtxt('../scottrun/results_direct/KplusKminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKminus=results[1]
NK=NK+results[2]
errordirect_K+=results[3]*dNdphi_K
bfdirect_KK=dNdphi_K*(cfdirect_KplusKminus-cfdirect_KplusKplus)

results = np.loadtxt('../scottrun/results_allwfs/KplusKplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKplus=results[1]
errorallwfs_K=results[3]*dNdphi_K
results = np.loadtxt('../scottrun/results_allwfs/KplusKminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKminus=results[1]
errordirect_K+=results[3]*dNdphi_K
bfallwfs_KK=dNdphi_K*(cfallwfs_KplusKminus-cfallwfs_KplusKplus)

#--------PROTONS--------

results = np.loadtxt('../scottrun/results_bal/pp/bf0_phi.dat',skiprows=1,unpack=True)
xbal=results[0]
bfbal_pp=results[1]

results = np.loadtxt('../scottrun/results_direct/pp/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspplus=results[1]
Np=results[2]
errordirect_pp=results[3]*dNdphi_p
results = np.loadtxt('../scottrun/results_direct/ppbar/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspminus=results[1]
Np=Np+results[2]
errordirect_pp+=results[3]*dNdphi_p
bfdirect_pp=dNdphi_p*(cfdirect_ppluspminus-cfdirect_ppluspplus)


results = np.loadtxt('../scottrun/results_allwfs/pp/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pp=results[1]
errorallwfs_pp=results[3]*dNdphi_p
results = np.loadtxt('../scottrun/results_allwfs/ppbar/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_ppbar=results[1]
errorallwfs_pp+=results[3]*dNdphi_p
bfallwfs_pp=dNdphi_p*(cfallwfs_ppbar-cfallwfs_pp)

#bfbal_pipi=bfbal_pipi*180.0/np.pi
#bfbal_KK=bfbal_KK*180.0/np.pi
#bfbal_pp=bfbal_pp*180.0/np.pi

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
      ymax=0.65
      ax.set_yticks(np.arange(-1,1,0.2),minor=False)
      ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans')
      ax.set_yticks(np.arange(-1,1,0.1),minor=True)
      plt.xlabel('$\Delta\phi$', fontsize=28, weight='normal')
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
      text(5,0.05,'final-state interactions',color='blue',fontsize='22',family='sans')
      text(5,0.32,'charge balance',color='red',fontsize='22',family='sans')
      text(30,0.52,'sum',color='black',fontsize='22',family='sans')
   if jpanel==1:
      type='$KK$'
      ymin=-0.05
      ymax=0.65
      plt.errorbar(x,bfallwfs_KK+bfbal_KK,errordirect_K,linestyle='-',linewidth=2,color='k',markersize=10,marker='o',label=type)
      plt.errorbar(x,bfallwfs_KK,errorallwfs_K,linestyle='-',linewidth=2,color='b',markersize=10,marker='s',label=type)
      plt.plot(xbal,bfbal_KK,linestyle='-',linewidth=4,color='r',label=type)
      ax.set_yticks(np.arange(-1,1,0.2),minor=False)
      ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans')
      ax.set_yticks(np.arange(-1,1,0.1),minor=True)
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      plt.ylabel('$B(\Delta\phi)$',fontsize='24')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
   if jpanel==2:
     type='$pp$'
     ymin=-0.05
     ymax=0.45
     plt.errorbar(x,bfallwfs_pp+bfbal_pp,errordirect_pp,linestyle='-',linewidth=2,color='k',markersize=10,marker='o',label=type)
     plt.errorbar(x,bfallwfs_pp,errorallwfs_pp,linestyle='-',linewidth=2,color='b',markersize=10,marker='s',label=type)
     plt.plot(xbal,bfbal_pp,linestyle='-',linewidth=4,color='r',label=type)
     ax.set_yticks(np.arange(-1,1,0.2),minor=False)
     ax.set_yticklabels(np.arange(-1,1,0.2),minor=False,family='sans')
     ax.set_yticks(np.arange(-1,1,0.1),minor=True)
     plt.xlim(xmin,xmax)
     plt.ylim(ymin,ymax)
     plt.ylabel(None)
     ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
     
   ax.set_xticks(np.arange(xmin,xmax,30), minor=False)
   ax.set_xticks(np.arange(xmin,xmax,10), minor=True)
   if jpanel==0:
      ax.set_xticklabels(np.arange(xmin,xmax,30), minor=False, family='sans')
      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
   else:
      ax.set_xticklabels([])

   if jpanel==0:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
   if jpanel==1:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
   if jpanel==2:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
   if jpanel==3:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right') 
  
plt.savefig('bf_phi.pdf',format='pdf')
os.system('open -a Preview bf_phi.pdf')
#plt.show()
quit()
