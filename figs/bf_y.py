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

xmin=0.0
xmax=1.8

dNdY_pi=1000.0
dNdY_K=100.0
dNdY_p=50.0

#--------PIONS--------

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiplus=results[1]
Npi=results[2]
results = np.loadtxt('../scottrun/results_direct/pipluspiminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiminus=results[1]
Npi=Npi+results[2]
bfdirect_pipi=dNdY_pi*(cfdirect_pipluspiminus-cfdirect_pipluspiplus)
error_pipi=results[3]*dNdY_pi

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiplus=results[1]
results = np.loadtxt('../scottrun/results_allwfs/pipluspiminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiminus=results[1]
bfallwfs_pipi=dNdY_pi*(cfallwfs_pipluspiminus-cfallwfs_pipluspiplus)

#--------KAONS--------

results = np.loadtxt('../scottrun/results_direct/KplusKplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKplus=results[1]
NK=results[2]
results = np.loadtxt('../scottrun/results_direct/KplusKminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKminus=results[1]
NK=NK+results[2]
bfdirect_KK=dNdY_K*(cfdirect_KplusKminus-cfdirect_KplusKplus)
error_K=results[3]*dNdY_K

results = np.loadtxt('../scottrun/results_allwfs/KplusKplus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKplus=results[1]
results = np.loadtxt('../scottrun/results_allwfs/KplusKminus/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKminus=results[1]
bfallwfs_KK=dNdY_K*(cfallwfs_KplusKminus-cfallwfs_KplusKplus)

#--------PROTONS--------

results = np.loadtxt('../scottrun/results_direct/pp/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspplus=results[1]
Np=results[2]
results = np.loadtxt('../scottrun/results_direct/ppbar/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspminus=results[1]
Np=Np+results[2]
bfdirect_pp=dNdY_p*(cfdirect_ppluspminus-cfdirect_ppluspplus)
error_p=results[3]*dNdY_p

results = np.loadtxt('../scottrun/results_allwfs/pp/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pp=results[1]
results = np.loadtxt('../scottrun/results_allwfs/ppbar/cf_y.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_ppbar=results[1]
bfallwfs_pp=dNdY_p*(cfallwfs_ppbar-cfallwfs_pp)

for jpanel in range(0,3):
   ax = fig.add_axes([x0,y0+jpanel*height,width,height])
   #ax.tick_params(axis='both', which='major', labelsize=14)
   plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')

   if jpanel==0:
      type='$\pi\pi$'
      #plt.plot(x,bfdirect_pipi,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
      plt.errorbar(x,bfdirect_pipi,error_pipi,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
      #plt.plot(x,bfallwfs_pipi,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
      plt.errorbar(x,bfallwfs_pipi,error_pipi,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
      ymin=-0.3
      ymax=0.3
      ax.set_yticks(np.arange(-0.2,0.4,0.05),minor=False)
      ax.set_yticklabels(np.arange(-0.2,0.4,0.05),minor=False)
      ax.set_yticks(np.arange(-0.2,0.4,0.01),minor=True)
      plt.xlabel('$\Delta y$', fontsize=24, weight='normal')
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
   if jpanel==1:
      type='$KK$'
      ymin=-0.05
      ymax=0.05
      #plt.plot(x,bfdirect_KK,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
      plt.errorbar(x,bfdirect_KK,error_K,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
      #plt.plot(x,bfallwfs_KK,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
      plt.errorbar(x,bfallwfs_KK,error_K,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
      ax.set_yticks(np.arange(-0.1,0.2,0.01),minor=False)
      ax.set_yticklabels(np.arange(-0.10,0.2,0.01),minor=False)
      ax.set_yticks(np.arange(-0.10,0.2,0.005),minor=True)
      plt.xlim(xmin,xmax)
      plt.ylim(ymin,ymax)
      plt.ylabel('$B(y)$')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
   if jpanel==2:
     type='$pp$'
     ymin=-0.05
     ymax=0.05
     #plt.plot(x,bfdirect_pp,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
     plt.errorbar(x,bfdirect_pp,error_p,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
     #plt.plot(x,bfallwfs_pp,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
     plt.errorbar(x,bfallwfs_pp,error_p,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
     ax.set_yticks(np.arange(-0.1,0.2,0.01),minor=False)
     ax.set_yticklabels(np.arange(-0.10,0.2,0.01),minor=False)
     ax.set_yticks(np.arange(-0.10,0.2,0.005),minor=True)
     plt.xlim(xmin,xmax)
     plt.ylim(ymin,ymax)
     plt.ylabel(None)
     ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
     
   ax.set_xticks(np.arange(xmin,xmax,0.2), minor=False)
   ax.set_xticks(np.arange(xmin,xmax,0.1), minor=True)
   if jpanel==0:
      ax.set_xticklabels(np.arange(xmin,xmax,0.2), minor=False, family='serif')
      ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
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
  
plt.savefig('bf_y.pdf',format='pdf')
os.system('open -a Preview bf_y.pdf')
#plt.show()
quit()
