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

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(9,9))
fig = plt.figure(1)

x0=0.12
width=(1.0-x0-0.1)
y0=0.08
height=(1.0-y0-0.04)/3.0

xmin=0.0
xmax=180.0
ymin=-0.1
ymax=0.2

dNdY_pi=500.0
dNdY_K=100.0
dNdY_p=50.0

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiplus=results[1]
results = np.loadtxt('../scottrun/results_direct/pipluspiminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_pipluspiminus=results[1]
bfdirect_pipi=dNdY_pi*(cfdirect_pipluspiminus-cfdirect_pipluspiplus)

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiplus=results[1]
results = np.loadtxt('../scottrun/results_allwfs/pipluspiminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pipluspiminus=results[1]
bfallwfs_pipi=dNdY_pi*(cfallwfs_pipluspiminus-cfallwfs_pipluspiplus)

results = np.loadtxt('../scottrun/results_direct/KplusKplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKplus=results[1]
results = np.loadtxt('../scottrun/results_direct/KplusKminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_KplusKminus=results[1]
bfdirect_KK=dNdY_K*(cfdirect_KplusKminus-cfdirect_KplusKplus)

results = np.loadtxt('../scottrun/results_allwfs/KplusKplus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKplus=results[1]
results = np.loadtxt('../scottrun/results_allwfs/KplusKminus/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_KplusKminus=results[1]
bfallwfs_KK=dNdY_K*(cfallwfs_KplusKminus-cfallwfs_KplusKplus)

results = np.loadtxt('../scottrun/results_direct/pp/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspplus=results[1]
results = np.loadtxt('../scottrun/results_direct/ppbar/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfdirect_ppluspminus=results[1]
bfdirect_pp=dNdY_p*(cfdirect_ppluspminus-cfdirect_ppluspplus)

results = np.loadtxt('../scottrun/results_allwfs/pp/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_pp=results[1]
results = np.loadtxt('../scottrun/results_allwfs/ppbar/cf_phi.dat',skiprows=0,unpack=True)
x=results[0]
cfallwfs_ppbar=results[1]
bfallwfs_pp=dNdY_p*(cfallwfs_ppbar-cfallwfs_pp)


for jpanel in range(0,3):
  ax = fig.add_axes([x0,y0+jpanel*height,width,height])
  plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')
  if jpanel==0:
    type='$\pi\pi$'
  if jpanel==1:
    type='$KK$'
  if jpanel==2:
    type='$pp$'
    
  if jpanel==0:
    plt.plot(x,bfdirect_pipi,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
    plt.plot(x,bfallwfs_pipi,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
  if jpanel==1:
    plt.plot(x,100*bfdirect_KK,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
    plt.plot(x,100*bfallwfs_KK,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
  if jpanel==2:
    plt.plot(x,100*bfdirect_pp,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
    plt.plot(x,100*bfallwfs_pp,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)

  ax.tick_params(axis='both', which='major', labelsize=14)

  ax.set_xticks(np.arange(xmin,xmax+0.01,60), minor=False)
  ax.set_xticks(np.arange(xmin,xmax,30), minor=True)
  if jpanel==0:
    ax.set_xticklabels(np.arange(xmin,xmax+0.01,60), minor=False, family='serif')
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

  if jpanel==0:
    plt.xlabel('$\Delta y$', fontsize=24, weight='normal')
  else:
    plt.xlabel(None)
  if jpanel==1:
    plt.ylabel('$B(\Delta y)$',fontsize=24)
  else:
    plt.ylabel(None)
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)

plt.savefig('bf_phi.pdf',format='pdf')
os.system('open -a Preview bf_phi.pdf')
#plt.show()
quit()
