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
plt.figure(figsize=(6,9))
fig = plt.figure(1)

x0=0.12
width=0.5*(1.0-x0-0.1)
y0=0.08
height=(1.0-y0-0.04)/3.0

results = np.loadtxt('cfresults_pp_allwfs/pp/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
q=results[0]
cfout=results[1]
cfside=results[3]
cflong=results[5]
cfqinv=results[7]

ymin=-0.5
ymax=0.1
xmin=0
xmax=100.0
height=0.43

for jpanel in range(0,4):
  if jpanel < 2:
    ax = fig.add_axes([0.08,0.08+jpanel*height,0.90,height])
  plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')
  if jpanel==1:
    type='out'
  if jpanel==2:
    type='side'
  if jpanel==3:
    type='long'
  if jpanel==0:
    type='qinv'

  if jpanel==0:
    plt.plot(q,cfqinv,linestyle='-',linewidth=3,color='k',markersize=6,marker='o',label=type)
    #plt.plot(q,cfall_qinv,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
  if jpanel==1:
    plt.plot(q,cfout,linestyle='-',linewidth=3,color='r',markersize=6,marker='o',label=type)
    #plt.plot(q,cfall_out,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
  if jpanel==2:
    plt.plot(q,cfside,linestyle='-',linewidth=3,color='g',markersize=6,marker='o',label=type)
    #plt.plot(q,cfall_side,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
  if jpanel==3:
    plt.plot(q,cflong,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)
    #plt.plot(q,cfall_long,linestyle='-',linewidth=3,color='b',markersize=6,marker='o',label=type)

  ax.tick_params(axis='both', which='major', labelsize=14)

  ax.set_xticks(np.arange(xmin,xmax,40), minor=False)
  ax.set_xticks(np.arange(xmin,xmax,20), minor=True)
  if jpanel==0:
    ax.set_xticklabels(np.arange(xmin,xmax,40), minor=False, family='serif')
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
    plt.xlabel('$q$', fontsize=24, weight='normal')
  else:
    plt.xlabel(None)
  if jpanel==1:
    plt.ylabel('$C(q)-1.0$',fontsize=24)
  else:
    plt.ylabel(None)
  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)

  
plt.savefig('cf2panels_outsidelong.pdf',format='pdf')
os.system('open -a Preview cf2panels_outsidelong.pdf')
#plt.show()
quit()
