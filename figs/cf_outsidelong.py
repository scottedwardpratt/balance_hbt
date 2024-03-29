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

x0=0.2
width=(1.0-x0-0.1)
y0=0.08
height=(1.0-y0-0.04)/3.0

results = np.loadtxt('../scottrun/results_direct/pipluspiplus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
q=results[0]
cf1wf_out=results[1]
cf1wf_side=results[4]
cf1wf_long=results[7]
cf1wf_qinv=results[10]

results = np.loadtxt('../scottrun/results_allwfs/pipluspiplus/cf_outsidelong.dat',skiprows=0,unpack=True)
x=results[0]
q=results[0]
cfall_out=results[1]
cfall_side=results[4]
cfall_long=results[7]
cfall_qinv=results[10]

expresults = np.loadtxt('cf_exp.dat',skiprows=0,unpack=True)
qexp=0.5*expresults[0]*1000.0
cfexp_out=expresults[1]-1.0
cfexp_side=expresults[2]-1.0
cfexp_long=expresults[3]-1.0

ymin=-0.02
ymax=0.13
xmin=0
xmax=100.0
height=0.224

A=0.615

for jpanel in range(0,4):
  ax = fig.add_axes([x0,y0+jpanel*height,width,height])
  plt.plot([xmin,xmax],[0,0],linestyle='dashed',color='grey')
  if jpanel==0:
    type='out'
  if jpanel==1:
    type='side'
  if jpanel==2:
    type='long'
  if jpanel==3:
    type='qinv'
    ymax=0.23
    ymin=-0.23

  if jpanel==0:
    plt.plot(q,A*cf1wf_out,linestyle='-',linewidth=2,color='r',markersize=6,marker='o',label=type)
    plt.plot(q,A*cfall_out,linestyle='-',linewidth=2,color='b',markersize=6,marker='s',label=type)
    plt.plot(qexp,cfexp_out,linestyle='-',linewidth=2,color='g',markersize=6,marker='o',label=type)
  if jpanel==1:
    plt.plot(q,A*cf1wf_side,linestyle='-',linewidth=2,color='r',markersize=6,marker='o',label=type)
    plt.plot(q,A*cfall_side,linestyle='-',linewidth=2,color='b',markersize=6,marker='s',label=type)
    plt.plot(qexp,cfexp_side,linestyle='-',linewidth=2,color='g',markersize=6,marker='o',label=type)
  if jpanel==2:
    plt.plot(q,A*cf1wf_long,linestyle='-',linewidth=2,color='r',markersize=6,marker='o',label=type)
    plt.plot(q,A*cfall_long,linestyle='-',linewidth=2,color='b',markersize=6,marker='s',label=type)
    plt.plot(qexp,cfexp_long,linestyle='-',linewidth=2,color='g',markersize=6,marker='o',label=type)
  if jpanel==3:
    plt.plot(q,cf1wf_qinv,linestyle='-',linewidth=2,color='r',markersize=6,marker='o',label=type)
    plt.plot(q,cfall_qinv,linestyle='-',linewidth=2,color='b',markersize=6,marker='s',label=type)

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

#plt.legend()
  
plt.savefig('cf_outlongside.pdf',format='pdf')
os.system('open -a Preview cf_outlongside.pdf')
#plt.show()
quit()
