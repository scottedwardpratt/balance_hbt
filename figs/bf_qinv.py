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
width=0.5*(1.0-x0-0.1)
y0=0.08
height=(1.0-y0-0.04)/3.0
xmax=500

results = np.loadtxt('results/bf_qinv.dat',skiprows=0,unpack=True)
x=results[0]
bf_pipi=1000*results[1]
bf_piK=1000*results[2]
bf_pip=1000*results[3]
bf_KK=1000*results[4]
bf_Kp=1000*results[5]
bf_pp=1000*results[6]

for ipanel in range(0,2):
  for jpanel in range(0,3):
    ax = fig.add_axes([x0+ipanel*width,y0+jpanel*height,width,height])
    plt.plot([0,500],[0,0],linestyle='dashed',color='grey')
    if jpanel==0:
      type='pipi'
    if jpanel==1:
      type='KK'
    if jpanel==2:
      type='pp'

    #####################################
    # cascade
    
    
    if ipanel==0 and jpanel==0:
      type='$\pi\pi$'
      ymin=-2.5
      ymax=1
      plt.plot(x,bf_pipi,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==0 and jpanel==1:
      type='$\pi K$'
      ymin=-1
      ymax=1
      plt.plot(x,bf_piK,linestyle='-',linewidth=3,color='r',label=type) 
    if ipanel==0 and jpanel==2:
      type='$\pi p$'
      ymin=-0.5
      ymax=0.5
      plt.plot(x,bf_pip,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==0:
      type='$KK$'
      ymin=-0.2
      ymax=0.2
      plt.plot(x,bf_KK,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==1:
      type='$Kp$'
      ymin=-0.05
      ymax=0.3
      plt.plot(x,bf_Kp,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==2:
      type='$pp$'
      ymin=-0.01
      ymax=0.04
      plt.plot(x,bf_pp,linestyle='-',linewidth=3,color='r',label=type)

    ax.tick_params(axis='both', which='major', labelsize=14)

    ax.set_xticks(np.arange(0,500,100), minor=False)
    ax.set_xticks(np.arange(0,500,25), minor=True)
    if jpanel==0:
      ax.set_xticklabels(np.arange(0,500,100), minor=False, family='serif')
    else:
      ax.set_xticklabels([])
    plt.xlim(0.0,xmax)
  
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

    if ipanel==0 and jpanel==0:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-5,5,1), minor=False)
      ax.set_yticks(np.arange(-5,5,0.25), minor=True)
      ax.set_yticklabels(np.arange(-5,5,1), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==0 and jpanel==1:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-5,5,0.5), minor=False)
      ax.set_yticks(np.arange(-5,5,0.1), minor=True)
      ax.set_yticklabels(np.arange(-5,5,0.5), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==0 and jpanel==2:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-5,5,0.2), minor=False)
      ax.set_yticks(np.arange(-5,5,0.05), minor=True)
      ax.set_yticklabels(np.arange(-5,5,0.2), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==0:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-5,5,0.1), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.025), minor=True)
      ax.set_yticklabels(np.arange(-5,5,0.1), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==1:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-5,5,0.1), minor=False)
      ax.set_yticks(np.arange(-5,5,0.025), minor=True)
      ax.set_yticklabels(np.arange(-5,5,0.1), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==2:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-5,5,0.02), minor=False)
      ax.set_yticks(np.arange(-5,5,0.005), minor=True)
      ax.set_yticklabels(np.arange(-5,5,0.02), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))

    plt.ylim(ymin,ymax)
    
    if ipanel==0 and jpanel==0:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==0 and jpanel==1:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==0 and jpanel==2:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==0:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==1:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==2:
      text(0.95*xmax,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    

    if jpanel==0:
      plt.xlabel('$q_{\\rm inv}$', fontsize=24, weight='normal')
    else:
      plt.xlabel(None)
    if jpanel==1 and ipanel ==0:
      plt.ylabel('$B_{HBT}(q_{\\rm inv})$  (GeV/$c$)$^{-1}$',fontsize=24)
    else:
      plt.ylabel(None)

  
plt.savefig('bf_qinv.pdf',format='pdf')
os.system('open -a Preview bf_qinv.pdf')
#plt.show()
quit()
