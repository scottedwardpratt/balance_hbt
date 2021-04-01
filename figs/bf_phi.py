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

results = np.loadtxt('results/bf_phi.dat',skiprows=0,unpack=True)
x=results[0]
bf_pipi=results[1]
bf_piK=results[2]
bf_pip=results[3]
bf_KK=results[4]
bf_Kp=results[5]
bf_pp=results[6]

for ipanel in range(0,2):
  for jpanel in range(0,3):
    ax = fig.add_axes([x0+ipanel*width,y0+jpanel*height,width,height])
    plt.plot([0,180],[0,0],linestyle='dashed',color='grey')
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
      ymin=-0.15
      ymax=0.05
      plt.plot(x,bf_pipi,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==0 and jpanel==1:
      type='$\pi K$'
      ymin=-0.05
      ymax=0.3
      plt.plot(x,bf_piK,linestyle='-',linewidth=3,color='r',label=type) 
    if ipanel==0 and jpanel==2:
      type='$\pi p$'
      ymin=-0.05
      ymax=0.16
      plt.plot(x,bf_pip,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==0:
      type='$KK$'
      ymin=-0.01
      ymax=0.02
      plt.plot(x,bf_KK,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==1:
      type='$Kp$'
      ymin=-0.02
      ymax=0.09
      plt.plot(x,bf_Kp,linestyle='-',linewidth=3,color='r',label=type)
    if ipanel==1 and jpanel==2:
      type='$pp$'
      ymin=-0.005
      ymax=0.01
      plt.plot(x,bf_pp,linestyle='-',linewidth=3,color='r',label=type)

    ax.tick_params(axis='both', which='major', labelsize=14)

    if ipanel == 0:
      ax.set_xticks(np.arange(0,170,30), minor=False)
      ax.set_xticks(np.arange(0,200,10), minor=True)
    else:
      ax.set_xticks(np.arange(0,200,30), minor=False)
      ax.set_xticks(np.arange(0,200,10), minor=True)
    if jpanel==0:
      if ipanel==0:
        ax.set_xticklabels(np.arange(0,170,30), minor=False, family='serif')
      else:
        ax.set_xticklabels(np.arange(0,200,30), minor=False, family='serif')
    else:
      ax.set_xticklabels([])
    plt.xlim(0.0,180)
  
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    
    if ipanel==0 and jpanel==0:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-0.5,0.5,0.05), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.01), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.5,0.05), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==0 and jpanel==1:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-0.5,0.29,0.1), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.02), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.29,0.1), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==0 and jpanel==2:
      ax.yaxis.tick_left()
      ax.set_yticks(np.arange(-0.5,0.5,0.05), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.01), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.5,0.05), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==0:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-0.5,0.5,0.01), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.002), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.5,0.01), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==1:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-0.5,0.5,0.05), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.01), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.5,0.05), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    if ipanel==1 and jpanel==2:
      ax.yaxis.tick_right()
      ax.set_yticks(np.arange(-0.5,0.5,0.005), minor=False)
      ax.set_yticks(np.arange(-0.5,0.5,0.001), minor=True)
      ax.set_yticklabels(np.arange(-0.5,0.5,0.005), minor=False, family='serif')
      ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%3.2f'))
    

    plt.ylim(ymin,ymax)
    
    if ipanel==0 and jpanel==0:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==0 and jpanel==1:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==0 and jpanel==2:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==0:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==1:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')
    if ipanel==1 and jpanel==2:
      text(175,ymin+0.85*(ymax-ymin),type,size=24,color='black',ha='right')

    if jpanel==0:
      plt.xlabel('$\Delta\phi$', fontsize=24, weight='normal')
    else:
      plt.xlabel(None)
    if jpanel==1 and ipanel ==0:
      plt.ylabel('$B_{HBT}(\Delta\phi)$',fontsize=24)
    else:
      plt.ylabel(None)

  
plt.savefig('bf_phi.pdf',format='pdf')
os.system('open -a Preview bf_phi.pdf')
#plt.show()
quit()
