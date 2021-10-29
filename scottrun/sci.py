#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

import numpy as np
import math
import scipy 
import os
from scipy import optimize

os.system('rm -f logfiles/sci.txt')
logfile=open("logfiles/sci.txt","a")

x0 = np.empty(2, dtype=float)
x0[0]=12.185068 #14.79443084    #tau
x0[1]=13.815719 #12.33820392    #R

def chisquare(x0):
   subp=subprocess.Popen(['chisquare', str(x0[0]), str(x0[1])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
   stdout,stderror = subp.communicate()
   split = stdout.split()
   length=len(split)-1
   chi2=float(split[length])
   A=float(split[length-1])
   logfile.write("tau=%f, " % x0[0])
   logfile.write("R=%f, " % x0[1])
   logfile.write("A=%f, " % A)
   logfile.write("chi^2=%f\n" % chi2)
   logfile.flush()
   return float(chi2)

    
res = scipy.optimize.minimize(chisquare, x0, method='powell', options={'disp': True, 'maxiter': 3, 'maxfev': 40})
logfile.write("tau=%f," % res.x[0])
logfile.write("R=%f\n" % res.x[1])
logfile.write("sci.py calculating best values\n")
x0[0]=res.x[0]
x0[1]=res.x[1]
chisquare(x0)
logfile.flush()
logfile.close()
quit()