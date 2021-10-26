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
x0[0]=10.98130514 #14.79443084    #tau
x0[1]=6.74439321 #12.33820392    #R

def chisquare(x0):
   subp=subprocess.Popen(['chisquare', str(x0[0]), str(x0[1])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
   stdout,stderror = subp.communicate()
   split = stdout.split()
   length=len(split)-1
   logfile.write("tau=%f, " % x0[0])
   logfile.write("R=%f, " % x0[1])
   logfile.write("chi^2=%f\n" % float(split[length]))
   logfile.flush()
   return float(split[length])

    
res = scipy.optimize.minimize(chisquare, x0, method='powell', options={'disp': True, 'maxiter': 100, 'maxfev': 100})
logfile.write("best values=")
logfile.write(res.x)
logfile.close()
quit()