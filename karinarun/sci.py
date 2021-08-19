#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

import numpy as np
import math
import scipy 
from scipy import optimize

x0 = np.empty(3, dtype=float)
x0[0]=2     #A
x0[1]=10    #tau
x0[2]=10    #R

def chisquare(x0):
    subp=subprocess.Popen(['chisquare', str(x0[0]), str(x0[1]), str(x0[2])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    stdout,stderror = subp.communicate()
    split = stdout.split()
    length=len(split)
    print(split[length-1])
    return float(split[length-1])
    
res = scipy.optimize.minimize(chisquare, x0, method='nelder-mead', options={'disp':True})
print(res.x)
print(res.nit)
quit()