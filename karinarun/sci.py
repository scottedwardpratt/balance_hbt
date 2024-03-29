#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

import numpy as np
import math
import scipy 
from scipy import optimize


x0 = np.empty(3, dtype=float)
x0[0]=0.40754638 #0.4173118     #A
x0[1]=10.98130514 #14.79443084    #tau
x0[2]=6.74439321 #12.33820392    #R


def chisquare(x0):
    subp=subprocess.Popen(['chisquare', str(x0[0]), str(x0[1]), str(x0[2])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    stdout,stderror = subp.communicate()
    split = stdout.split()
    length=len(split)
    print(x0[0])
    print(x0[1])
    print(x0[2])
    print(split[length-1])
    return float(split[length-1])
    
res = scipy.optimize.minimize(chisquare, x0, method='powell', options={'disp': True, 'maxiter': 4, 'maxfev': 100})
print(res.x)
quit()