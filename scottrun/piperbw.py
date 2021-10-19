#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

import numpy as np
import math
import scipy 
from scipy import optimize

par = np.empty(2, dtype=float)
#par[0]=1000
#par[1]=800
#par[2]=750
par[0]=98.5
par[1]=1.09

def chisquare(x0):
    subp=subprocess.Popen(['bwspectra',str(x0[0]),str(x0[1])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    stdout,stderror = subp.communicate()
    print(x0[0])
    print(x0[1])
#   print(x0[2])
#    print(x0[3])
#    print(x0[4])
    print(stdout)
    return float(stdout)

res = scipy.optimize.minimize(chisquare, par, method='powell', options={'disp': True,'return_all': True})
print(res.x)
quit()
