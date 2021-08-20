#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

import numpy as np
import math
import scipy 
from scipy import optimize

par = np.empty(5, dtype=float)
par[0]=386
par[1]=260
par[2]=170
par[3]=86
par[4]=1.2

def chisquare(x0):
    subp=subprocess.Popen(['bwspectra',str(x0[0]),str(x0[1]),str(x0[2]), str(x0[3]), str(x0[4])], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    stdout,stderror = subp.communicate()
    return float(stdout)

res = scipy.optimize.minimize(chisquare, par, method='nelder-mead')
print(res.x)
quit()
