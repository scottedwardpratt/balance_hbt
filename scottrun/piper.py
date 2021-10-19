#!/usr/bin/env python
import subprocess
from subprocess import Popen, PIPE

x=str(3)
y=str(11)

subp=subprocess.Popen(['crap',x,y], stdout=PIPE,stderr=PIPE,universal_newlines=True)
stdout,stderror = subp.communicate()
print(stdout)
quit()
