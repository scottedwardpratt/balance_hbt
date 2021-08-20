import subprocess
import multiprocessing
import os
import numpy as np
from subprocess import Popen, PIPE

#subprocess.call(['sh', './runner_cf.sh'])

par = np.arange(12)

def balhbt(x1,x2,x3):
    subp=subprocess.Popen(['balhbt',str(x1), str(x2), str(x3)], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    print("ID of process running balhbt: {}".format(os.getpid()))
    stdout,stderror = subp.communicate()
    return float(stdout)
        
def summer():
    subp1=subprocess.Popen(['summer.exe'], stdout=PIPE,stderr=PIPE,universal_newlines=True)
    print("ID of process running summer: {}".format(os.getpid()))
    stdout,stderror = subp1.communicate()
 #   return float(stdout)

if __name__ =='__main__':
        p1=multiprocessing.Process(target=balhbt, args=(10,10,par[0]))
        p2=multiprocessing.Process(target=balhbt, args=(10,10,par[1]))
        p3=multiprocessing.Process(target=balhbt, args=(10,10,par[2]))
        p4=multiprocessing.Process(target=balhbt, args=(10,10,par[3]))
        p5=multiprocessing.Process(target=balhbt, args=(10,10,par[4]))
        p6=multiprocessing.Process(target=balhbt, args=(10,10,par[5]))
        p7=multiprocessing.Process(target=balhbt, args=(10,10,par[6]))
        p8=multiprocessing.Process(target=balhbt, args=(10,10,par[7]))
        p9=multiprocessing.Process(target=balhbt, args=(10,10,par[8]))
        p10=multiprocessing.Process(target=balhbt, args=(10,10,par[9]))
        p11=multiprocessing.Process(target=balhbt, args=(10,10,par[10]))
        p12=multiprocessing.Process(target=balhbt, args=(10,10,par[11]))
        p13=multiprocessing.Process(target=summer)
        
        p1.start()
        p2.start()
        p3.start()
        p4.start()
        p5.start()
        p6.start()
        p7.start()
        p8.start()
        p9.start()
        p10.start()
        p11.start()
        p12.start()
        
        p1.join()
        p2.join()
        p3.join()
        p4.join()
        p5.join()
        p6.join()
        p7.join()
        p8.join()
        p9.join()
        p10.join()
        p11.join()
        p12.join()
        
#    for i in range(12):    
#        p[i].start()
#        print("ID of process running balhbt: {}".format(os.getpid()))
    
#    for i in range(12):
#        p[i].join()

        p13.start()

        p13.join()
        