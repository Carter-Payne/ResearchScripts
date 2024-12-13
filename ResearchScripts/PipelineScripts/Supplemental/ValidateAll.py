import os
import sys
import getopt
import time
import subprocess
from threading import *

def usage():
    print("usage: python3 ValidateAll.py /path/to/bam/files")
starttime=time.time()
limit=Semaphore(64)
def ValidateSam(limit,file,run):
    curr=time.time()
    try:
        check=(subprocess.check_output(["java","-jar","-Xmx40G","$PICARD","ValidateSamFile","-I",file,"-MODE","SUMMARY"],stderr=subprocess.DEVNULL))
        check=(check.decode("ascii"))
        if(check=="No errors found\n"):
            print(run + " validated in "+str(time.time()-curr)+" seconds")
        else:
            print("error in run "+run)
    finally:
        limit.release()

# Directory path
def Main():
    if(len(sys.argv)!=2):
        usage()
        exit(-1)
    directory = sys.argv[1]
    if (directory[-1]!="/"):
        directory+="/"
    Ts=[]
    files = os.listdir(directory)
    i=0
    run=''
    while i< 70:
        limit.acquire()
        file=directory+files[i]
        run=files[i][0:10]
        T=Thread(target=ValidateSam,args=(limit,file,run))
        T.start()
        Ts.append(T)
        i+=1

    for T in Ts:
        T.join()

if __name__=="__main__":
    Main()
