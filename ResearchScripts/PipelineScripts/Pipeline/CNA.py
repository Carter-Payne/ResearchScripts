import os
import time

def SeCNV(input,output,reference,genome,pattern,bin,prefix):
    os.system("python3 $SECNV "+input+" "+output+" "+reference+" -r "+genome+" --pattern "+pattern+" -b "+str(bin)+" && mv "+output+"cnv_matrix.csv "+output+prefix+".csv")
