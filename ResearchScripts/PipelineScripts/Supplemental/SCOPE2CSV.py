import os
import sys
import argparse

def Parse():
    parser = argparse.ArgumentParser(prog="SCOPE2CSV",description="converts the output of either the single or paired SCOPE Rscript into a .csv file in the same style as one generated by SeCNV.")
    parser.add_argument("-y",'--Y', help="path to Y.txt file created by the SCOPE Rscript", type=str)
    parser.add_argument("-CN",'--CN', help="path to CN.txt file created by the SCOPE Rscript", type=str)
    parser.add_argument("-o",'--OUTPUT', help="Output file in the form of path/name.csv", type=str)
    args = parser.parse_args()
    return args

def CSVConvert(Y, CN, output):
    csv=list()
    locations=list()
    csv.append("")
    i=0
    with open(Y,'r') as f:
        Yfile=f.readlines()
        f.close()
    with open(CN,'r') as f:
        CNfile=f.readlines()
        f.close()   
    #print(len(file))
    for i in range(1,len(Yfile)): #Yfile needed for chromosome regions
        csv[0]+=(','+Yfile[i].split(',')[0][1:-1])
    for i in range(len(CNfile)):
        if i==0:
            for j in range(1,len(CNfile[0].split(','))):
                csv.append(CNfile[0].split(',')[j].strip()[1:-1])
        else:
            for j in range(1,len(CNfile[0].split(','))):
                csv[j]+=(','+CNfile[i].split(',')[j].strip())
    with open(output, 'w') as f:
        for i in csv:
            f.write(i+'\n')
        f.close()

if __name__=="__main__":
    args=Parse()
    CSVConvert(args.Y,args.CN,args.OUTPUT)
