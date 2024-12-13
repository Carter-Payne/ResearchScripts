import os
import sys
import argparse

def Parse():
    parser = argparse.ArgumentParser(prog="FormatTSV",description="Takes the profiles.tsv file produced by CNAsim and formats it in a way to be properly used by CNDistance.py to measure the distance")
    parser.add_argument("-i",'--INPUT', help="path to .tsv file created by CNASim", type=str)
    parser.add_argument("-o",'--OUTPUT', help="Output file in the form of path/name.tsv", type=str)
    args = parser.parse_args()
    return args

def TSVConvert(input, output):
    with open(input,'r') as f:
        file=f.readlines()
    for i in range(1,len(file)):
        temp=file[i].split('\t')
        states=temp[4][:-1].split(',')
        sum=float(states[0])+float(states[1])
        file[i]=temp[0]+"\t"+temp[1]+"\t"+str(int(temp[2])+1)+"\t"+temp[3]+"\t"+str(sum)+"\n"
    file.sort()
    file[-1]=file[-1][:-1]
    with open(output, 'w') as f:
        for i in file:
            f.write(i)
        f.close()
if __name__=="__main__":
    args=Parse()
    TSVConvert(args.INPUT,args.OUTPUT)
