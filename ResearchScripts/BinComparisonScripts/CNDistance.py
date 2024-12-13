import os
import sys
import argparse

def Parse():
    parser = argparse.ArgumentParser(prog="CNDistance",description="Returns the  distance of two different cell by location copy number matrices of equal bin length")
    parser.add_argument("-i1",'--INPUT1', help="path to the first DICE formatted tsv file", type=str)
    parser.add_argument("-i2",'--INPUT2', help="path to the second DICE formatted tsv file", type=str)
    args = parser.parse_args()
    return args

def FindDist(i1, i2):
    with open(i1,'r') as f:
        file1=f.readlines()
    with open(i2,'r') as f:
        file2=f.readlines()
    if len(file1)>=len(file2):#whichever file is larger is going to be the true tree created by CNASim
        Truth=sorted(file1[1:], key=lambda x: (int(x.split("\t")[0][4:]),int(x.split("\t")[1][3:]), int(x.split("\t")[3])))
        Infer=sorted(file2[1:], key=lambda x: (int(x.split("\t")[0][4:]),int(x.split("\t")[1][3:]), int(x.split("\t")[3])))
    else:
        Truth=sorted(file2[1:], key=lambda x: (int(x.split("\t")[0][4:]),int(x.split("\t")[1][3:]), int(x.split("\t")[3])))
        Infer=sorted(file1[1:], key=lambda x: (int(x.split("\t")[0][4:]),int(x.split("\t")[1][3:]), int(x.split("\t")[3]))) 
    j=0
    dist=0
    non=0
    i=0
    while(True):
            if(int(Truth[i].split('\t')[2])==int(Infer[j].split('\t')[2])):#if they are the same bin:
                #if round(float(Truth[i].split('\t')[4][:-1]))!=round(float(Infer[j].split('\t')[4][:-1])):dist+=1
                #if abs(float(Truth[i].split('\t')[4][:-1])-float(Infer[j].split('\t')[4][:-1])) >1 :dist+=1
                dist+=abs(float(Truth[i].split('\t')[4][:-1])-float(Infer[j].split('\t')[4][:-1]))
                j+=1
                i+=1
            else:
                if(int(Infer[j].split('\t')[3])-int(Infer[j].split('\t')[2]))!=999999:#edge case because the true profiles can sometimes go over 1000000 at the end of a chr, where in that case we just cut off the end of the inferred profiles so the starting regions still line up.
                    j+=1
                else: i+=1
                non+=1

            if i==len(Truth): break
    print(non)
    return dist
if __name__=="__main__":
    args=Parse()
    print(FindDist(args.INPUT1,args.INPUT2))

#TODO:
#Build on this code to compare two files now that they are in the sameish format
#FIRST: Find a way to sort each file, first by cell, then by chrom, then by position
#AFTER: iterate through now that they are a 1-1 comparison and get each 
