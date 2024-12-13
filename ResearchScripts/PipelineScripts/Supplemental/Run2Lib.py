import os
import sys
import argparse

parser = argparse.ArgumentParser("Run2Lib")
parser.add_argument("--input", help="path to file containing the strings to be replaced", type=str)
parser.add_argument("--map", help="path to text file with each line consisting of a:b, where a will be replaced with b in the input", type=str)
parser.add_argument("--type", help="\'file\' if you want to replace the text within a file or \'directory\' if you want to replace the names of files in a directory", type=str)
args = parser.parse_args()
if args.type=='file':
    y=list()
    with open(args.map,"r") as f:
        for x in f:
            y.append(x.split(':'))
    for x in y:
        os.system('sed -i \'s/'+x[0].strip()+'/'+x[1].strip()+'/1\' '+args.input)
elif args.type=='directory':
    y=list()
    with open(args.map,"r") as f:
        for x in f:
            y.append(x.split(':'))
    for x in y:
        os.system('ls '+args.input+x[0].strip()+'* | cut -c '+str(len(args.input)+len(x[0].strip())+1)+'- | parallel \'mv '+args.input+x[0].strip()+'{} '+args.input+x[1].strip()+'{}\'')