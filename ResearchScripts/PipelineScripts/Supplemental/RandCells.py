import random
import os
import argparse
import sys
import subprocess
def ParseArgs():
    parser = argparse.ArgumentParser(prog='RandCells', description="Program that reduces the size of a vcf file using random sampling of snv cells")
    parser.add_argument("-i", "--INPUT", metavar='Input', required=True,help="Path to the vcf file", type=str)
    parser.add_argument("-o", "--OUTPUT", metavar='Output', required=True,help="file path of output file in the form of /path/to/output.vcf", type=str)
    parser.add_argument("-n","--NUMBER", metavar="number",required=True, help="Number of SNV cells you want", type=int)
    return parser.parse_args()
def Validate(Args):
    check=True
    if Args.INPUT.endswith(".vcf") is False:
        print("incorect input file")
        check=False
    if Args.OUTPUT.endswith(".vcf") is False:
        print("incorect output file")
        check=False
    if check is False:
        exit(-1)
    samples = subprocess.check_output(['bcftools', 'query' ,'-l',Args.INPUT]).decode(sys.stdout.encoding).strip().split('\n')
    if (Args.NUMBER>len(samples)):
        print("random sampling number too large")
        exit(-1)
    return samples
def Main():
    Args=ParseArgs()
    samples=Validate(Args)
    randnums=list()
    keep=""
    i=Args.NUMBER
    while(i>0):
        x=random.randint(0,len(samples)-1)
        if x not in randnums:
            randnums.append(x)
            i-=1
    randnums.sort()
    with open('tmprs12321.txt','w') as f:
        for i in randnums:
            f.write(samples[i]+"\n")
    f.close()
    os.system('bcftools view -S tmprs12321.txt -o '+Args.OUTPUT+" "+ Args.INPUT)
    os.remove('tmprs12321.txt')
if __name__=="__main__":
    Main()
