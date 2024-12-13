import random
import subprocess
import argparse
import sys
def ParseArgs():
    parser = argparse.ArgumentParser(prog='RandSample', description="Program that reduces the size of a vcf file using random sampling of snv loci")
    parser.add_argument("-i", "--INPUT", metavar='Input', required=True,help="Path to the vcf file", type=str)
    parser.add_argument("-o", "--OUTPUT", metavar='Output', required=True,help="file path of output file in the form of /path/to/output.vcf", type=str)
    parser.add_argument("-n","--NUMBER", metavar="number",required=True, help="Number of SNV loci you want", type=int)
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
    locations = subprocess.check_output(['bcftools', 'view' ,'-H',Args.INPUT]).decode(sys.stdout.encoding).strip().split('\n')
    if (Args.NUMBER>len(locations)):
        print("random sampling number too large")
        exit(-1)
    return locations
def Main():
    Args=ParseArgs()
    location=Validate(Args)
    randnums=list()
    i=Args.NUMBER
    while(i>0):
        x=random.randint(0,len(location)-1)
        if x not in randnums:
            randnums.append(x)
            i-=1
    randnums.sort()
    f=open(Args.OUTPUT,'w')
    f.write(subprocess.check_output(['bcftools', 'head' ,Args.INPUT]).decode(sys.stdout.encoding).strip()+'\n')
    for i in range(len(randnums)):
        f.write(location[randnums[i]])
        if i!=len(randnums): f.write('\n')
    f.close()
if __name__=="__main__":
    Main()