import os
import time
from threading import *
def VCFpara(input,temp,reference,output,regions,t,Q,minD,maxD,keep=True,extension="*.fixed.mq.bam"):
    filter=''
    if maxD is not None: filter="\'QUAL>="+str(Q)+" && DP>"+str(minD)+" && DP<"+str(maxD)+"\'"
    else: filter="\'QUAL>="+str(Q)+" && DP>"+str(minD)+"\'"
    if regions is not None:
        os.system("cat "+regions+" | parallel -j "+str(t)+" \"bcftools mpileup -Ou --fasta-ref "+reference+" -r {} --max-depth 10000 "+input+extension+" | bcftools call -mv -Ov -V indels -o "+temp+"{}.vcf\"")
        os.system("bcftools concat "+temp+"*.vcf -Ov | bcftools filter -i "+filter+" -o "+output+"filtered.vcf")
        if keep is False: os.system('rm '+temp+"*.vcf")
    else:
        os.system("cat "+reference+".fai | cut -f 1 | parallel -j "+str(t)+" \"bcftools mpileup -Ou --fasta-ref "+reference+" -r {} --max-depth 10000 "+input+extension+" | bcftools call -mv -Ov -V indels -o "+temp+"{}.vcf\"")
        os.system("bcftools concat "+temp+"*.vcf -Ov | bcftools filter -i "+filter+" -o "+output+"filtered.vcf")
        if keep is False: os.system('rm '+temp+"*.vcf")     
def Main():
    starttime=time.time()
    # Directory path
    directory = '/home/ctp21002/data/FASTQFiles/Preprocess/bam/SNVnoBQSR/'
    #Output folder path
    output = "/home/ctp21002/data/FASTQFiles/vcf/"
    #hg38 reference path
    reference ="/home/ctp21002/GRCh38/hg38.fa"
    #Folder of the temp file
    temp= "/home/ctp21002/data/FASTQFiles/vcf/unmergedNoBQSR/"
    #Directory for the region txt files
    regions="/home/ctp21002/data/FASTQFiles/vcf/regions.txt"
    VCFpara(directory,temp,reference,output,regions, 192)
    print("Finished preparations for all files in "+str(time.time()-starttime) +" seconds.")
if __name__=="__main__":
    Main()
