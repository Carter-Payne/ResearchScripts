import os
import time
from threading import *
def bwa(input,reference, output,t,keep=True):
    '''Method used to convert the fastq files into sam files'''
    files = os.listdir(input)
    i = 0
    run=""
    #Creation of the key:value pairs for the PU used in the read groups; have to do it here because FASTQ files are the only place that hold the information and they get deleted during the next step.
    mapping={}
    '''
    for file in files:
        if file.endswith("_2.fastq") is not True:
            key=file.split('.')[0].split('_')[0]
            data=os.popen('head -1 '+input+file).read().strip().split(" ")[1]
            if(":" in data):
                data=data.split(":")
                value=data[0]+data[1]+data[2]+"."+data[3]
            else: value="machine"
            mapping[key]=value
    '''
    while i < len(files):
        #curr=time.time()
        file = files[i]
        # Check if the file ends with _1.fastq
        if file.endswith('_1.fastq'):
            # Print both _1.fastq and the next file
            run=file[:-8]
            filecomb=input+file + " " + input+run+"_2.fastq"
            #print("starting run "+run)
            os.system("bwa mem -M -t "+str(t)+" "+reference+" " + filecomb +" > "+output+run+".sam -v 0")
            if keep is False:
                os.system("rm "+filecomb)
            #print(run +" successfully prepared in "+str(time.time()-curr) +" seconds.")
            i+=1   
        elif file.endswith('_2.fastq'):
            i+=1
            continue
        else:
            run=file[:-6]    
            file1=input+file
            #print("starting run "+run)
            os.system("bwa mem -M -t "+str(t)+" "+reference+" " + file1 +" > "+output+run+".sam -v 0")
            if keep is False:
                os.system("rm "+file1)
            #print(run +" successfully prepared in "+str(time.time()-curr) +" seconds.")   
            i+=1  
    return mapping

def SamToBamAddRG(input, output, t, PL,LB,SAMN=None,PU=None ,keep=True):
    '''Converts the sam files obtained with bwa into bamfiles and adds read groups to them. deletes sam files'''
    if PU is not None: #In the case that they want to have individual Platform Units(due to the interaction between linux and python I can't parallelize the platform units :( )
        files=os.listdir(input)
        for file in files:
            run=file.split(".")[0]
            if SAMN==None:SAMN=run
            if keep:os.system("samtools view -bS "+input+file+" | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+run+".rg.bam -RGID "+run+" -RGPL "+PL+" -RGSM "+SAMN+" -RGLB "+LB+" -RGPU "+PU[run]+" -QUIET TRUE > /dev/null 2>&1 ")
            else:os.system("samtools view -bS "+input+file+" | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+run+".rg.bam -RGID "+run+" -RGPL "+PL+" -RGSM "+SAMN+" -RGLB "+LB+" -RGPU "+PU[run]+" -QUIET TRUE > /dev/null 2>&1 && rm "+input+file)
    else:#PU=ID for all inputs
        if SAMN==None:
            if keep:os.system("ls "+input+"*.sam | parallel -j "+str(t)+" \"samtools view -bS {} | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+"{/.}.rg.bam -RGID {/.} -RGPL "+PL+" -RGSM {/.} -RGLB "+LB+" -RGPU {/.} -QUIET TRUE > /dev/null 2>&1 \"")#&& rm "+temp+"{}\"")
            else: os.system("ls "+input+"*.sam | parallel -j "+str(t)+" \"samtools view -bS {} | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+"{/.}.rg.bam -RGID {/.} -RGPL "+PL+" -RGSM {/.} -RGLB "+LB+" -RGPU {/.} -QUIET TRUE > /dev/null 2>&1 && rm "+input+"{}\"")
        else:
            if keep:
                os.system("ls "+input+"*.sam | parallel -j "+str(t)+" \"samtools view -bS {} | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+"{/.}.rg.bam -RGID {/.} -RGPL "+PL+" -RGSM "+SAMN+" -RGLB "+LB+" -RGPU {/.} -QUIET TRUE > /dev/null 2>&1 \"")#&& rm "+temp+"{}\"")
            else: os.system("ls "+input+"*.sam | parallel -j "+str(t)+" \"samtools view -bS {} | java -jar $PICARD AddOrReplaceReadGroups -I /dev/stdin -O "+output+"{/.}.rg.bam -RGID {/.} -RGPL "+PL+" -RGSM "+SAMN+" -RGLB "+LB+" -RGPU {/.} -QUIET TRUE > /dev/null 2>&1 && rm "+input+"{}\"")

def MD(input,output, t, keep=True,remove=False):
    '''MarksDuplicates and sorts by coordinates'''
    tempo=os.listdir(input)
    files=list()
    for f in tempo:
        if f.endswith(".rg.bam"):
            files.append(f)
    i = 0
    run= ""
    while i< len(files):
        file=files[i]
        run=file.split('.')[0]
        #curr=time.time()
        #print("starting run "+run)        
        if keep: os.system("gatk MarkDuplicatesSpark -I "+input+file+" -O "+output+run+".sorted.rg.dedup.bam --conf \'spark.executor.cores="+str(t)+"\' --remove-all-duplicates "+str(remove)+" --create-output-bam-splitting-index false --create-output-bam-index false --QUIET true > /dev/null 2>&1")# && rm "+temp+file)
        else:os.system("gatk MarkDuplicatesSpark -I "+input+file+" -O "+output+run+".sorted.rg.dedup.bam --conf \'spark.executor.cores="+str(t)+"\' --remove-all-duplicates "+str(remove)+" --create-output-bam-splitting-index false --create-output-bam-index false --QUIET true > /dev/null 2>&1 && rm "+input+file)
        #print(run +" successfully prepared in "+str(time.time()-curr) +" seconds.")        
        i+=1

def FixTags(input,output,reference, t, keep=True):
    '''Fixes the tags that can be erred by MarkDuplicates and removes the files that are from MarkDuplicates'''
    if keep: os.system("ls "+input+"*.dedup.bam | parallel -j "+str(t)+" \"java -jar $PICARD SetNmMdAndUqTags -R "+reference+" -I {} -O "+output+"{/.}.fixed.bam -QUIET TRUE > /dev/null 2>&1\"")# && rm {}\"")
    else: os.system("ls "+input+"*.dedup.bam | parallel -j "+str(t)+" \"java -jar $PICARD SetNmMdAndUqTags -R "+reference+" -I {} -O "+output+"{/.}.fixed.bam -QUIET TRUE > /dev/null 2>&1 && rm {}\"")

def MapQ(input, output, t, Q, keep=True):
    '''Removes reads from files in input with a mapping quality <Q and send to output'''
    if keep: os.system("ls "+input+"*.fixed.bam | parallel -j "+str(t)+" \"samtools view {} -o "+output+"{/.}.mq.bam -q "+str(Q)+"\"")
    else:os.system("ls "+input+"*.fixed.bam | parallel -j "+str(t)+" \"samtools view {} -o "+output+"{/.}.mq.bam -q "+str(Q)+" && rm {}\"")

def Merge(input,output,t):
    '''Merges the bam files into one'''
    os.system("samtools merge "+output+"merged.bam "+input+"*.bam --threads "+str(t)+" > /dev/null 2>&1")

def Index(input,t,extension='*.fixed.mq.bam'):
    '''indexes the file(s) specified'''
        #print("starting run "+run)
    os.system("ls "+input+extension+" | parallel -j "+str(t)+" \"samtools index {} \"")
        #print(run +" successfully prepared in "+str(time.time()-curr) +" seconds.")

def BQSRMult(input,output, reference, SNP,keep=True):
    '''Use if you want to run BQSR on the pre merged files'''
    temp = os.listdir(input)
    files=list()
    for file in temp:
        if file.endswith(".fixed.bam"):
            files.append(file)
    i = 0
    run=""
    while i < len(files) :
        file = files[i]
        run=file.split('.')[0]
        #curr=time.time()
        #print("starting run "+run)
        if keep:
            os.system("gatk --java-options \"-Xmx150G\" BaseRecalibrator  -I "+input+file+" -R "+reference+SNP+" -O "+output+run+".table --QUIET true > /dev/null 2>&1")
            os.system("gatk --java-options \"-Xmx150G\" ApplyBQSR  -I "+input+file+" -R "+reference+" --bqsr-recal-file "+output+run+".table -O "+output+run+".sorted.rg.dedup.fixed.BQSR.mq.bam --create-output-bam-index false --QUIET true > /dev/null 2>&1")   
        else:
            os.system("gatk --java-options \"-Xmx150G\" BaseRecalibrator  -I "+input+file+" -R "+reference+SNP+" -O "+output+run+".table --QUIET true > /dev/null 2>&1")
            os.system("gatk --java-options \"-Xmx150G\" ApplyBQSR  -I "+input+file+" -R "+reference+" --bqsr-recal-file "+output+run+".table -O "+output+run+".sorted.rg.dedup.fixed.BQSR.mq.bam --create-output-bam-index false --QUIET true > /dev/null 2>&1 && rm "+output+run+".table "+input+file)
        #print(run +" successfully prepared in "+str(time.time()-curr) +" seconds.")
        i+=1   
def BQSR(merged, reference, SNP,dbSNP,gold,known):
    '''Calculates and applies the BQSR to a single file'''
    os.system("gatk --java-options \"-Xmx150G\" BaseRecalibrator  -I "+merged+"merged.bam -R "+reference+" --known-sites "+SNP+gold+" --known-sites "+SNP+dbSNP+" --known-sites "+SNP+known+" -O "+merged+"merged.table --QUIET true > /dev/null 2>&1")
    os.system("gatk --java-options \"-Xmx150G\" ApplyBQSR  -I "+merged+"merged.bam -R "+reference+" --bqsr-recal-file "+merged+"merged.table -O "+merged+"mergedBQSR.bam --QUIET true > /dev/null 2>&1 && rm "+merged+"merged.bam")




def Main():
    starttime=time.time()
    #Number of desired bash threads
    t= str(192)
    #Name of sample, can later on use dict for multisample preprocessing
    SAMN="SAMN04893778"
    #Name of the library that created the sample
    LB="HiSeq-2500"
    #Name of the sample's platform
    Platform="ILLUMINA"
    # Directory path
    directory = '/home/ctp21002/data/FASTQ/'
    #Output folder path
    output = "/home/ctp21002/data/Preprocess/bam/"
    #Output folder path for BQSR
    BQSR = "/home/ctp21002/data/Preprocess/bam/BQSR/"
    #hg38 reference path
    reference ="/home/ctp21002/GRCh38/hg38.fa"
    #SNP Command for BQSR:
    SNP="--known-sites /home/ctp21002/GRCh38/SNP/hg38_dbSNP.vcf --known-sites /home/ctp21002/GRCh38/SNP/hg38_gold.vcf.gz --known-sites /home/ctp21002/GRCh38/SNP/hg38_indels.vcf.gz"
    #Folder of the merged bam file and BQSR table
    merged="/home/ctp21002/data/Preprocess/mergedbam/"
    #Metrics folder path
    metrics="/home/ctp21002/data/Preprocess/Metrics/"
    #temporary folder for holding items mid processing (All other temp files will be deleted when the pipeline is cleaned up)
    temp= "/home/ctp21002/data/Preprocess/temp/"
    #temp file for storing Dedup files
    tempDedup="/home/ctp21002/data/Preprocess/tempDedup/"
    #temp file for storing Fixed files
    tempFixed="/home/ctp21002/data/Preprocess/tempFixed/"
    #File names of the known indels, gold standard indels, and the dbsnp used for BQSR
    dbSNP="hg38_dbSNP.vcf"
    gold_standard="hg38_gold.vcf.gz"
    known_indels="hg38_indels.vcf.gz"
    print(len(bwa(directory, reference, temp,t)))
    #print("Finished indexing for all files in "+str((time.time()-starttime)/3600) +" hours.")
    '''
    newtime=time.time()
    SamToBamAddRG(temp,temp,SAMN,Platform,LB)
    print("Finished converting to bam and adding read groups for all files in "+str((time.time()-newtime)/3600) +" hours.")
    
    newtime=time.time()
    MD(temp,tempDedup)
    print("Finished marking duplicates and sorting coordinates for all files in "+str((time.time()-newtime)/3600) +" hours.")

    newtime=time.time()
    FixTags(tempDedup,tempFixed,reference)
    print("Finished fixing tags for all files in "+str((time.time()-newtime)/3600) +" hours.")

    newtime=time.time()
    MapQ(tempFixed,output,40)
    print("Finished adjusting mapping quality for all files in "+str((time.time()-newtime)/3600) +" hours.")
    
    newtime=time.time()
    Index(output,t)
    print("Finished indexing the file(s) in "+str((time.time()-newtime)/3600) +" hours.")
    
    newtime=time.time()
    BQSRMult(output,BQSR,reference,SNP)
    print("Finished recalibrating the file(s) in "+str((time.time()-newtime)/3600) +" hours.")

    newtime=time.time()
    Index(output+"BQSR/",t)
    print("Finished indexing the file(s) in "+str((time.time()-newtime)/3600) +" hours.")

    print("Finished preparations for all files in "+str(time.time()-starttime) +" seconds.")'''
if __name__=="__main__":
    Main()
