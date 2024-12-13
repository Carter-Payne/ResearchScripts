import os
import time
from Preprocess import *
from SNV import *
from S2M2 import *
from CNA import *
import argparse
def ParseArgs():
    parser = argparse.ArgumentParser(prog='Pipeline', description="Accession List to SNV and CNA analysis files pipeline")
    parser.add_argument("-a", "--ACCESSION", metavar='Accession List', required=False,help="Path to the accession list of the sample(s) you want to process, if you already have the samples, specify their location using -ao", type=str)
    parser.add_argument("-ao","--ACCESSIONOUTPUT", metavar='FASTQ folder', required=True,help="Path to the output folder for the FASTQ files", type=str)
    parser.add_argument("-pd","--PREFETCHDIR", metavar='Prefetch folder', required=False, help="Path to the output folder of the prefetch directory, set by the SRA toolkit", type=str)
    parser.add_argument("-r", "--REFERENCE", metavar="Reference file", required=True, help="Path to the reference file",type=str)
    parser.add_argument("-i","--INDEX", help="Set if you need to index the reference file", default=False,action='store_true')
    parser.add_argument("-t", "--THREADS", metavar="Threads",help="Number of threads to use (default: Auto)",type=int,default=int(os.popen('nproc --all').read().strip()))
    parser.add_argument("-rg","--READGROUP", metavar="ReadGroup data", nargs='+', required=True, help = "Read group data for each run in the format of PL LB SAMN. If only the PL and LB are given, each run will instead be assigned as its own sample",type=str)
    parser.add_argument("-pu","--PLATFORMUNIT", help="select if you to add the specific platform unit found in the FASTQ file to each corresponding run. Only recommended to select if each run has a specific barcode, and will take much longer if selected. If not selected, each platform unit will be the Run ID", default=False, action= 'store_true')
    parser.add_argument("-temp", "--TEMP", metavar='Temporary folder', required=True,help="path to temporary folder of the outputs of the pipeline", type=str)
    parser.add_argument("-rem","--REMOVE", help="select if you wish to remove duplicates instead of marking them", default=False, action= 'store_true')
    parser.add_argument("-o","--OUTPUT",metavar='Bam output', required=True, help="Path to the finalized .bam files",type=str)
    parser.add_argument("-snv","--SNVPATH", metavar='SNV folder', help="Path to the output folder for the finished .vcf file. If not given, .vcf file won't be created", type=str)
    parser.add_argument("-cna","--CNAPATH", metavar='CNA folder', help="Path to the output folder for the .csv file. If not given, .csv file won't be created", type=str)
    parser.add_argument("-hg", "--HG", metavar="Genome type", choices=['hg38','hg19'], help="The type of reference genome being used, either hg38 or hg19(default: hg38)",type=str, default='hg38')
    parser.add_argument("-bqsr","--BASEQUALITY", metavar="Base Quality Recalibration", nargs='+', help = "Paths to each file to be used for BQSR. If none are given, BQSR will not be done.",type=str)
    parser.add_argument("-regions", "--REGIONS", metavar='VCFRegions', help="Path to a text file containing which chromosomes to perform variant calling on(default: all chromosomes in reference file)",type=str)
    parser.add_argument("-mq", "--MQ", metavar="Mapping Quality",help="Minimum mapping quality for read filtering(default: 40)",type=int,default=40)
    parser.add_argument("-pq", "--PQ", metavar="Phred Quality",help="Minimum phred quality for variant filtering (default: 40)",type=int,default=40)
    parser.add_argument("-mind", "--MINDEPTH", metavar="Minimum Read Depth",help="Minimum read depth for variant filtering(default: 10)",type=int,default=10)
    parser.add_argument("-maxd", "--MAXDEPTH", metavar="Maximum Read Depth",help="Maximum read depth for variant filtering(default: None)",type=int)
    parser.add_argument("-p","--PREFIX",metavar="File Prefix", required=False, help="Prefix for the finished CNA and SNV files.",type=str)
    parser.add_argument("-b","--BIN",metavar="Bin sizes", help="Bin sizes to be used for CNA analysis(default: 500kb)",type=int,default=500000)
    parser.add_argument("-d","--DEBUG", help="Select this option if you need to keep each intermitten file in the pipeline(only select this option if you have enough space to do so.)",default=False,action='store_true')
    parser.add_argument("-l","--LOG",metavar="Log file", help="Path for the file storing pipeline running times.(default: current directory)",type=str,default=os.getcwd())
    return parser.parse_args()

def ValidateInputs(args):
    '''Just a bunch of if statements to check if each folder and file exists. Also one to check if the number of threads <= number of cores and if the regions file is valid. Also converts any folder path that doesn't contain a / at the end to have one.'''
    check=True
    if args.ACCESSIONOUTPUT[-1]!='/':
        args.ACCESSIONOUTPUT+="/"
    if args.PREFETCHDIR is not None:
        if args.PREFETCHDIR[-1]!='/':
            args.PREFETCHDIR+="/"
    if args.TEMP[-1]!='/':
        args.TEMP+="/"
    if args.LOG[-1]!='/':
        args.LOG+="/"
    if args.OUTPUT[-1]!='/':
        args.OUTPUT+="/"
    if args.ACCESSION is not None:
        if os.path.exists(args.ACCESSION) is False:
            check=False
            print(args.ACCESSION +" does not exist")
    if os.path.exists(args.OUTPUT) is False:
        check=False
        print(args.OUTPUT +" does not exist")
    if os.path.exists(args.ACCESSIONOUTPUT) is False:
        check=False
        print(args.ACCESSIONOUTPUT +" does not exist")
    if args.PREFETCHDIR is not None:
        if os.path.exists(args.PREFETCHDIR) is False:
            check=False
            print(args.PREFETCHDIR +" does not exist")      
    if os.path.exists(args.REFERENCE) is False:
        check=False
        print(args.REFERENCE +" does not exist")     
    if os.path.exists(args.TEMP) is False:
        check=False
        print(args.TEMP +" does not exist")
    if args.SNVPATH is not None:
        if args.SNVPATH[-1]!='/':
            args.SNVPATH+="/"
        if os.path.exists(args.SNVPATH) is False:
            check=False
            print(args.SNVPATH +" does not exist")
    if args.CNAPATH is not None:
        if args.CNAPATH[-1]!='/':
            args.CNAPATH+="/"
        if os.path.exists(args.CNAPATH) is False:
            check=False
            print(args.CNAPATH +" does not exist")
        if args.PREFIX is None:
            check=False
            print("No Prefix Given")
    if len(args.READGROUP)!=2 and len(args.READGROUP)!=3:
        check=False
        print("Incorrect number of Read Groups arguments")     
    else:
        if len(args.READGROUP)==2:
            args.READGROUP.append(None)
    
    if args.REGIONS is not None:#might as well read the regions file to make sure it's valid(each region in the reference)
        checkRegions=True
        if os.path.isfile(args.REGIONS) is False:
            check=False
            print("Specified regions file does not exist")    
        else:
            if args.INDEX is True:
                print('Unable to validate the regions file, as the reference index has not been built. If you want to confirm its validity, index the reference file separately.')
            else:
                file=list()
                with open(args.REFERENCE+".fai",'r') as reference:
                    for f in reference:
                        file.append(f.split('\t')[0])
                with open(args.REGIONS,'r') as regions:
                    for x in regions:
                        if checkRegions==True:
                            if x.strip() not in file:
                                checkRegions=False
                                check=False
                                print("Not a valid regions file.")
                        else: pass
    if args.THREADS is not None:
        if int(os.popen('nproc --all').read().strip()) < args.THREADS:
            check=False
            print("You selected more threads than your available cores.")

    if args.BASEQUALITY is not None:
        for x in range(len(args.BASEQUALITY)):
            if os.path.isfile(args.BASEQUALITY[x]) is False:
                check=False
                print(args.BASEQUALITY[x] +" does not exist")

    if os.path.exists(args.LOG) is False:
        check=False
        print(args.LOG+" does not exist.")

    return check

def Pipeline(args):
    with open(args.LOG+"log.txt",'w+') as f:
        starttime=time.time()
        #FASTQ retrieval
        if args.ACCESSION is not None:
            print("Retrieving FASTQ files...")
            curr=time.time()        
            os.system("./FASTERQ.sh "+args.ACCESSION+" "+args.ACCESSIONOUTPUT+" "+args.PREFETCHDIR)
            run=str((time.time()-curr)/3600)
            print("Finished retreiving FASTQ files in "+run+" hours.")       
            f.write("Retrieving FASTQ files: "+run+" hours\n")     

        #BWA index 
        if args.INDEX:
            print("Indexing reference file...")
            curr=time.time()
            os.system("bwa index "+args.REFERENCE)                            
            run=str((time.time()-curr)/3600)
            print("Finished indexing reference in "+run+" hours.")       
            f.write("Indexing reference file: "+run+" hours\n")   
        
        #indexing FASTQ files
        print("Indexing runs...")
        curr=time.time()
        PU=bwa(args.ACCESSIONOUTPUT,args.REFERENCE,args.TEMP,args.THREADS,args.DEBUG)
        run=str((time.time()-curr)/3600)
        print("Finished indexing FASTQ files in "+run+" hours")
        f.write("Indexing FASTQ files: "+run+" hours\n")         
        
        #Samtools view > AddRG

        print("Adding Read Groups...")
        curr=time.time()
        if(args.PLATFORMUNIT):
            SamToBamAddRG(args.TEMP,args.TEMP,args.THREADS, args.READGROUP[0],args.READGROUP[1],args.READGROUP[2],PU,args.DEBUG)
        else: SamToBamAddRG(args.TEMP,args.TEMP,args.THREADS, args.READGROUP[0],args.READGROUP[1],args.READGROUP[2],None,args.DEBUG)
        run=str((time.time()-curr)/3600)
        print("Finished adding Read Groups in "+run+" hours")
        f.write("Converting to BAM and adding Read Groups: "+run+" hours\n")        

        #MarkDuplicates
        print("Marking Duplicates...")
        curr=time.time()
        MD(args.TEMP,args.TEMP,args.THREADS, args.DEBUG,args.REMOVE)
        run=str((time.time()-curr)/3600)
        print("Finished Marking Duplicates in "+run+" hours")
        f.write("Marking Duplicates and Sorting: "+run+" hours\n")     

        #FixTags
        print("Fixing tags...")
        curr=time.time()
        FixTags(args.TEMP,args.TEMP,args.REFERENCE,args.THREADS,args.DEBUG)
        run=str((time.time()-curr)/3600)
        print("Finished fixing tags in "+run+" hours")
        f.write("Fixing tags: "+run+" hours\n")

        #MQ
        print("Filtering mapping quality...")           
        curr=time.time()
        MapQ(args.TEMP,args.TEMP,args.THREADS,args.MQ,args.DEBUG)
        run=str((time.time()-curr)/3600)
        print("Finished filtering mapping quality in "+run+" hours")
        f.write("Filtering MQ: "+run+" hours\n")     
        
        #BQSR
        if args.BASEQUALITY is not None:
            SNP=''
            for x in args.BASEQUALITY:
                SNP+=" --known-sites "+x
            print("Applying BQSR...")
            curr=time.time()
            BQSRMult(args.TEMP,args.OUTPUT,args.REFERENCE,SNP,args.DEBUG)
            run=str((time.time()-curr)/3600)
            print("Finished BQSR in "+run+" hours")
            f.write("Fixing tags: "+run+" hours\n")      
        else:#If no BQSR, move the mapped bam files to the finished directory
            os.system("mv "+args.TEMP+"*.mq.bam "+args.OUTPUT)  
        
        #Indexing
        print("Indexing files...")
        curr=time.time()
        if args.BASEQUALITY is not None:
            Index(args.OUTPUT,args.THREADS,"*.BQSR.mq.bam")
        else:
            Index(args.OUTPUT,args.THREADS)
        run=str((time.time()-curr)/3600)
        print("Finished indexing files in "+run+" hours")
        f.write("Indexing: "+run+" hours\n")

        #SNV
        if args.SNVPATH is not None:
            print("Producing .vcf file...")
            curr=time.time()
            if args.BASEQUALITY is not None:
                VCFpara(args.OUTPUT,args.TEMP,args.REFERENCE,args.SNVPATH, args.REGIONS, args.THREADS, args.PQ, args.MINDEPTH, args.MAXDEPTH, args.DEBUG, "*.BQSR.mq.bam")
            else:
                VCFpara(args.OUTPUT,args.TEMP,args.REFERENCE,args.SNVPATH, args.REGIONS, args.THREADS, args.PQ, args.MINDEPTH, args.MAXDEPTH, args.DEBUG)
            run=str((time.time()-curr)/3600)
            print("Finished variant calling in "+run+" hours")
            f.write("Variant calling: "+run+" hours\n")

        #CNA(and creating MEDICC and DICE tsv files)
        if args.CNAPATH is not None:
            print("Producing .csv file")
            curr=time.time()
            if args.BASEQUALITY is not None:
                SeCNV(args.OUTPUT,args.CNAPATH,args.REFERENCE,args.HG,"*.BQSR.mq.bam",args.BIN,args.PREFIX)
            else:
                SeCNV(args.OUTPUT,args.CNAPATH,args.REFERENCE,args.HG,"*.fixed.mq.bam",args.BIN,args.PREFIX)
            run=str((time.time()-curr)/3600)
            CSVConvert(args.CNAPATH+args.PREFIX+".csv",args.CNAPATH+args.PREFIX+"DICE.tsv",mode="dice")
            CSVConvert(args.CNAPATH+args.PREFIX+".csv",args.CNAPATH+args.PREFIX+"MEDICC2.tsv",mode="medicc2")
            print("Finished creating CN matrix in "+run+" hours")
            f.write("Creating CN matrix: "+run+" hours\n")
        run=str((time.time()-starttime)/86400)
        f.write("Finished the entire pipeline in "+run+" days")
        f.close()
def Main():
    args=ParseArgs()
    print("Validating files...")
    check=ValidateInputs(args)
    if check==False:
        print("Validation of inputs failed.")
        exit(-1)
    print("Confirmed Validity of inputs. Starting pipeline.")
    Pipeline(args)
if __name__=="__main__":
    Main()


#TODO:
    #Implement Prefix in SNV and CNA
    #Create User manual(describing how it works, how all of the script files should be placed in the same directoried, which programs should be installed and set to path/environmental variables(and what to set them as))
    #Test that the validate all files works
    #Test that it works on current data set(Test up to SNV/CNA, then split and try to run SNV/CNA portions separately)
    #Note that SeCNV had to be slightly altered in how it called its python files to work with being called in a place outside of its own directory
    #Test that it works on a different run(Try to find one that can be successfully ran through both SNV and CNA methods.)
    #Move around files and folders in VM for better user experience and write a short manual denoting placement of everything(and PATH and Environmental variables)
