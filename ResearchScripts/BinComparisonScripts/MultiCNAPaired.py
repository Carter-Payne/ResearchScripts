import os
import time
import argparse
import numpy as np
import pandas as pd
import pysam
from scipy.signal import find_peaks

def Parse_Args():
    parser = argparse.ArgumentParser(prog='MultiCNAPaired', description="Creates multiple CNA matrices of varying bin sizes for a given set of bam files that are paired end")
    parser.add_argument("-csvo","--CSVOUTPUT",metavar='CSV output', required=True, help="Path to the folder to store all of the heatmaps",type=str)
    parser.add_argument("-tsvo","--TSVOUTPUT",metavar='TSV output', required=True, help="Path to the folder to store all of the converted heatmaps into the analysis ready tsv files",type=str)
    parser.add_argument("-treeo","--TREEOUTPUT",metavar='Tree output', required=True, help="Path to the folder to store all of the converted trees",type=str)
    parser.add_argument("-to","--TEMPOUTPUT",metavar='Temporary Folder output', required=True, help="Path to the folder to store the intermediary files created by SECNV or SCOPE; everything in this folder will be deleted.",type=str)
    parser.add_argument("-i","--INPUT",metavar='Input', required=True, help="Path to the folder containing the .bam files",type=str)
    parser.add_argument("-l","--LIST",metavar='Bin size list',required=True, help= 'Path to the text file containing the values of the bin sizes, in bp, each separated by a new line',type=str)
    parser.add_argument("-r","--REFERENCE",metavar='Reference',required=True,type=str,help='path to hg19/38 reference file')
    parser.add_argument("-p","--PATTERN",metavar='.bam file pattern',help='pattern of each bam file (default:*.bam)',type=str,default="*.bam")
    parser.add_argument("-hg", "--HG", metavar="Genome type", choices=['hg38','hg19'], help="The type of reference genome being used, either hg38 or hg19(default: hg38)",type=str, default='hg38')
    parser.add_argument('-pre','--PREFIX',metavar='Prefix',required=True,help='Prefix to be used for each matrix, which will take the form of (Prefix)(bin size).csv',type=str)
    parser.add_argument('-tree','--TREE',metavar='tree reconstruction algorithm',required=True,help='Tree reconstruction algorithm you wish to use and have the conda env activate for, currently supporting medicc2 & dice',type=str,choices=['medicc2','dice'])
    parser.add_argument('-mode','--MODE',metavar='Copy Number reconstruction',required=True,help='which copy number method you wish to use, either SECNV or SCOPE',type=str,choices=['SECNV','SCOPE'])
    return parser.parse_args()
def Validate(args):
    check=True
    if args.TSVOUTPUT[-1]!='/':
        args.TSVOUTPUT+="/"
    if args.CSVOUTPUT[-1]!='/':
        args.CSVOUTPUT+="/"
    if args.TREEOUTPUT[-1]!='/':
        args.TREEOUTPUT+="/"
    if args.INPUT[-1]!='/':
        args.INPUT+="/"
    if args.TEMPOUTPUT[-1]!='/':
        args.TEMPOUTPUT+="/"
    if os.path.exists(args.INPUT) is False:
        check=False
        print(args.INPUT +" does not exist")
    if os.path.exists(args.CSVOUTPUT) is False:
        check=False
        print(args.CSVOUTPUT +" does not exist")        
    if os.path.exists(args.TEMPOUTPUT) is False:
        check=False
        print(args.TEMPOUTPUT +" does not exist")    
    if os.path.exists(args.TSVOUTPUT) is False:
        check=False
        print(args.TSVOUTPUT +" does not exist")        
    if os.path.exists(args.TREEOUTPUT) is False:
        check=False
        print(args.TREEOUTPUT +" does not exist")    
    bins=list()
    if args.LIST.endswith('.txt') is False:
        check=False
        print("Incorrect bin size format")
    else:
        with open(args.LIST,'r') as f:
            for i in f.readlines():
                bins.append(i.strip())
            f.close()
    if check: return bins
    else: exit(-1)
def Create_Matrices(input,seoutput,csvoutput,tsvoutput,treeoutput,reference,genome,pattern,bins,prefix,tree,mode):
    if tree == 'dice':
        for i in bins:
            try:
                if os.path.isfile(csvoutput+str(i)+'.csv') is False:
                    if mode=='SECNV':
                        SeCNV(input,seoutput,reference,genome,pattern,i)
                        os.system('mv '+seoutput+'cnv_matrix.csv '+csvoutput+str(i)+'.csv')
                        os.system('rm '+seoutput+"*")
                    else:
                        os.system('Rscript $GINIPaired '+str(i)[:-3]+' '+seoutput)
                        Gini=GetGini(seoutput+'gini.csv')
                        os.system('rm '+seoutput+'gini.csv && Rscript $SCOPEPaired '+str(i)[:-3]+' '+str(Gini))
                        os.system('python3 $SCOPE2CSV -y '+str(i)[:-3]+'_Y.txt -CN '+str(i)[:-3]+'_CN.txt -o '+csvoutput+str(i)+'.csv && rm '+str(i)[:-3]+'_Y.txt '+str(i)[:-3]+'_CN.txt')

                os.system('python3 $S2M2.py -i '+csvoutput+str(i)+'.csv -o '+tsvoutput+prefix+str(i)[:-3]+'.tsv -mode '+tree)
                os.system('dice -i '+tsvoutput+prefix+str(i)[:-3]+'.tsv -o '+treeoutput+' -p '+prefix+str(i)[:-3]+' -m balME -t')
            except KeyError:
                pass
    else:
        for i in bins:
            try:
                if os.path.isfile(csvoutput+str(i)+'.csv') is False:
                    if mode=='SECNV':
                        SeCNV(input,seoutput,reference,genome,pattern,i)
                        os.system('mv '+seoutput+'cnv_matrix.csv '+csvoutput+str(i)+'.csv')
                        os.system('rm '+seoutput+"*")
                    else:
                        os.system('Rscript $GINIPaired '+str(i)[:-3]+' '+seoutput)
                        Gini=GetGini(seoutput+'gini.csv')
                        os.system('rm '+seoutput+'gini.csv && Rscript $SCOPEPaired '+str(i)[:-3]+' '+str(Gini))
                        os.system('python3 $SCOPE2CSV -y '+str(i)[:-3]+'_Y.txt -CN '+str(i)[:-3]+'_CN.txt -o '+csvoutput+str(i)+'.csv && rm '+str(i)[:-3]+'_Y.txt '+str(i)[:-3]+'_CN.txt')
                os.system('python3 $S2M2.py -i '+csvoutput+str(i)+'.csv -o '+tsvoutput+prefix+str(i)[:-3]+'.tsv -mode '+tree)
                os.system('medicc2 '+tsvoutput+prefix+str(i)[:-3]+'.tsv  '+treeoutput+' --total-copy-numbers -j '+str(int(os.popen('nproc --all').read().strip()))+' -a \'CN\' -p '+prefix+str(i)[:-3])
            except KeyError:
                pass
def SeCNV(input,output,reference,genome,pattern,bin):
    os.system("python3 $SECNV "+input+" "+output+" "+reference+" -r "+genome+" --pattern "+pattern+" -b "+str(bin))
def get_normal_from_gini_auto(gini, num_boxes = None, min_count = 10, max_gini = 0.3):#Obtained from Samson, slightly modified
    gini_dict = gini.to_dict()
    gini_flattened = list(gini_dict.items())
    gini_flattened.sort(key = lambda x: x[1])
    [cell_names, g_vals] = map(list, zip(*gini_flattened))


    minVal, maxVal = min(g_vals), max(g_vals)
    num_cells = len(cell_names)

    if num_boxes == None:
        num_boxes = 100
    interval_len = 1 / num_boxes
    intervals = [interval_len*i for i in range(1, num_boxes+1)]
    boxes = [[] for i in range(num_boxes)]

    cur_idx, cur_box = 0, 0
    while cur_idx < num_cells:
        while cur_box < len(intervals) and g_vals[cur_idx] > intervals[cur_box]:
            cur_box += 1
        if cur_box == len(intervals):
            boxes[cur_box].extend([i for i in range(cur_idx, num_cells)])
            break
        else:
            boxes[cur_box].append(cur_idx)
        cur_idx += 1
    
    box_counts = [len(x) for x in boxes]
    min_height = max(round(num_cells*0.05), 1)
    peaks = find_peaks(box_counts, height=min_height)
    if len(peaks[0]) < 2:
         peaks = find_peaks(box_counts)
    if len(peaks[0]) == 0:
        print('Error')
        return

    if len(peaks[0]) == 1:
        end_idx = peaks[0][0] + 1
        num_normal = sum(box_counts[:end_idx])
        while num_normal < min_count:
            end_idx += 1
            num_normal = sum(box_counts[:end_idx])
    else:
        i = 1
        end_idx = int(np.floor((peaks[0][i-1] + peaks[0][i])/2))
        num_normal = sum(box_counts[:end_idx])
        while i+1 < len(peaks[0]) and num_normal <= min_count: 
            i += 1
            end_idx = int(np.floor((peaks[0][i-1] + peaks[0][i])/2))
            num_normal = sum(box_counts[:end_idx])
    normal_cells = []
    for i in range(end_idx):
        for j in boxes[i]:
            if g_vals[j] <= max_gini:
                normal_cells.append(cell_names[j])

    return normal_cells
def GetGini(gini):
    with open(gini,'r') as f:
        file=f.readlines()
        f.close()
    file.pop(0)
    inde=list()
    for i in range(len(file)):
        file[i]=file[i][:-1]
    inde.append('Gini')
    df=pd.Series({(i+1):float(file[i]) for i in range(len(file))})
    #print(df)
    normal=(get_normal_from_gini_auto(df))
    max=-1
    for i in normal:
        if(float(file[i-1])>max):max=float(file[i-1])
    return(max)
def Main():
    args=Parse_Args()
    bin_sizes=Validate(args)
    Create_Matrices(args.INPUT, args.TEMPOUTPUT, args.CSVOUTPUT,args.TSVOUTPUT,args.TREEOUTPUT, args.REFERENCE, args.HG, args.PATTERN, bin_sizes, args.PREFIX,args.TREE,args.MODE )
if __name__=="__main__":
    Main()

