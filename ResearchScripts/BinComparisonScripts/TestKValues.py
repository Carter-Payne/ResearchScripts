import os
import argparse
def Parse_Args():
    parser = argparse.ArgumentParser(prog='MultiK', description="Tests multiple K values for SCOPE")
    parser.add_argument("-o","--OUTPUT",metavar='Folder output', required=True, help="Path to the folder to store all of the heatmaps and trees",type=str)
    parser.add_argument("-k","--K",metavar='K values', nargs='+',required=True, help= 'All K values you wish to try',type=str)
    parser.add_argument('-pre','--PREFIX',metavar='Prefix',required=True,help='Prefix to be used for each matrix, which will take the form of (Prefix)(bin size).csv',type=str)
    parser.add_argument('-mode','--MODE',metavar='mode',required=True,help='Tree reconstruction algorithm you wish to use and have the conda env activate for, currently supporting medicc2 & dice',type=str,choices=['medicc2','dice'])
    return parser.parse_args()
def Validate(args):
    check=True
    if args.OUTPUT[-1]!='/':
        args.OUTPUT+="/"
    if os.path.exists(args.OUTPUT) is False:
        check=False
        print(args.OUTPUT +" does not exist")       
    if check: return
    else: exit(-1)
def Create_Matrices(args):
    for i in args.K:
        os.system('Rscript $SCOPESIM '+str(i))
        os.system('python3 $SCOPE2CSV -y K'+str(i)+'_Y.txt -CN K'+str(i)+'_CN.txt -o '+args.OUTPUT+args.PREFIX+'K'+str(i)+'.csv')
        os.system('python3 $S2M2 -i '+args.OUTPUT+args.PREFIX+'K'+str(i)+'.csv -o '+args.OUTPUT+args.PREFIX+'K'+str(i)+'.tsv -mode medicc2')
        os.system('medicc2 '+args.OUTPUT+args.PREFIX+'K'+str(i)+'.tsv  '+args.OUTPUT+' --total-copy-numbers -j '+str(int(os.popen('nproc --all').read().strip()))+' -a \'CN\' -p '+args.PREFIX+'K'+str(i))
def Main():
    args=Parse_Args()
    Validate(args)
    Create_Matrices(args)

if __name__=='__main__':
    Main()
