import os
import argparse
def Parse_Args():
    parser = argparse.ArgumentParser(prog='BinComparisons', description="Takes two folders and compares each tree inside them. Assumes trees inside the folders are of the format BinAlgorithm_TreeAlgorithm_kBP_tree.(nwk/new), and sorts by the kBP")
    parser.add_argument("-o","--OUTPUT",metavar='output', required=True, help="Path to the csv file to store all the results",type=str)
    parser.add_argument("-f1","--F1",metavar='Folder 1',required=True, help= 'The first folder of newick trees to compare',type=str)
    parser.add_argument("-f2","--F2",metavar='Folder 2',required=True, help= 'The second folder of newick trees to compare',type=str)    
    return parser.parse_args()
def Validate(args):
    check=True
    if args.F1[-1]!='/':
        args.F1+="/"
    if args.F2[-1]!='/':
        args.F2+="/"
    if os.path.exists(args.F1) is False:
        check=False
        print(args.F1 +" does not exist")
    if os.path.exists(args.F2) is False:
        check=False
        print(args.F2 +" does not exist")    
    return check
def sort(files):
    '''Given a list of files of X_X_Y, with Y being a number, sort the files in the list by Y.'''
    temp=list()
    for i in files: 
        try:
            temp.append(int(i.split("_")[2]))
        except IndexError:
            return files
    temp.sort()
    file=[None]*len(files)
    for i in files:file[temp.index(int(i.split("_")[2]))] =i
    return file
def compare(args):
    f1=os.listdir(args.F1)
    f2=os.listdir(args.F2)
    Trees1=list()
    Trees2=list()
    for i in f1:
        if i.endswith('.new') or i.endswith('.nwk'):
            Trees1.append(i)
    for i in f2:
        if i.endswith('.new') or i.endswith('.nwk'):
            Trees2.append(i)
    Trees1=sort(Trees1)
    Trees2=sort(Trees2)
    matrix=","
    for i in Trees1:
        matrix+=i[:-4]+","
    for i in Trees2:
        matrix+="\n"+i[:-4] +","
        for j in Trees1:
            x=os.popen('ete3 compare -r '+args.F1+j+' -t '+args.F2+i+' --unrooted | cut -d "|" -f 4').read()
            matrix+=x[21:26]+","
    with open(args.OUTPUT,'w') as f:
        f.write(matrix)
    f.close()

def Main():
    args=Parse_Args()
    if Validate(args) is False:
        exit(-1)
    else:
        compare(args)
if __name__ == "__main__":
    Main()
