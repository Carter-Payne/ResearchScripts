import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

def skim(n,len):
    fin=[]
    fin.append(0)
    i=1
    while i<len:
        fin.append(i)
        i+=n
    return fin

def rename(data):
    for x in data.columns:
        if(x[4]==':'):
            data.rename({x:x[0:4]},axis=1,inplace=True)
        else:
            data.rename({x:x[0:5]},axis=1,inplace=True)
    return data
def getPlotLabels(data):
    #print(len(data.columns))
    fin=[]
    prev=''
    curr=''
    for x in data.columns:
        if(x[4]==':'):
            curr=x[0:4]
        else:
            curr=x[0:5]
        if(prev==curr):
            fin.append('')
        elif(prev==''):
            fin.append(curr+" start")
            prev=curr
        elif(prev!=curr):
            #if (prevx not in fin):
                #fin.pop()
                #fin.append(prev+" end")
            fin.append(curr+" start")
            prev=curr
    #print(len(fin))
    return fin
def extract_csv_gen_plot(csv_path):
    
    data = pd.read_csv(csv_path,index_col=0)# usecols=skim(50,5075))
    data = data.drop(data.columns[[0, 1]], axis=1)
    cmap_dict = {0: '#152d6e', 1: '#5cbcd1', 2: '#FFFFFF', 3: '#e0de55', 4: '#e68932', 5: '#e64a32', 6: '#541910', 7: '#541910', 8: '#541910', 9: '#541910', 10: '#541910'}
    cmap = ListedColormap([cmap_dict[i] for i in range(10)])
    #data=rename(data)
    #data.index.names = ['Name']
    g = sns.heatmap(data,xticklabels=getPlotLabels(data),cmap=cmap,vmin=-0.5, vmax=9.5, yticklabels=True)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_title('1500kb bins')
    #plt.tight_layout()
    plt.show()
    '''data = data.drop(data.columns[[0, 1]], axis=1)
    data.index.names = ['Name']
    g = sns.heatmap(data)
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_title('Heatmap')
    plt.tight_layout()
    plt.show()'''


extract_csv_gen_plot("/Users/carterpayne/Desktop/MatrixBQSRDedup1500kb.csv")