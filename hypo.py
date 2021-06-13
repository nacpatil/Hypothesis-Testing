
from dict import *
import numpy as np
import os
import pandas as pd
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import matplotlib.pyplot as splt
import seaborn as sb
import string
import pandas as pd
  
#=====INPUTS===========================================================================================
peak_file="peaks/sample_peak.bed"
sample_file="binned_data/sample_.bedgraph"
noise_file="binned_data/noise_.bedgraph"
out_file = "pvalues_peak_file.bed"
#======READING=============================================================================================




class bedgraph:
    def __init__(self,filename ):
        self.filename=filename
        data=read_pdata(filename, 4, True)
        [self.chrs, self.st, self.en, self.enr]= data
        [self.chr_list, self.start_ind, self.end_ind]=chr_starts(self.chrs)
        
        self.bin_width = data[2][0]-data[1][0]
        print(self.bin_width)

        
class peak_class:
    def __init__(self,filename ):
        data=read_pdata(filename, 5, False)
        [self.chrs, self.st, self.en, self.width, self.log]= data
        self.num_peaks=len(self.chrs)
        
    def write_pvalues(self,out_file,stat,pvalue):
        f=open(out_file,'w')
        f.write("#chr start end width log ks-statistic pvalue\n")
        for i in range(0,self.num_peaks):
            out_str =   str(self.chrs[i]) + " " +str(self.st[i])+ " "+str(self.en[i])+ " " + \
            str(self.width[i])+" "+str(self.log[i])+" "+str(stat[i]) + " "+str(pvalue[i]) +"\n"  
            f.write(out_str)
        f.close()
        
    
chrdict = {
	"chr1": "1",
	"chr2": "2",
	"chr3": "3",
	"chr4": "4",
	"chr5": "5",
	"chr6": "6",
	"chr7": "7",
	"chr8": "8",
	"chr9": "9",
	"chr10": "10",
	"chr11": "11",
	"chr12": "12",
	"chr13": "13",
	"chr14": "14",
	"chr15": "15",
	"chr16": "16",
	"chr17": "17",
	"chr18": "18",
	"chr19": "19",
	"chrX": "20",
	"chrY": "21",
	"chrM": "22"
}



def read_pdata(filename, n_cols, header_flag):
    rows = 1000
    read_n_th_line=1
    skip = np.arange(rows)
    skip = np.delete(skip, np.arange(0, rows, read_n_th_line))
    col_to_read=set(range(0,n_cols))
    if header_flag:
            b=pd.read_table(filename, delim_whitespace=True, skiprows=skip, usecols=col_to_read, engine='c')# 
    else:
            b=pd.read_table(filename, delim_whitespace=True, skiprows=skip, usecols=col_to_read, engine='c', header=None)# 

    b.iloc[:, 0]=b.iloc[:, 0].map(chrdict) 

    n_rows=len(b.index)

    for i in b.columns:
        try:
            b[[i]] = b[[i]].astype(int)
        except:
            print("Error in converting data pandas dataframe to numbers")
            exit()

    n_array=b.to_numpy()
 
    return_list=[]
    for i in range(0,n_cols):
        tmp=n_array[:,i]
        return_list.append(tmp)
    return return_list
    

def chr_starts(peak_chr_ind):
    start_ind=[]
    end_ind=[]
    chr_bin_seq=[]

    prev=0
    ct=-1
    for indd in peak_chr_ind:
        ct+=1
        if prev!=indd:
            start_ind.append(ct)
            chr_bin_seq.append(indd)
        prev=indd
    end_ind=start_ind[1:]
    end_ind.append(len(peak_chr_ind)-1)
    return [chr_bin_seq, start_ind, end_ind]

def get_bin_Data(bed, peak, pind, pst, pen ):
    
    try:
        tmp_ind=bed.chr_list.index(pind)
    except:
        print("chr not found")
        exit()
        
    sv=bed.start_ind[tmp_ind]
    ev=bed.end_ind[tmp_ind]
    temp_numpy_st=bed.st[sv:ev]
    temp_numpy_enr = bed.enr[sv:ev]
    ps1=position = np.argmax(temp_numpy_st[sv:ev] >= pst)
    ps2=position = np.argmax(temp_numpy_st[sv:ev] >= pen)
    num_bins=round((pen-pst)/bed.bin_width)+1
    sample_vals=np.zeros([num_bins,1],dtype=int)
    sample_vals[0:ps2-ps1,0]=temp_numpy_enr[ps1:ps2]
    return sample_vals



==MAIN=====================================================================================
==MAIN=====================================================================================
  
sample = bedgraph(sample_file)
noise = bedgraph(noise_file)
peak = peak_class(peak_file)
print(peak.num_peaks)
stat_var=[]
pvalue=[]

for pi in range(0,peak.num_peaks):
    pind=peak.chrs[pi]
    pst=peak.st[pi]
    pen=peak.en[pi]
    print('peak',pind,pst,pen)
    
    sample_vals = get_bin_Data(sample, peak, pind, pst, pen )
    noise_vals = get_bin_Data(noise, peak, pind, pst, pen )
    '''
    a=range(0,len(sample_vals))
    b=range(0,len(noise_vals))

    plt.plot(a,sample_vals)
    plt.plot(b,noise_vals)
    plt.show()

    sb.distplot(sample_vals, hist=False, label = "x")
    sb.distplot(noise_vals, hist=False, label = "z")
    plt.legend()
    splt.show()
    '''
    relt=ks_2samp(sample_vals[:,0], noise_vals[:,0])
    stat_var.append(relt.statistic)
    pvalue.append(relt.pvalue)
    
    print(relt.statistic,relt.pvalue)
    
peak.write_pvalues(out_file,stat_var,pvalue)
    
    
    
    
    
    








