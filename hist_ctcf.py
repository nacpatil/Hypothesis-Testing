import numpy as np
import os
import pandas as pd
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import matplotlib.pyplot as splt
import seaborn as sb
from scipy import stats
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
#======Supporting Functions============================================================================

class bedgraph:
    def __init__(self,filename ):
        self.filename=filename
        data=read_pdata(filename, 4, True)
        [self.chrs, self.st, self.en, self.enr]= data
        [self.chr_list, self.start_ind, self.end_ind]=chr_starts(self.chrs)
        self.bin_width = data[2][0]-data[1][0]
        print(self.bin_width)
        self.chr_split_st=[[]]*25
        self.chr_split_enr=[[]]*25
        self.split_chr_data()
        
    def split_chr_data(self):
        ct=-1
        for chh in self.chr_list:
            ct+=1
            self.chr_split_st[chh]  = self.st[self.start_ind[ct] : self.end_ind[ct]] 
            self.chr_split_enr[chh] =self.enr[self.start_ind[ct] : self.end_ind[ct] ]


    def get_full_baspair(self):
        #nm=self.start_ind[1]
        #x=range(0,CTCFd.st[0:nm-1])
        full_list=[]
        for chrn in range(0,len(self.chr_list)):
            print('chrn',chrn,self.chr_list[chrn])
            stt=CTCFd.start_ind[chrn]
            enn=CTCFd.end_ind[chrn]-1

            x=[0]*(CTCFd.st[enn]+4)
            print(self.start_ind)
            print(self.chr_list)
            old=-1
            #acc=[]
            print(self.end_ind)
            for i in range(stt,enn+1):
                x[self.st[i]]=int(self.enr[i])
                if old>self.st[i]:
                    exit('Error: Issue with indexes')
                old=self.st[i]
                #cc.append(self.st[i])

            #if len(acc)<100:    
            #    print(acc)
            full_list.extend(x)
        return full_list


        
class peak_class:
    def __init__(self,filename ):
        data=read_pdata(filename, 5, False)
        [self.chrs, self.st, self.en, self.width, self.log]= data
        self.num_peaks=len(self.chrs)
        
    def write_pvalues(self,out_file,stat,pvalue, s_ratio):
        f=open(out_file,'w')
        f.write("#chr start end width log ks-statistic pvalue\n")
        for i in range(0,self.num_peaks):
            out_str =   str(self.chrs[i]) + " " +str(self.st[i])+ " "+str(self.en[i])+ " " + \
            str(self.width[i])+" "+str(self.log[i])+" "+str(stat[i]) + " "+str(pvalue[i]) +" "+str(s_ratio[i])+"\n"  
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
    "chr20": "20",
    "chr21": "21",
	"chrX": "22",
	"chrY": "23",
	"chrM": "24"
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
    #emp_numpy_st=bed.st[sv:ev]
    #temp_numpy_enr = bed.enr[sv:ev]
    temp_numpy_st=bed.chr_split_st[pind] 
    temp_numpy_enr=bed.chr_split_enr[pind] 
    
    #ps1 = np.argmax(temp_numpy_st >= pst)
    #ps2 = np.argmax(temp_numpy_st >= pen)
    ps1=np.searchsorted(temp_numpy_st, pst)
    ps2=np.searchsorted(temp_numpy_st, pen)

    num_bins=round((pen-pst)/bed.bin_width)+1
    sample_vals=np.zeros([num_bins,1],dtype=int)
    sample_vals[0:ps2-ps1,0]=temp_numpy_enr[ps1:ps2]
    return sample_vals



#=====INPUTS===========================================================================================
CTCF = "bin_new/JT_MmSc_WTCTCF912_.bedgraph"
CTCF = "bin_new/trial.bedgraph"

K27 = "bin_new/JT_MmSc_WTK27m912_.bedgraph"
K4 = "bin_new/JT_MmSc_WTK4m912_.bedgraph"
K9 = "bin_new/JT_MmSc_WTK9m912_.bedgraph"


#===MAIN===============================================================================================
#===MAIN===============================================================================================
CTCFd = bedgraph(CTCF)
#K27d = bedgraph(K27)
#K4d = bedgraph(K4)
#K9d = bedgraph(K9)

v=False

'''
sb.distplot(CTCFd.enr,fill=False, bins=5500, hist=True, label = "CTCF",kde=v)
sb.distplot(K27d.enr, bins=5500,fill=False, hist=True, label = "K27",kde=v)
sb.distplot(K4d.enr, bins=5500,fill=False, hist=True, label = "K4",kde=v)
sb.distplot(K9d.enr,  bins=5500,fill=False, hist=True, label = "K9",kde=v)
'''

'''
sb.histplot(data=CTCFd.enr,  log_scale=True, element="step", fill=False,color='red')
sb.histplot(data=K27d.enr, bins=5500, log_scale=True, element="step", fill=False,color="green")
sb.histplot(data=K9d.enr,  log_scale=True, element="step", fill=False,color="blue")
sb.histplot(data=K4d.enr, bins=5500, log_scale=True, element="step", fill=False,color="yellow")

plt.show()

'''

'''
nm=len(CTCFd.enr)
x=range(0,nm)
plt.plot(x,CTCFd.enr)
plt.show()
exit()
'''

full_list=CTCFd.get_full_baspair()
exit()
print(CTCFd.start_ind)
nm=CTCFd.start_ind[1]
print('nm',nm)
x=range(0,max(CTCFd.st[0:nm])+5)
y=np.zeros([max(CTCFd.st[0:nm])+5,1],dtype=int)
for i in range(0,nm-1):
    print('i',i)
    if True:#CTCFd.st[i]<nm:
        print(i)
        y[CTCFd.st[i]]=float(CTCFd.enr[i])
    else:
        break
print('ctcf',(CTCFd.start_ind))
#plt.plot(x,y)
#plt.show()
#splt.show()

std=np.std(y)
print('std',std)
mn=np.mean(y)
print('mean',mn)



N = int(nm/100)
T = 1.0 
#x = np.linspace(0.0, N*T, N, endpoint=False)
#y = np.sin(50.0 * 2.0*np.pi*x) + 0*0.5*np.sin(80.0 * 2.0*np.pi*x)
for i in range(0,len(y)):
    y[i]=i
yf = fft(y)
xf = fftfreq(N, T)[:N//2]
#>>> import matplotlib.pyplot as plt
print(yf)
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()


