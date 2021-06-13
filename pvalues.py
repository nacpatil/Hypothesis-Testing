
import numpy as np
import os
import pandas as pd
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import matplotlib.pyplot as splt
import seaborn as sb
  
#=====INPUTS===========================================================================================
peak_file="peaks/sample_peak.bed"
sample_file="binned_data/sample_.bedgraph"
noise_file="binned_data/noise_.bedgraph"
output_file = "pvalues_peak_file.bed"



#======Supporting Functions============================================================================

  
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
    
peak.write_pvalues(output_file,stat_var,pvalue)
    
    
    
    
    
    








