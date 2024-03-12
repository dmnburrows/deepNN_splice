#CONDA ENV pytorch1 (python 3.9.7)
#Import packages
#---------------------------------------
import sys
import os
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
import process_functions as pf
import pyBigWig as bw


#Import your modules
#---------------------------------------
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm

# Define paths
#----------------------------------------------------------------------
l_code = '/Users/dominicburrows/Dropbox/PhD/Analysis/my_scripts/GitHub/'
l_data = '/Users/dominicburrows/Dropbox/PhD/analysis/Project/'
l_fig = '/Users/dominicburrows/Dropbox/PhD/figures/'

s_code = '/cndd3/dburrows/CODE/'
s_data = '/cndd3/dburrows/DATA/'
s_fig = '/cndd3/dburrows/FIGS/'


#PARS
window_size = 8 #8 -> change???
min_entries = 1 #minimum values in window to compute mean



#ATAC
# baseline -> merge.bingraph.raw
#============================================

#ATAC -> normalise, smoothe, take average
#============================================

gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
name_l = glu_l + gaba_l
data_l = ['ATAC']
#PARS
path = '/cndd2/dburrows/DATA/splice/snatacseq-atac_MOp_biccn/pseudobulk/gran1/'
for cell in name_l[:1]:
    for d in data_l[:1]:
        print(glob.glob(path+cell+'/*'+d+'*merge.bingraph.raw')[0])
        curr = pd.read_csv(glob.glob(path+cell+'/*'+d+'*merge.bingraph.raw')[0], sep='\t', 
                           low_memory=False)
        rmv_quot = [i[1:-1] for i in curr.columns[3:]]
        curr.columns =np.append(['chr', 'start', 'end'],rmv_quot)
        lib_size = np.asarray([curr[i].sum() for i in curr.columns[3:] ]) #library size
        chr_sort = np.append(np.asarray(sorted(np.asarray(sorted(curr['chr'].unique())[:-1]).astype(int))).astype(str),'X')
        out_df = pd.DataFrame() #output for new values
        
        #Smoothe -> centered, sliding window, take mean, ignore nans and impute neighbouring values
        #split by chromosome!
        for ch in chr_sort:
            curr_chr = curr[curr['chr'] == ch]
            temp = curr_chr.copy()
        
            #compute sliding window over each sub-cluster first
            for num,c in enumerate(curr_chr.columns[3:]):
                curr_norm = (curr_chr[c] /(lib_size[num])) * 1e6 #compute CPMs
                temp[c] = curr_norm.rolling(window=window_size, min_periods=min_entries, center=True).mean()
                
            out_df = pd.concat([out_df, temp])
        
        #take mean mc across subclusters -> impute nans by ignoring nans again
        out_df['mean_smoothed'] = np.nanmean(out_df[curr.columns[3:]],axis=1)
        out_df.index = np.arange(0,len(out_df))
        out_df.to_csv(f'{path}{cell}/{cell}_{d}.norm-CPM.smoothe8.bingraph.raw', sep='\t', index=False)

        
        