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

print(f'Python version {sys.version}')
print(f'Python path {sys.executable}')


#load in SJs
tot = pd.read_csv('/cndd2/dburrows/DATA/splice/model/processed_data/SJ.comb.bed', 
                 sep='\t', index_col=0)

#SMOOTHE + MEAN
#=====================
#PARS
window_size = 8 #8 -> change???
min_entries = 1 #minimum values in window to compute mean

flank = 5000 #flanking context either side
binsize=128


#STORE all epi data into dataframe
#======================================

gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
name_l = glu_l + gaba_l
data_l = ['CG', 'CAC', 'ATAC']


#set number of bins to include 
n_bins = int(np.ceil((2*flank)/binsize)) #always round up, bin number must be consistent across all SJs

#set positional encodings for dfs
pos_enc = np.arange(-1*(flank//binsize),(n_bins - flank//binsize))
pos_enc = pos_enc.tolist()
assert len(pos_enc) == n_bins, 'incorrect nbins in positional encoding'

#declare dfs
cg_, ch_, atac_ = [pd.DataFrame(columns=pos_enc) for _ in range(3)]


#for each celltype
for cell in name_l: 
    print(cell)
    #subset by celltype to reduce loading time on epi files
    sub = tot[tot['cell'] == cell]
    assert len(sub) > 0, 'wrong cell name'

    for d in data_l:
        if 'CG' in d or 'CAC' in d: 
            path = '/cndd2/dburrows/DATA/splice/snmcseq-mc_MOp_biccn/pseudobulk/gran1/' + cell + '/' + cell + '_' + d + '.norm.smoothe8.bingraph.raw'
        elif 'ATAC' in d: 
            path = '/cndd2/dburrows/DATA/splice/snatacseq-atac_MOp_biccn/pseudobulk/gran1/' + cell + '/' + cell + '_' + d + '.norm-CPM.smoothe8.bingraph.raw'
        
        #load in bingraph
        curr = pd.read_csv(path, sep='\t', low_memory=False)
        curr['chr'] = curr['chr'].astype(str)

        #loop through each chr for current celltype
        for ch in curr['chr'].unique():
            curr_chr = curr[curr['chr'] == ch] #access chromosome bingraph
            curr_chr.index = np.arange(0,len(curr_chr))
            
            chr_sub = sub[sub['chromosome'] == 'chr'+ch] #access chromosome SJ
            chr_sub.index = np.arange(0,len(chr_sub))
            
            #loop through each SJ for current chr
            for i in range(len(chr_sub)):
                #subset SJs by celltype
                sj = chr_sub.loc[i]

                ind = sj['id-unique'] #extract unique index

                #define start and end of flanking regions
                flank_start, flank_end = sj['SJ-point_start'] - flank, sj['SJ-point_end'] + flank
                bind_fs,  bind_fe = flank_start//binsize, ((flank_start//binsize) + n_bins).astype(int)
                vec = curr_chr[bind_fs:bind_fe]['mean_smoothed'].values #vector of smoothed values
                assert len(vec) == n_bins, 'incorrect shape of extracted epi'

                #strand swap positions
                if sj['strand'] == '-': vec = vec[::-1]
                # Creating the DataFrame
                vec = pd.DataFrame(data=[vec], index=[ind], columns=pos_enc)
                
                if 'CG' in d: cg_=pd.concat([cg_, vec])
                elif 'CAC' in d: ch_=pd.concat([ch_, vec])
                elif 'ATAC' in d: atac_=pd.concat([atac_, vec])
                    
    
cg_.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/CG.bingraph.raw', sep='\t', index=True)
ch_.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/CAC.bingraph.raw', sep='\t', index=True)
atac_.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/ATAC.bingraph.raw', sep='\t', index=True)

                  
assert len(cg_.index) == len(tot), 'some SJs missing in CG file'
assert len(ch_.index) == len(tot), 'some SJs missing in CH file'
assert len(atac_.index) == len(tot), 'some SJs missing in ATAC file'

assert sum(np.in1d(cg_.index, tot['id-unique'])) == len(cg_.index), 'some SJs missing in CG file'
assert sum(np.in1d(ch_.index, tot['id-unique'])) == len(ch_.index), 'some SJs missing in CH file'
assert sum(np.in1d(atac_.index, tot['id-unique'])) == len(atac_.index), 'some SJs missing in ATAC file'
        
