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
import torch
from enformer_pytorch import Enformer, GenomeIntervalDataset
import polars as pl

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

print('running')

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
data_l = ['CG', 'CAC', 'ATAC']
gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
name_l = glu_l + gaba_l


#Store flanking epi for each SJ in one file
#==========================================
pf.run_epiflank_store_pool(flank, binsize, data_l, tot, name_l, cores=len(name_l))
