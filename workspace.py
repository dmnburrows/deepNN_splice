#CONDA ENV base_conda (python 3.9.7)
#Import packages
#---------------------------------------
import sys
import os
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import scanpy as sc
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



#Load metadata from pachter paper
#load in metadata 
meta = sc.read_h5ad('/cndd3/dburrows/DATA/public_datasets/splice//smartseq_mouse_p53-59_MOp_pachter21/metadata/isoform.h5ad')

glu_cell = meta.obs['cell_id'][meta.obs['class_label'] == 'Glutamatergic'].values
gaba_cell = meta.obs['cell_id'][meta.obs['class_label'] == 'GABAergic'].values
cell_l = glu_cell, gaba_cell
name_l = ['GLU', 'GABA']

# access each file, untar -> append contents of r1 to file and r2 to file
in_path = '/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/fastq/'
out_path = '/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/coarse/'
for y,cell in enumerate(cell_l):
    for x,g in enumerate(cell):
        run = f""" 
        tar -xf {in_path}{g}.fastq.tar --to-stdout {g}_R1.fastq.gz | zcat >> {out_path}{name_l[y]}/{name_l[y]}_R1.fastq
        tar -xf {in_path}{g}.fastq.tar --to-stdout {g}_R2.fastq.gz | zcat >> {out_path}{name_l[y]}/{name_l[y]}_R2.fastq
        """
        
        subprocess.run(run, shell=True)
        if int(x%(len(cell)/10)) == 0: print(f'Done {int(x/len(cell)*100)} ')

    print("All tasks completed.")
    