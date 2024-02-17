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
import subprocess
import multiprocessing
import scanpy as sc

#Import your modules
#---------------------------------------
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm


#======================================
#pseudobulk granular celltypes
#NB only run once! or remove out fastqs
#======================================
def merge(y, cell):
    # access each file, untar -> append contents of r1 to file and r2 to file
    for x,g in enumerate(cell):
        run = f""" 
        tar -xf {in_path}{g}.fastq.tar --to-stdout {g}_R1.fastq.gz | zcat >> {out_path}{name_l[y]}/{name_l[y]}_R1.fastq
        tar -xf {in_path}{g}.fastq.tar --to-stdout {g}_R2.fastq.gz | zcat >> {out_path}{name_l[y]}/{name_l[y]}_R2.fastq
        """

        subprocess.run(run, shell=True)
        if int(x%(len(cell)/10)) == 0: print(f'Done {int(x/len(cell)*100)} of {name_l[y]}')
    gzip= f"""
    gzip {out_path}{name_l[y]}/{name_l[y]}_R1.fastq
    gzip {out_path}{name_l[y]}/{name_l[y]}_R2.fastq
    """
    subprocess.run(gzip, shell=True)
    print(f'Done {name_l[y]}')    
        
#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    save=f"""
    scp /cndd3/dburrows/CODE/deepNN_splice/workspace.py /cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/log.workspace.pseudobulk
    """
    subprocess.run(save, shell=True)
    in_path = '/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/fastq/'
    out_path = '/cndd3/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1/'    
    #Load metadata from pachter paper
    #load in metadata 
    meta = sc.read_h5ad('/cndd3/dburrows/DATA/public_datasets/splice//smartseq_mouse_p53-59_MOp_pachter21/metadata/isoform.h5ad')

    gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
    glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
    name_l = glu_l + gaba_l
    cell_l = []
    #change names to match 
    new_names = []
    for i in meta.obs['subclass_label'].values:
        if '/' in i or ' ' in i: 
            new = i.replace('/','-')
            new = new.replace(' ', '_')
            new_names.append(new)
        else: new_names.append(i)
    meta.obs['new_names'] = new_names

    for n in name_l:
        d = meta[meta.obs['new_names'] == n].obs['cell_id'].values
        cell_l.append(d)

    args_list = list(enumerate(cell_l))  # Create a list of tuples (index, cell_list)

    with multiprocessing.Pool(10) as pool:
        pool.starmap(merge, args_list)

    print("All tasks completed.")

    