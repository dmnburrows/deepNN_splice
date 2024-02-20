import sys
import os
import glob
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import process_functions as pf
import multiprocessing
import subprocess

#Import your modules

#---------------------------------------
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm

def copy(filename, mode, cell):
    
    run=f"""
    scp /cndd3/dburrows/CODE/deepNN_splice/batch_copy-bw.py {out_path}/log.workspace.copy_bw
    echo Running {cell} for datatype {mode}
    scp {in_path}/{filename} {out_path}/{mode}/pseudobulk/gran1/{cell}/{filename}
    echo Done 
    """
    subprocess.run(run, shell=True)

if __name__ == '__main__':

    # filter reads quality? // trim_galore
    out_path='/cndd2/dburrows/DATA/splice/'
    in_path = '/cndd/emukamel/enformer/CEMBA_MOp_mc_bigwigs/bigwig'

    #group mc data
    gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
    glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
    name_l = glu_l + gaba_l

    #function to grab bw files
    cg = pf.group_bw(in_path, 'CG', name_l)
    cg_l = np.concatenate(np.asarray([cg[i] for x,i in enumerate(cg.keys()) if x>0]))
    ch = pf.group_bw(in_path, 'CAC', name_l)
    ch_l = np.concatenate(np.asarray([ch[i] for x,i in enumerate(ch.keys()) if x>0]))
    atac = pf.group_bw(in_path, 'ATAC', name_l)
    atac_l = np.concatenate(np.asarray([atac[i] for x,i in enumerate(atac.keys()) if x>0]))


    filename_l = np.append(cg_l, np.append(ch_l, atac_l))
    mode_l, cell_l = [],[]
    for x,i in enumerate(filename_l):
        if 'CG' in (i.split('.bw')[0]).split('_')[-1] or 'CAC' in (i.split('.bw')[0]).split('_')[-1] : mode_l.append('snmcseq-mc_MOp_biccn')
        elif 'ATAC' in (i.split('.bw')[0]).split('_')[-1] : mode_l.append('snatacseq_MOp_biccn')
        if 'Vip_' in i: cell_l.append('Vip')
        if 'Pvalb' in i: cell_l.append('Pvalb')
        if 'Sst' in i: cell_l.append('Sst')
        if 'Lamp5' in i: cell_l.append('Lamp5')
        if 'L2-3_IT' in i: cell_l.append('L2-3_IT')
        if 'L5_IT' in i: cell_l.append('L5_IT')
        if 'L6_IT' in i: cell_l.append('L6_IT')
        if 'L6_CT' in i: cell_l.append('L6_CT')
        if 'L6b' in i: cell_l.append('L6b')
    args_list = list(zip(filename_l, mode_l, cell_l))
    
    with multiprocessing.Pool(10) as pool:
        pool.starmap(copy, args_list)

    print("All tasks completed.")
    