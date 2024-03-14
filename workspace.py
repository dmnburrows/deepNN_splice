# intersect with bed annotation 

import sys
import os
import glob
import pandas as pd
import numpy as np

#use pool multiprocessing
import multiprocessing
import glob


def process_directory(c):
    print(f'Processing {c}')
    run = f"""   
    samtools view {out_path}/{c}/Aligned.sortedByCoord.out.bam | wc -l > {out_path}/{c}/library_size.txt
    """
    get_ipython().run_cell_magic('bash', '', run)
    print(f'Finished {c}')


#this section shows what default arguments will be run if just executing the script
if __name__ == "__main__":
    gaba_l = [ 'Vip', 'Pvalb', 'Sst', 'Lamp5']
    glu_l = ['L2-3_IT', 'L5_IT', 'L6_IT', 'L6_CT', 'L6b']
    name_l = glu_l + gaba_l

    out_path='/cndd2/dburrows/DATA/splice/smartseq-rna_MOp_biccn/pseudobulk/gran1'
    
    # Create a pool of 10 worker processes
    with multiprocessing.Pool(10) as pool:
        # Map the process_directory function to each item in dir_list
        pool.map(process_directory, name_l)

    print("All tasks completed.")
   