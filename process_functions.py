import numpy as np
import pandas as pd
import os
import glob
#==================================
def group_bw(path , mode, name_l):
#==================================

    """
    This function groups bigwigs by celltype:
    
    Inputs:
    path = path to bw
    mode = 'CAC', 'CG' or 'ATAC'
    name_l = list of names
    
    Outputs:
    dic
    """
    di = {'Data': mode}
    for n in name_l:
        if '_' in n: search = n.split('_')[0] + '*' +  n.split('_')[1]
        else: search = n
        curr_l = np.asarray([os.path.basename(i) for i in glob.glob(f'{path}/*{search}_*{mode}.bw')])
        
        #ignore C10_1_Sst_Chodl_CAC.bw
        if mode == 'CAC' and n=='Sst': curr_l = curr_l['C10_1_Sst_Chodl_CAC.bw'!=np.asarray(curr_l)]
        di.update({n: curr_l})
    return(di)
    
#=====================================================================
def STARSJ_to_bed(sj, cell, min_unq=5, min_overhang=5, min_sjdist=30):
#=====================================================================

    """
    This function converts STARs SJ.out.tab file to a bed file, and filters out SJs based on mapping statistics.
    
    Inputs:
        sj (pd dataframe): SJ.out.tab pd dataframe
        cell (str): string of celltype
        min_unq (int): minimum number of uniquely mappiung reads for inclusion
        min_overhang (int): minimum number of max-bp-overhang for inclusion
        min_sjdist (int): minimum bp distance between SJ donor and acceptor for inclusion
    
    Outputs:
        sj_comb (df): output bed file
    
    """
    
    #naming convention in SJ.out.tab
    # column 2: first base of the intron (1-based)
    # column 3: last base of the intron (1-based)

    #naming convention here
    # SJ-point_start: first base of the intron (0-based)
    # SJ-point_end: last base of the intron (0-based)

    #0 base for bed format -> start is inclusive, end is exclusive: e.g. AGCT, range A->C: 0,1,2,3; range A->T: 0,1,2,3,4
    #original bam -> start is inclusive, end is inclusive: A->C: 1,2,3
    #my new ranges will only contain a single bp; the first and last base of the intron


    
    #Convert strand to symbol
    #1 = +; 2 = - 
    sj.loc[sj['strand'] == 1, 'strand'] = '+'
    sj.loc[sj['strand'] == 2, 'strand'] = '-'

    #Filter Sjs with =< n reads
    sj_filt = sj.loc[sj['n_unq'] > min_unq]

    #SJ overhang filter > n
    sj_filt = sj_filt.loc[sj_filt['overhang'] > min_overhang]

    # Add in SJ distance, and filter minimum sjdist (minimum intron size)
    sj_filt['SJ_distance'] = sj_filt['SJ_end'] - sj_filt['SJ_start'].values
    sj_filt = sj_filt.loc[sj_filt['SJ_distance'] > min_sjdist]

    #set cell type
    sj_filt['cell'] = cell

    # Add in entry for donor + acceptor
    new_index = np.repeat(sj_filt.index.values, 2)
    # Use loc to duplicate the rows
    sj_filt = sj_filt.loc[new_index].reset_index(drop=True)

    #Label donor + acceptr
    sj_filt.loc[sj_filt.loc[::2].index, 'type'] = 'donor'
    sj_filt.loc[sj_filt.loc[::2].index + 1, 'type'] = 'acceptor'

    #set unique id and SJ id
    sj_filt['id-unique'] = sj_filt['cell'].astype(str).values + sj_filt['chromosome'].astype(str).values + sj_filt['SJ_start'].astype(str).values + sj_filt['SJ_end'].astype(str).values + sj_filt['strand'].astype(str).values + sj_filt['type'].astype(str).values
    sj_filt['id-SJ'] = sj_filt['chromosome'].astype(str).values + sj_filt['SJ_start'].astype(str).values + sj_filt['SJ_end'].astype(str).values + sj_filt['strand'].astype(str).values + sj_filt['type'].astype(str).values

    #Add in new SJ positions at each unique junction
    #====================================
    sj_filt['SJ-point_start'] = 'NA'
    sj_filt['SJ-point_end'] = 'NA'

    #Split into plus
    pl = sj_filt.loc[sj_filt['strand'] == '+'].copy()
    #splice donor
    pl.loc[pl['type'] == 'donor', 'SJ-point_start'] = pl.loc[pl['type'] == 'donor', 'SJ_start']-1
    pl.loc[pl['type'] == 'donor', 'SJ-point_end'] = pl.loc[pl['type'] == 'donor', 'SJ_start'] 
    #splice acceptor
    pl.loc[pl['type'] == 'acceptor', 'SJ-point_start'] = pl.loc[pl['type'] == 'acceptor', 'SJ_end']-1
    pl.loc[pl['type'] == 'acceptor', 'SJ-point_end'] = pl.loc[pl['type'] == 'acceptor', 'SJ_end'] 

    #Split into minus
    mi = sj_filt.loc[sj_filt['strand'] == '-'].copy()
    #splice donor
    mi.loc[mi['type'] == 'donor', 'SJ-point_start']  = mi.loc[mi['type'] == 'donor', 'SJ_end']-1
    mi.loc[mi['type'] == 'donor', 'SJ-point_end']  = mi.loc[mi['type'] == 'donor', 'SJ_end']

    #splice acceptor
    mi.loc[mi['type'] == 'acceptor', 'SJ-point_start']  = mi.loc[mi['type'] == 'acceptor', 'SJ_start']-1
    mi.loc[mi['type'] == 'acceptor', 'SJ-point_end']  = mi.loc[mi['type'] == 'acceptor', 'SJ_start']

    #Cat and re-order cols for bed
    sj_comb = pd.concat([pl, mi])
    sj_comb = sj_comb[['chromosome', 'SJ-point_start', 'SJ-point_end', 'id-SJ', 'id-unique', 'strand', 'cell',
             'SJ_start', 'SJ_end','intron_motif', 'annotated', 'overhang', 'SJ_distance', 
             'n_unq', 'n_mult']]
    return(sj_comb)