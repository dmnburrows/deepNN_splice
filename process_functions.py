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


#=====================================================================
def SJ_gene_anot(rna_path, cell):
#=====================================================================

    """
    This function performs bedtools intersect on SJ.bed to give SJ.bed including the gene annotations. 
    
    Inputs:
    rna_path (str): path to rna SJ directory
    cell (str): celltype directory name in rna_path
    
    """
    # intersect with bed annotation to get gene names for each SJ
    run = f"""    
    tail -n +2 {rna_path}/{cell}/SJ.bed >  {rna_path}/{cell}/temp.bed
    bedtools intersect -s -wo -a {rna_path}/{cell}/temp.bed -b /cndd3/dburrows/DATA/annotations/genome/grcm38.p3/splice/genes.chr_1-19-X.bed > {rna_path}/{cell}/SJ.annotated.bed
    rm  {rna_path}/{cell}/temp.bed
    """
    get_ipython().run_cell_magic('bash', '', run)
    return()
 


#=====================================================================
def SJ_norm(rna_path, cell):
#=====================================================================

    """
    This function normalises SJ usage, removes SJs further apart than their genes and averages over SJs found across multiple genes.
    
    Inputs:
    rna_path (str): path to rna SJ directory
    cell (str): celltype directory name in rna_path
    
    Outputs:
    comb (df): directory of normalised SJs
    
    """

    curr = pd.read_csv(rna_path+cell + '/SJ.annotated.bed',sep='\t', header=None, index_col=False,low_memory=False,
                     names =['chromosome', 'SJ-point_start', 'SJ-point_end', 'id-SJ', 'id-unique', 'strand', 'cell',
                 'SJ_start', 'SJ_end','intron_motif', 'annotated', 'overhang', 'SJ_distance', 
                 'n_unq', 'n_mult', 'chr-anot', 'start-anot', 'end-anot', 'name-anot', 'geneid-anot', 'strand-anot', 'NA'])

    #remove any SJs that are further apart than their genes
    new = curr.loc[(curr['end-anot']-curr['start-anot']) - curr['SJ_distance'] > 0].copy()

    # Read in library size + gene annotation
    gene_mat = pd.read_csv(rna_path+cell+'/ReadsPerGene.out.tab', sep='\t', header=None, index_col=0,
                          names = ['counts', 'read1', 'read2'])
    wc = pd.read_csv(rna_path+cell+'/library_size.txt', header=None).loc[0].values[0]

    #Confirm all genes present 
    assert len(gene_mat.loc[new['geneid-anot'].unique()]) == len(new['geneid-anot'].unique()), 'Some genes found in SJ.bed not in RNA data'

    # Normalise -> 
    # Normalised SJ usage (NSJU) = (unique SJ reads) / RPKM of gene -> SJ usage relative to gene expression + library size
    # RPKM = (gene mapped reads x 10**9) / (gene length in bp * library size)
    ge =  (gene_mat.loc[new['geneid-anot']]['counts'].values) #raw gene expression
    le = (new['end-anot']-new['start-anot']).values #gene lengths

    #NSJU_rpkm
    rpkm =  (ge * 10**9) / (le  * wc)    
    new['NSJU_rpkm']=new['n_unq'] / rpkm
    #set all infs to 0
    new.loc[np.isinf(new['NSJU_rpkm']),'NSJU_rpkm'] = 0
    #CPMs
    new['CPM']=((new['n_unq'])/wc)*1e6
    #NSJU
    new['NSJU'] = (new['n_unq'] / ((ge/(le/1000))*wc))*1e9
    new.loc[np.isinf(new['NSJU']),'NSJU'] = 0


    #deal with SJs that overlap multiple genes -> take mean
    #=========================================================
    unq_ = np.unique(new['id-SJ'].values, return_counts=True)
    rep_ind = unq_[0][unq_[1] > 1] #repeated indeces
    unq_ind =  unq_[0][unq_[1] == 1] 
    new.index = new['id-SJ']
    new_unq = new.loc[unq_ind].copy()

    new_rep = new.loc[rep_ind].copy()
    assert len(new_rep) > len(rep_ind), 'repeated indeces are not repeated'

    #loop and merge
    newnew_rep = pd.DataFrame()
    for t,ind in enumerate(rep_ind):
        curr_ = new_rep.loc[ind]

        nsju = np.mean(curr_['NSJU'].values)
        rpkm = np.mean(curr_['NSJU_rpkm'].values)
        cpm = np.mean(curr_['CPM'].values)
        geneid_anot = "_".join(curr_['geneid-anot'].values)
        start_anot = "_".join(curr_['start-anot'].astype(str).values)
        end_anot = "_".join(curr_['end-anot'].astype(str).values)
        name_anot = "_".join(curr_['name-anot'].astype(str).values)

        hed = curr_.head(1).copy()
        hed['NSJU'] = nsju
        hed['CPM'] = cpm
        hed['NSJU_rpkm'] = rpkm
        hed['geneid-anot'] = geneid_anot
        hed['start-anot'] = end_anot
        hed['end-anot'] = end_anot
        hed['name-anot'] = name_anot
        newnew_rep = pd.concat([newnew_rep, hed])    


    #check that all rep_inds accounted for
    assert len(newnew_rep) == len(rep_ind), 'some repeat indeces not accounted for'

    comb = pd.concat([new_unq,newnew_rep])

    #check that there are no repeats
    assert len(np.unique(comb.index)) == len(comb), 'some repeat indeces remain'
    comb.reset_index(drop=True, inplace=True)
    return(comb)



#=====================================================================
def mc_process(path, cell, d, window_size, min_entries):
#=====================================================================

    """
    This function smoothes each subcluster trace and then takes the average across traces, to combine a mean smoothed trace for mc. 
    
    Inputs:
    path (str): path to mc bingraph directory
    cell (str): celltype directory name in path
    d (str): datatype
    window_size (int): bins to smoothe over by taking mean
    min_entries (int): minimum values in window required to be able to compute mean

    
    Outputs:
    out_df (df): dataframe of mean smoothed binned mc data
    
    """

    curr = pd.read_csv(glob.glob(path+cell+'/*'+d+'*merge.bingraph.raw')[0], sep='\t', 
                       low_memory=False)
    rmv_quot = [i[1:-1] for i in curr.columns[3:]]
    curr.columns =np.append(['chr', 'start', 'end'],rmv_quot)
    chr_sort = np.append(np.asarray(sorted(np.asarray(sorted(curr['chr'].unique())[:-1]).astype(int))).astype(str),'X')
    out_df = pd.DataFrame() #output for new values

    #Smoothe -> centered, sliding window, take mean, ignore nans and impute neighbouring values
    #split by chromosome!
    for ch in chr_sort:
        curr_chr = curr[curr['chr'] == ch]
        temp = curr_chr.copy()

        #compute sliding window over each sub-cluster first
        for c in curr_chr.columns[3:]:
            temp[c] = curr_chr[c].rolling(window=window_size, min_periods=min_entries, center=True).mean()

        out_df = pd.concat([out_df, temp])

    #take mean mc across subclusters -> impute nans by ignoring nans again
    out_df['mean_smoothed'] = np.nanmean(out_df[curr.columns[3:]],axis=1)
    out_df.index = np.arange(0,len(out_df))
    return(out_df)       


#=====================================================================
def atac_process(path, cell, d, window_size, min_entries):
#=====================================================================

    """
    This function smoothes each subcluster trace and then takes the average across traces, to combine a mean smoothed trace for mc. 
    
    Inputs:
    path (str): path to atac bingraph directory
    cell (str): celltype directory name in path
    d (str): datatype
    window_size (int): bins to smoothe over by taking mean
    min_entries (int): minimum values in window required to be able to compute mean

    
    Outputs:
    out_df (df): dataframe of mean smoothed binned atac data
    
    """

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
    return(out_df)
    
    

#=====================================================================
def epiflank_store(path, sub, binsize, n_bins, pos_enc):
#=====================================================================
    """
    This function computes the flanking epi context for a dataframe of SJs and stores them with the corresponding index.
    
    Inputs:
    path (str): path to bingraph directory
    sub (dataframe): cell dataframe of SJs
    binsize (int): size of each bin
    n_bins (int): number of bins covering entire flanking context
    pos_enc (list): list of positional encodings

    
    Outputs:
    li (list): list of SJ flanking contexts
    
    """
    li = []
    #load in bingraph
    curr = pd.read_csv(path, sep='\t', low_memory=False)
    curr['chr'] = curr['chr'].astype(str)

    #loop through each chr for current celltype
    for ch in curr['chr'].unique():
        curr_chr = curr[curr['chr'] == ch] #access chromosome bingraph
        curr_chr.index = np.arange(0,len(curr_chr)) #reset index
        assert np.unique(np.diff(curr_chr['start'])) == binsize, 'some bins not correct spacing'

        chr_sub = sub[sub['chromosome'] == 'chr'+ch] #access chromosome SJ
        chr_sub.index = np.arange(0,len(chr_sub)) #reset index

        #loop through each SJ for current chr
        for i in range(len(chr_sub)):
            #subset SJs by celltype
            sj = chr_sub.loc[i]
            ind = sj['id-unique'] #extract unique index

            #define start and end of flanking regions
            sj_point = sj['SJ-point_start']
            sj_bind = int(np.floor(sj_point/binsize))# find bin of SJ
            #check bin contains SJ
            assert sj_point <= curr_chr.loc[sj_bind]['end'] and sj_point >= curr_chr.loc[sj_bind]['start'], 'SJ bin position does not contain SJ'

            bind_start = int(sj_bind - (np.floor(n_bins/2))) #start bin
            bind_end = int(sj_bind + (np.floor(n_bins/2))) #end bin
            assert bind_end - bind_start == (n_bins-1), 'incorrect number of bins'
            vec = curr_chr.loc[bind_start:bind_end]['mean_smoothed'].values #vector of smoothed values
            assert len(vec) == n_bins, 'incorrect shape of extracted epi'

            #strand swap positions
            if sj['strand'] == '-': vec = vec[::-1]
            vec = pd.DataFrame(data=[vec], index=[ind], columns=pos_enc)
            li.append(vec)
    return(li)


#================================================================================
def epiflank_store_multidata(cell, data_l, sub, binsize, n_bins, pos_enc):
#================================================================================

    """
    This function computes the computes flanking epi context over a list of epi datasets and outputs them as separate lists.
    
    Inputs:
    cell (str): celltype
    data_l (list): list of datatypers
    sub (dataframe): dataframe of cell SJs
    binsize (int): size of each bin
    n_bins (int): number of bins covering entire flanking context
    pos_enc (list): list of positional encodings

    
    Outputs:
    cg_ (df): dataframe of CG SJ flanks
    ch_ (df): dataframe of CH SJ flanks
    atac_ (df): dataframe of ATAC SJ flanks
    
    """
    print(f'Running {cell}')
    assert len(sub) > 0, 'wrong cell name'
    
    #cellspecific lists
    cg_list, ch_list, atac_list, ind_l = [], [], [],[]

    #loop over each datatype
    for d in data_l:
        if 'CG' in d or 'CAC' in d: path = '/cndd2/dburrows/DATA/splice/snmcseq-mc_MOp_biccn/pseudobulk/gran1/' + cell + '/' + cell + '_' + d + '.norm.smoothe8.bingraph.raw'
        elif 'ATAC' in d: path = '/cndd2/dburrows/DATA/splice/snatacseq-atac_MOp_biccn/pseudobulk/gran1/' + cell + '/' + cell + '_' + d + '.norm-CPM.smoothe8.bingraph.raw'

        #compute flank list for each cell + datatype
        li = epiflank_store(path, sub, binsize, n_bins, pos_enc) 
        if 'CG' in d: cg_list = li
        elif 'CAC' in d: ch_list = li
        elif 'ATAC' in d: atac_list = li

    #cat into dataframes
    cg_ = pd.concat(cg_list)
    ch_ = pd.concat(ch_list)
    atac_ = pd.concat(atac_list)
    return(cg_, ch_, atac_)

#================================================================================
def run_epiflank_store_pool(flank, binsize, data_l, tot, name_l, cores):
#================================================================================

    """
    This function pools the celltype epiflank information together. 
    
    Inputs:
    flank (int): flanking bp length on one side, ie. total_flank = flank*2
    binsize (int): size of each bin
    data_l (list): list of datatypes
    tot (dataframe): dataframe of all SJs
    name_l (list): list of cellnames
    cores (int): number of cores
    
    """
    
    import sys
    import os
    import glob
    import pandas as pd
    from matplotlib import pyplot as plt
    import numpy as np
    from multiprocessing import Pool
    import multiprocessing

    #if flank not equally divisible, deal with remaining sequence
    #NEVER CLIP DATA, always add extra bin -> consistency across SJs
    n_bins = int(np.ceil((2*flank)/binsize))
    if n_bins%2 ==0: n_bins+=1 #make sure bin number is odd, to allow even number of bins either side

    #make sure 0 surrounded by even numbers on both sides
    pos_enc = np.arange(-1*(np.floor(n_bins/2)),(np.floor(n_bins/2)+1))
    assert len(pos_enc) == n_bins, 'incorrect nbins in positional encoding'
    assert sum(pos_enc > 0) == sum(pos_enc < 0), 'uneven numbers of bins for + or -'
    pos_enc = pos_enc.tolist()

    #declare output lists
    cg_list, ch_list, atac_list, ind_l = [], [], [],[]
    
    #define arglist
    sub_l = [tot[tot['cell'] == f] for f in name_l]
    arglist = [(cell, data_l, sub_l[f], binsize, n_bins, pos_enc) for f,cell in enumerate(name_l)]
    with multiprocessing.Pool(cores) as pool:
        results = pool.starmap(epiflank_store_multidata, arglist)
        
    cg_list, ch_list, atac_list = zip(*results)

    cg_df = pd.concat(cg_list)
    ch_df = pd.concat(ch_list)
    atac_df = pd.concat(atac_list)
    
    print("All data grouped.")
    
    assert len(cg_df.index) == len(tot), 'some SJs missing in CG file'
    assert len(ch_df.index) == len(tot), 'some SJs missing in CH file'
    assert len(atac_df.index) == len(tot), 'some SJs missing in ATAC file'
    
    tot.index = tot['id-unique'].values
    #make sure same order as SJ file
    cg_df = cg_df.loc[tot.index.values]
    ch_df = ch_df.loc[tot.index.values]
    atac_df = atac_df.loc[tot.index.values]

    assert sum(tot.index.values == cg_df.index.values) == len(tot.index), 'some SJs missing in CG file'
    assert sum(tot.index.values == ch_df.index.values) == len(tot.index), 'some SJs missing in CH file'
    assert sum(tot.index.values == atac_df.index.values) == len(tot.index), 'some SJs missing in ATAC file'
    print("All data pass QC.")
                    
    
    cg_df.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/CG.bingraph.raw', sep='\t', index=True)
    ch_df.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/CAC.bingraph.raw', sep='\t', index=True)
    atac_df.to_csv('/cndd2/dburrows/DATA/splice/model/processed_data/ATAC.bingraph.raw', sep='\t', index=True)

    print("All data saved.")
    

        

#==================================================================
def remove_nans_impute(curr, max_nan=10):
#==================================================================

    """
    Remove rows with more than `max_nan` NaN values and impute remaining NaNs by taking the mean of values 
    to the immediate left and right in the DataFrame `curr`. This process repeats, increasing the distance 
    from which the mean is calculated until no NaNs remain or all possible imputations have been attempted.

    Parameters:
    - curr: pd.DataFrame
        The input DataFrame containing epi info from splice junctions (SJs) with potential NaN values.
    - max_nan: int, default = 10
        The maximum number of NaN values allowed per row (SJ). Rows exceeding this threshold are removed.

    Returns:
    - tuple: (pd.DataFrame, np.ndarray)
        - The modified DataFrame with NaNs imputed where possible and rows with excess NaNs removed.
        - An array of indices for rows (SJs) removed due to exceeding the `max_nan` threshold.

    Raises:
    - AssertionError: If NaNs still exist after attempting imputation.

    """
    # Find number of NaNs per SJ and identify SJs to remove
    sum_nan = np.sum(np.isnan(curr.values), axis=1)
    remove_id = curr[sum_nan > max_nan].index.values
    sub_curr = curr[sum_nan <= max_nan].copy()  # Extract only SJs with fewer NaNs than max_nan

    # Keep looping until no NaNs remain
    count = 1
    while np.sum(np.isnan(sub_curr.values)) > 0:
        loc = np.where(np.isnan(sub_curr.values) == True)
        left = loc[0], loc[1] - count
        left[1][left[1] < 0] = 0
        right = loc[0], loc[1] + count
        right[1][right[1] >= curr.shape[1]] = curr.shape[1] - 1

        new_v = np.nanmean(np.vstack((sub_curr.values[left], sub_curr.values[right])), axis=0)  # Take mean either side, ignore NaNs
        sub_curr.values[loc] = new_v
        count += 1

    assert np.sum(np.isnan(sub_curr.values)) == 0, 'NaNs still left after imputation'

    return sub_curr, remove_id



        