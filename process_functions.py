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
        print(ch, path)
        curr_chr = curr[curr['chr'] == ch] #access chromosome bingraph
        curr_chr.index = np.arange(0,len(curr_chr)) #reset index

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


        