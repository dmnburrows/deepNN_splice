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
    