import pandas as pd
import numpy as np
import pyranges as pr

def calc_roh_overlap(genes, genes_pos, rohs_dir, populations, N_probands, probands_populations, qual_filter=20):
    """
    This function calculates the proportion of probands with a ROH overlapping each gene
    input:
    genes - list of genes
    genes_pos - dictionary of genes and gene positions 
    rohs_dir - directory with ROH files
    populations - dictionary of populations with samples in those populations
    N_probands - dictionary with number of number of probands per population
    probands_populations - dataframe of probands and their population assignment
    qual_filter - quality filter to use, default 20
    output:
    a - dictionary of genes with the proportion of probands with ROH overlapping per population
    """
    lds=[0.2, 0.4, 0.6, 0.8]
    
    a_tmp={gene:{pop:{ld:0 for ld in lds} for pop in populations} for gene in genes}
    
    for pop in populations:
        for ld in lds:
            rohs=pd.read_csv(rohs_dir+'/'+populations[pop]+'_ld'+str(ld)+'.txt', header=None, sep=r'\s+', low_memory=False, comment='#')
            rohs.columns=['RG', 'individual_id', 'Chromosome', 'Start', 'End', 'length', 'n_variants', 'qual']
            rohs=rohs[rohs['qual']>=qual_filter]
            rohs = rohs[np.isin(rohs['individual_id'], probands_populations[probands_populations['subpop']==pop][0].tolist())]
            for gene in genes:
                gene_df=pd.DataFrame([[genes_pos[gene].split(':')[0], genes_pos[gene].split(':')[1].split('-')[0], genes_pos[gene].split(':')[1].split('-')[1]]], columns=['Chromosome', 'Start', 'End'])
                gr1, gr2 = pr.PyRanges(rohs), pr.PyRanges(gene_df)
                gr = gr1.intersect(gr2)
                rohs_subset = gr.df
                if rohs_subset.shape[0]>0:
                    a_tmp[gene][pop][ld]+=len(rohs_subset['individual_id'].unique())
    
    
    a={gene:{pop:{ld:a_tmp[gene][pop][ld]/N_probands[pop] if N_probands[pop]!=0 else 0 for ld in lds} for pop in populations} for gene in genes}
    return(a)
