import pandas as pd

def calc_roh_overlap(genes, gene_pos, rohs_dir, populations, qual_filter=20, N_probands)
    lds=[0.2, 0.4, 0.6, 0.8]
    
    a_tmp={gene:{pop:{ld:0 for ld in lds} for pop in populations} for gene in genes}
    
    for pop in populations:
        for ld in lds:
            rohs=pd.read_csv(rohs_dir+'/'+populations[pop]+'_ld'+str(ld)+'.txt', header=None, sep=r'\s+', low_memory=False)
            rohs.columns=['RG', 'individual_id', 'Chromosome', 'Start', 'End', 'length', 'n_variants', 'qual']
            rohs=rohs[rohs['qual']>=qual_filter]
            rohs = rohs[np.isin(rohs['individual_id'], probands_populations[probands_populations[2]==pop][0].tolist())]
            for gene in genes:
                gene_df=pd.DataFrame([[genes_pos[gene].split(':')[0], genes_pos[gene].split(':')[1].split('-')[0], genes_pos[gene].split(':')[1].split('-')[1]]], columns=['Chromosome', 'Start', 'End'])
                gr1, gr2 = pr.PyRanges(rohs), pr.PyRanges(gene_df)
                gr = gr1.intersect(gr2)
                rohs_subset = gr.df
                if rohs_subset.shape[0]>0:
                    a_tmp[gene][pop][ld]+=len(rohs_subset['individual_id'].unique())
    
    
    a={gene:{pop:{ld:a_tmp[gene][pop][ld]/N_probands[pop] if N_probands[pop]!=0 else 0 for ld in lds} for pop in populations} for gene in genes}
    return(a)
