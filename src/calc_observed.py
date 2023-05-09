import pandas as pd
import numpy as np

def calc_obs(x, unrel_parents, unrel_probands, populations, consequence_classes, genes, probands_populations, chets):
    x=x[(np.isin(x['individual_id'], unrel_parents[0].tolist())) | (x['is_proband']==True)]
    x=x[(np.isin(x['individual_id'], unrel_probands[0].tolist())) | (x['is_proband']==False)]
    
    OB={c:{gene:{pop:0 for pop in populations} for gene in genes} for c in consequence_classes}
    varIDs={c:{gene:{pop:{} for pop in populations} for gene in genes} for c in consequence_classes}
    for gene in genes:
        for pop in populations:
            chets_in_gene=[]
            prev_count=[]
            for c in consequence_classes:
                x_subset = x[(x['gene.stable.id']==gene) & (x['population']==pop) & (np.isin(x['canonical_vep_annotation_category'],consequence_classes[c])) & (x['is_proband']==True)]
                for i in range(0, x_subset.shape[0]):
                    if x_subset.iloc[i,0] not in probands_populations['individual_id'].tolist():
                        continue
                    if c!='synonymous/synonymous':
                        if x_subset.iloc[i,0] in prev_count:
                            continue
                    if x_subset.iloc[i, 1]==1:
                        if len(chets[(chets['individual_id']==x_subset.iloc[i,0]) & (chets['variant_id']==x_subset.iloc[i,13]) & (chets['functional_category']==c) & (chets['gene.stable.id']==gene)])>0:
                            if x_subset.iloc[i,0] not in chets_in_gene:
                                varIDs[c][gene][pop][x_subset.iloc[i,0]]={var:1 for var in chets[(chets['individual_id']==x_subset.iloc[i,0]) & (chets['functional_category']==c) & (chets['gene.stable.id']==gene)]['variant_id'].tolist()}
                                OB[c][gene][pop]+=1
                                chets_in_gene += [x_subset.iloc[i,0]]
                                prev_count+=[x_subset.iloc[i,0]]
                    if x_subset.iloc[i, 1]==2:
                        if c not in ['lof/missense']:
                            varIDs[c][gene][pop][x_subset.iloc[i,0]]={x_subset.iloc[i,13]:2}
                            OB[c][gene][pop]+=1
                            prev_count+=[x_subset.iloc[i,0]]
    return((OB, varIDs))
