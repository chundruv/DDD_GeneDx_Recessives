import pandas as pd

def calc_expected(x, lds, genes, classes, parents_populations, populations, N_haps, a, consequence_classes, unrel_parents):
    x=x[(np.isin(x['stable_id'], unrel_parents[0].tolist())) | (x['is_proband']==True)]

    for gene in genes:
        for c in classes:
    
            # subset to variants of parents in gene of class c
            x_subset = x[(x['gene.stable.id']==gene) & (x['canonical_vep_annotation_category']==c) & (x['is_proband']==False)]
    
            # loop over child_ids
            for i in x_subset['child_id'].unique():
                ind_subset=x_subset[x_subset['child_id']==i]
    
                # loop over parents
                for p in ind_subset['stable_id'].unique():
                    # get parent populations
                    # change column name for gdx, probably 0 if no header
                    if p in parents_populations['person_stable_id'].tolist():
                        parent_pop = parents_populations[parents_populations['person_stable_id']==p][2].tolist()[0]
                        # only continue if parent in a population of interest
                        if parent_pop in populations:
                            # if geno=2, both haplotypes have genotype
                            if ind_subset[ind_subset['stable_id']==p]['genotype'].max()==2:
                                h[c][gene][parent_pop]+=2
                            # if geno=1 and only 1 variant in parent is het => only 1 haplotype has genotype
                            elif len(ind_subset[ind_subset['stable_id']==p]['position'].unique())==1:
                                h[c][gene][parent_pop]+=1
                            # else multiple hets in parent
                            else:
                                # if child geno doesn't change => only 1 haplotype has genotype
                                if len(ind_subset[ind_subset['stable_id']==p]['child_genotype'].unique())==1:
                                    h[c][gene][parent_pop]+=1
                                # child geno changes. If other parent doesn't have multiple het or homalts => both haplotypes have genotype
                                elif len(ind_subset[ind_subset['stable_id']!=p]['child_genotype'].unique())<2 and ind_subset[ind_subset['stable_id']!=p]['genotype'].max()<2:
                                    h[c][gene][parent_pop]+=2
                                # else can't tell if both haplotypes have a genotype or not so we only count 1
                                else:
                                    h[c][gene][parent_pop]+=1
    
    f={c:{gene:{pop:h[c][gene][pop]/N_haps[pop] if N_haps[pop]!=0 else 0 for pop in populations} for gene in genes} for c in classes}
    f1={c:{gene:{pop:(h[c][gene][pop]+1)/(N_haps[pop]+2) if h[c][gene][pop]==0 else h[c][gene][pop]/N_haps[pop] if N_haps[pop]!=0 else 0 for pop in populations} for gene in genes} for c in classes}
    
    lmbda={c:{gene:{pop:{ld:0 for ld in lds} for pop in populations} for gene in genes} for c in consequence_classes}
    lmbda1={c:{gene:{pop:{ld:0 for ld in lds} for pop in populations} for gene in genes} for c in consequence_classes}
    for gene in genes:
        for pop in populations:
            for ld in lds:
                for c in ['lof/lof', 'missense/missense', 'synonymous/synonymous']:
                    lmbda[c][gene][pop][ld] = (1 - a[gene][pop][ld]) * (f[consequence_classes[c][0]][gene][pop] ** 2) + (a[gene][pop][ld] * f[consequence_classes[c][0]][gene][pop])
                    lmbda1[c][gene][pop][ld] = (1 - a[gene][pop][ld]) * (f1[consequence_classes[c][0]][gene][pop] ** 2) + (a[gene][pop][ld] * f1[consequence_classes[c][0]][gene][pop])
                for c in ['lof/missense']:
                    lmbda[c][gene][pop][ld] = (1 - a[gene][pop][ld]) * (2 * f[consequence_classes[c][0]][gene][pop] * f[consequence_classes[c][1]][gene][pop] * (1 - f[consequence_classes[c][0]][gene][pop]))
                    lmbda1[c][gene][pop][ld] = (1 - a[gene][pop][ld]) * (2 * f1[consequence_classes[c][0]][gene][pop] * f1[consequence_classes[c][1]][gene][pop] * (1 - f1[consequence_classes[c][0]][gene][pop]))
     return((lmbda, lmbda1))
