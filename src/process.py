import pandas as pd
import csv
import numpy as np

def process_ped(pop, unrelpar_path, ped_path, samples):
    
    population_assignment = pd.read_csv(pop,sep="\t", header=None)
    population_assignment.columns = ['sample_id', 'superpop', 'subpop']
    populations = population_assignment['subpop'].unique().tolist()
    samples_pop={pop:population_assignment[population_assignment['subpop']==pop]['sample_id'].tolist() for pop in populations}
    samples_pop_inds={population_assignment.iloc[i,0]:population_assignment.iloc[i,2] for i in range(population_assignment.shape[0])}

    unrelated_parents=pd.read_csv(unrelpar_path, header=None)
    unrelated_parents_pops_df=pd.merge(population_assignment, unrelated_parents[0], left_on='sample_id', right_on=0)
    unrelated_parents_pops={pop:unrelated_parents_pops_df[unrelated_parents_pops_df['subpop']==pop].iloc[:,0].tolist() for pop in populations}

    ped_df=pd.read_csv(ped_path, sep='\t', dtype='str')
    ped_df=ped_df.iloc[:,0:6]
    ped_df.columns = ['family_id', 'individual_id', 'dad_id', 'mum_id', 'sex', 'affected']
    ped_df=ped_df[(ped_df['dad_id']!='0') & (ped_df['mum_id']!='0') & (ped_df['dad_id'].isna()==False) & (ped_df['mum_id'].isna()==False)] 
    ped={ped_df.iloc[i,1]:[ped_df.iloc[i,2], ped_df.iloc[i,3]] for i in range(ped_df.shape[0])}
    parents={ped_df.iloc[i,2]:ped_df.iloc[i,1] for i in range(ped_df.shape[0])}
    parents.update({ped_df.iloc[i,3]:ped_df.iloc[i,1] for i in range(ped_df.shape[0])})

    unrelated_parents_pop_index={}
    for pop in populations:
        unrelated_parents_pop_index[pop] = np.where(np.isin(samples,unrelated_parents_pops[pop]))[0]
    
    return (ped, parents, samples_pop, samples_pop_inds, unrelated_parents_pop_index, populations)

def store_gts(variant, variant_genotypes, an):
    vcf_AF = [variant.INFO.get('AF')] if isinstance(variant.INFO.get('AF'), float) else list(variant.INFO.get('AF'))

    # ignore 1/2 variants
    if vcf_AF[an]>0.5:
        gt_types=np.array([2 if i==[0,0,False] else 1 if i==[0,an+1,False] else 0 if i==[an+1,an+1,False] else 3 if -1 in i else 0 for i in variant_genotypes])
    else:
        gt_types=np.array([0 if i==[0,0,False] else 1 if i==[0,an+1,False] else 2 if i==[an+1,an+1,False] else 3 if -1 in i else 0 for i in variant_genotypes])

    return (vcf_AF, gt_types)

def calc_af_pop(populations, gt_types, unrelated_parents_pop_index):
    af_pop = {}
    for pop in populations:
        unrel_parents_gt_types = gt_types[unrelated_parents_pop_index[pop]][np.where(gt_types[unrelated_parents_pop_index[pop]]!=3)[0]]
        if len(unrel_parents_gt_types)>200:
            af_pop[pop]=np.sum(unrel_parents_gt_types)/(2*len(unrel_parents_gt_types))
        else:
            af_pop[pop]=0.0
    return af_pop
