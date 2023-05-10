import pandas as pd
import numpy as np

def read_sample_lists(args):
    unrel_probands=pd.read_table(args.unrel_probands, header=None, sep=r'\s+')
    unrel_parents=pd.read_table(args.unrel_parents, header=None, sep=r'\s+')
    unaffected_parents=pd.read_table(args.unaff_parents, header=None, sep=r'\s+')
    fail_qc=pd.read_table(args.qcfail, header=None, sep=r'\s+')

    unrel_parents=unrel_parents[np.isin(unrel_parents[0], unaffected_parents[0])]
    unrel_parents=unrel_parents[np.isin(unrel_parents[0], fail_qc[0])==False]
    unrel_probands=unrel_probands[np.isin(unrel_probands[0], fail_qc[0])==False]

    population_table=pd.read_csv(args.popfile, header=None, sep=r'\s+',low_memory=False)
    population_table.columns=['individual_id','pop', 'subpop']

    parents_populations=pd.merge(population_table, unrel_parents, left_on='individual_id', right_on=0)
    probands_populations=pd.merge(population_table, unrel_probands, left_on='individual_id', right_on=0)

    return((unrel_probands, unrel_parents, parents_populations, probands_populations, population_table))

def read_and_filter_data(args, unrel_probands, unrel_parents, population_table):
    pedigree=pd.read_csv(args.pedfile, sep=r'\s+', low_memory=False)
    unrel_probands=unrel_probands[np.isin(unrel_probands[0],pedigree['individual_id'].tolist())]
    unrel_parents=unrel_parents[((np.isin(unrel_parents[0],pedigree['dad_id'].tolist())) | (np.isin(unrel_parents[1],pedigree['mum_id'].tolist())))]

    x=pd.read_csv(args.input_dir+'/chr'+str(args.chrom)+'_'+str(args.g1)+'_'+str(args.g2)+'_recessive_candidate.txt', sep='\t', low_memory=False)

    x=pd.merge(x, population_table[['individual_id', 'subpop']], right_on='individual_id', left_on='stable_id')
    x['population']=x['subpop']
    x['individual_id']=x['individual_id'].astype('str')
    unrel_parents[0]=unrel_parents[0].astype('str')
    unrel_probands[0]=unrel_probands[0].astype('str')

    x=x[( ((x['is_proband']==True) & (x['dad_genotype'].isna()==False) & (x['mum_genotype'].isna()==False)) | ((x['is_proband']==False) & (x['child_genotype'].isna()==False)) )]
    x['size']=x['ref'].str.len() - x['alt'].str.len()
    x=pd.concat([x[x['is_proband']==True], pd.merge(x[x['is_proband']==False], unrel_parents, left_on='individual_id', right_on=0)[x.columns]])
    return(x)

def consequence_filtering(x, cadd_filter, cadd_indel_filter, revel_filter, varity_filter, moipred_filter, clinpred_filter, polyphen_filter, synsplice_filter):
    consequence_categories = {"lof":["transcript_ablation", "splice_donor_variant","splice_acceptor_variant", "stop_gained","frameshift_variant", "stop_lost"],
    "missense/inframe":["start_lost", "inframe_insertion", "inframe_deletion", "missense_variant", "transcript_amplification", "protein_altering_variant", "splice_region_variant"],
    "synonymous_variant":["synonymous_variant"]}
    for i in consequence_categories:
        x.loc[np.isin(x['canonical_vep_annotation'], consequence_categories[i]),'canonical_vep_annotation_category']=i
        x.loc[np.isin(x['worst_vep_annotation'], consequence_categories[i]),'worst_vep_annotation_category']=i

    x.loc[x['canonical_vep_annotation'].isna(),'canonical_vep_annotation_category']=x.loc[x['canonical_vep_annotation'].isna(),'worst_vep_annotation_category']
    x.loc[x['canonical_vep_annotation'].isna(),'canonical_vep_annotation']=x.loc[x['canonical_vep_annotation'].isna(),'worst_vep_annotation']
    lof=x[(x['canonical_vep_annotation_category']=='lof') & (x['LOFTEE']=='HC')].copy()
    missense=x[(x['canonical_vep_annotation_category']=='missense/inframe') & (x['size']==0)].copy()
    inframe=x[((x['canonical_vep_annotation_category']=='missense/inframe') & (x['size']!=0))].copy()
    lclof=x[(x['canonical_vep_annotation_category']=='lof') & ((x['LOFTEE']=='LC')| (x['LOFTEE'].isna()))].copy()
    synonymous=x[(x['canonical_vep_annotation_category']=='synonymous_variant') & ((x['spliceai']<0.1) | (x['spliceai'].isna()))].copy()
    synsplice=x[(x['canonical_vep_annotation_category']=='synonymous_variant') & (x['spliceai']>=synsplice_filter)].copy()

    missense['CADD_pass']=False
    missense['REVEL_pass']=False
    missense['VARITYER_pass']=False
    missense['PolyPhen_pass']=False
    missense['ClinPred_pass']=False
    missense['MOIpredRP_pass']=False
    
    missense.loc[(missense['CADD_phred']>=cadd_filter),'CADD_pass']=True
    missense.loc[(missense['REVEL']>=revel_filter),'REVEL_pass']=True
    missense.loc[(missense['VARITYER_LOO']>=varity_filter),'VARITYER_pass']=True
    missense.loc[(missense['PolyPhen']>=polyphen_filter),'PolyPhen_pass']=True
    missense.loc[(missense['clinpred']>=clinpred_filter),'ClinPred_pass']=True
    missense.loc[(missense['moipred_rp']>=moipred_filter),'MOIpredRP_pass']=True
    
    missense['filter_count']=missense[['CADD_pass', 'REVEL_pass', 'VARITYER_pass', 'PolyPhen_pass', 'ClinPred_pass', 'MOIpredRP_pass']].sum(axis=1)
    missense['na_count']=missense[['CADD_phred','REVEL', 'VARITYER_LOO', 'PolyPhen', 'clinpred', 'moipred_rp']].isna().sum(axis=1)
    missense=missense[(((missense['filter_count']/(6-missense['na_count']))>=(0.7)) | (missense['filter_count']/(6-missense['na_count'])).isna())]
    missense=missense[inframe.columns]
    inframe=inframe[(inframe['CADD_phred']>=cadd_indel_filter) | (inframe['CADD_phred'].isna())]
    lclof=lclof[(lclof['CADD_phred']>=cadd_filter) | (lclof['CADD_phred'].isna())]
    lclof['canonical_vep_annotation_category']='missense/inframe'
    synsplice['canonical_vep_annotation_category']='missense/inframe'
    
    x=pd.concat([synonymous,missense])
    x=pd.concat([x,inframe])
    x=pd.concat([x,lclof])
    x=pd.concat([x,synsplice])
    x=pd.concat([x,lof])
    
    x=x.reset_index(drop=True)
    return(x)

def read_file(args):
    unrel_probands, unrel_parents, parents_populations, probands_populations, population_table = read_sample_lists(args)
    x = read_and_filter_data(args, unrel_probands, unrel_parents, population_table)
    x = consequence_filtering(x, args.cadd_threshold, args.cadd_indel_threshold, args.revel_threshold, args.varity_threshold, args.moipred_threshold, args.clinpred_threshold, args.polyphen_threshold, args.synsplice_threshold)
    populations={i[1]['subpop']:i[1]['pop'] for i in population_table[['pop','subpop']].drop_duplicates().iterrows() if i[1]['subpop'].endswith('OTH')==False}
    N_haps={pop:2*parents_populations[parents_populations['subpop']==pop].shape[0] for pop in populations}
    N_probands={pop:probands_populations[probands_populations['subpop']==pop].shape[0] for pop in populations}
    return((x, parents_populations, probands_populations, populations, N_haps, N_probands, unrel_parents, unrel_probands)) 
