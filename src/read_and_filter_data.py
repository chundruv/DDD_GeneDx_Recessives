import pandas as pd
import numpy as np

def read_sample_lists(args):
    """
    Read in lists of unrelated probands and parents, and the population assignments.

    input:
    * args
    output:
    * unrel_probands - list of unrelated probands
    * unrel_parents - list of unrelated unaffected parents
    * parents_populations - dataframe of unrelated parents and the population assignments
    * probands_populations - dataframe of unrelated probands and the population assignments
    * population_table - all population assignments
    """
    unrel_probands=pd.read_table(args.unrel_probands, header=None)
    unrel_probands[0]=unrel_probands[0].astype('str')
    unrel_parents=pd.read_table(args.unrel_parents, header=None)
    unrel_parents[0]=unrel_parents[0].astype('str')

    pedigree=pd.read_csv(args.pedfile, sep='\t', dtype=str)
    pedigree=pedigree.iloc[:,range(0,6)]
    pedigree.columns=["family_id", "individual_id", "dad_id", "mum_id", "sex", "affected"]
    unrel_probands=pd.merge(unrel_probands[0], pedigree['individual_id'].drop_duplicates(), left_on=0, right_on='individual_id')[[0]]
    unrel_parents=pd.concat([pd.merge(unrel_parents[0], pedigree['dad_id'].drop_duplicates(), left_on=0, right_on='dad_id')[0], pd.merge(unrel_parents[0], pedigree['mum_id'].drop_duplicates(), left_on=0, right_on='mum_id')[0]])
    population_table=pd.read_csv(args.popfile, header=None, sep='\t')
    population_table.columns=['individual_id','pop', 'subpop']
    population_table['individual_id']=population_table['individual_id'].astype('str')

    parents_populations=pd.merge(population_table, unrel_parents, left_on='individual_id', right_on=0)
    probands_populations=pd.merge(population_table, unrel_probands, left_on='individual_id', right_on=0)

    return((unrel_probands, unrel_parents, parents_populations, probands_populations, population_table))

def read_and_filter_data(args, unrel_parents, population_table):
    """
    read in output from parse.py with the hets and hom-alts
    input:
    * args
    * unrel_parents - list of unrelated unaffected parents
    * population_table - dataframe of population assignments
    output:
    * x - dataframe of genotypes of interest
    """
    x=pd.read_csv(args.input_dir+'/chr'+str(args.chrom)+'_'+str(args.g1)+'_'+str(args.g2)+'_recessive_candidate.txt', sep='\t', low_memory=False)
    x['stable_id']=x['stable_id'].astype('str')

    if args.idmap==None:  
        x=pd.merge(x, population_table[['individual_id', 'subpop']], right_on='individual_id', left_on='stable_id')
    else:
        idmap=pd.read_csv(args.idmap, sep='\t')
        idmap.columns=['stable_id', 'individual_id']
        x=pd.merge(x, idmap)
        x=pd.merge(x, population_table[['individual_id', 'subpop']], on='individual_id')

    x['population']=x['subpop']

    x=x[( ((x['is_proband']==True) & (x['dad_genotype'].isna()==False) & (x['mum_genotype'].isna()==False)) | ((x['is_proband']==False) & (x['child_genotype'].isna()==False)) )]
    x['size']=x['ref'].str.len() - x['alt'].str.len()
#    x=pd.concat([x[x['is_proband']==True], pd.merge(x[x['is_proband']==False], unrel_parents, left_on='individual_id', right_on=0)[x.columns]])
    return(x)

def consequence_filtering(x, cadd_filter, cadd_indel_filter, revel_filter, varity_filter, moipred_filter, clinpred_filter, polyphen_filter, synsplice_filter):
    """
    Filter genotypes
    - For synonymous variants we keep only variants with max_spliceai<0.1
    - For LoF variants we take only LOFTEE HC variants
    - For the LOFTEE LC variants, we filter on CADD and include it in the functional category
    - For Missense variants we filter on CADD, PolyPhen2, ClinPred, MOIpred recessive probability, REVEL, and VARITY_ER. We require the variant to pass >70% of the available annoations
    - For inframe indels we filter on CADD

    input:
    * x - dataframe of genotypes to filter
    * cadd_filter - SNV CADD filter threshold, default=24.18
    * cadd_indel_filter - Indel CADD filter threshold, default=17.34
    * revel_filter - REVEL filter threshold, default=0.36
    * varity_filter - VARITY_ER filter threshold, default=0.25
    * moipred_filter - MOIpred recessive probability filter threshold, default=0.11
    * clinpred_filter - ClinPred filter threshold, default=0.53
    * polyphen_filter - PolyPhen2 filter threshold, default=0.59
    * synsplice_filter - SpliceAI filter threshold, default=0.8

    output:
    * x - filtered dataframe
    """
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
    
    missense['PolyPhen']=missense['PolyPhen'].replace('[a-z]*[()_]', '', regex=True).astype('float64')

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
    """
    Read in data and filter
    input:
    * args
    output:
    * x - filtered dataframe of genotypes of interest
    * parents_populations - dataframe of unrelated unaffected parents and their population assignments
    * probands_populations - dataframe of unrelated probands and their population assignments
    * populations - list of populations in the data
    * N_haps - number of unrelated parent haplotypes = 2*number of unrelated parents
    * N_probands - number of unrelated probands
    * unrel_parents - list of unrelated unaffected parents
    * unrel_probands - list of unrelated probands
    """
    unrel_probands, unrel_parents, parents_populations, probands_populations, population_table = read_sample_lists(args)
    x = read_and_filter_data(args, unrel_parents, population_table)
    x = consequence_filtering(x, args.cadd_threshold, args.cadd_indel_threshold, args.revel_threshold, args.varity_threshold, args.moipred_threshold, args.clinpred_threshold, args.polyphen_threshold, args.synsplice_threshold)
    populations={i[1]['subpop']:i[1]['pop'] for i in population_table[['pop','subpop']].drop_duplicates().iterrows() if i[1]['subpop'].endswith('OTH')==False}
    N_haps={pop:2*parents_populations[parents_populations['subpop']==pop].shape[0] for pop in populations}
    N_probands={pop:probands_populations[probands_populations['subpop']==pop].shape[0] for pop in populations}
    return((x, parents_populations, probands_populations, populations, N_haps, N_probands, unrel_parents, unrel_probands)) 
