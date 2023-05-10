import pandas as pd
import numpy as np

def call_comphets(x, genes):
    chets= pd.DataFrame(columns=list(x.columns)+['functional_category'])
    chets=chets.astype({'is_proband':bool})
    chets=chets.astype({'is_proband':bool})
    # loop over genes
    for gene in genes:
        # subset rows which are in this gene, are hets, and the individual is a proband
        x_gene=x[(x['gene.stable.id']==gene) & (x['genotype']==1) & (x['is_proband']==True)]
        # loop over this subsetted dataframe
        for ind in np.unique(x_gene['individual_id']):
            # We can only pseudo-phase variants where one of the mum or dad is hom-ref. In all other cases it is uncertain which parent the allele was inherited from
            dad=np.where((x_gene['individual_id']==ind) & (x_gene['child_inheritance']=='dad'))[0]
            mum=np.where((x_gene['individual_id']==ind) & (x_gene['child_inheritance']=='mum'))[0]
    
            x_dad = x_gene.iloc[dad,]
            x_mum = x_gene.iloc[mum,]
            # Only count the compound hets if there is one variant with dad as het and mum as hom-ref, and another with vice versa
            if len(dad)>0 and len(mum)>0:
                # if both variants are LoFs
                if ('lof' in x_mum['canonical_vep_annotation_category'].tolist()) & ('lof' in x_dad['canonical_vep_annotation_category'].tolist()):
                    for xindex in range(x_mum[x_mum['canonical_vep_annotation_category']=='lof'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_mum[x_mum['canonical_vep_annotation_category']=='lof'].iloc[xindex])+['lof/lof']], columns=list(x.columns)+['functional_category'])])
                    for xindex in range(x_dad[x_dad['canonical_vep_annotation_category']=='lof'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_dad[x_dad['canonical_vep_annotation_category']=='lof'].iloc[xindex])+['lof/lof']], columns=list(x.columns)+['functional_category'])])
                # if one variant is LoF and the other missense or inframe
                if ('lof' in x_mum['canonical_vep_annotation_category'].tolist()) & ('missense/inframe' in x_dad['canonical_vep_annotation_category'].tolist()):
                    for xindex in range(x_mum[x_mum['canonical_vep_annotation_category']=='lof'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_mum[x_mum['canonical_vep_annotation_category']=='lof'].iloc[xindex])+['lof/missense']], columns=list(x.columns)+['functional_category'])])
                    for xindex in range(x_dad[x_dad['canonical_vep_annotation_category']=='missense/inframe'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_dad[x_dad['canonical_vep_annotation_category']=='missense/inframe'].iloc[xindex])+['lof/missense']], columns=list(x.columns)+['functional_category'])])
                if ('missense/inframe' in x_mum['canonical_vep_annotation_category'].tolist()) & ('lof' in x_dad['canonical_vep_annotation_category'].tolist()):
                    for xindex in range(x_mum[x_mum['canonical_vep_annotation_category']=='missense/inframe'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_mum[x_mum['canonical_vep_annotation_category']=='missense/inframe'].iloc[xindex])+['lof/missense']], columns=list(x.columns)+['functional_category'])])
                    for xindex in range(x_dad[x_dad['canonical_vep_annotation_category']=='lof'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_dad[x_dad['canonical_vep_annotation_category']=='lof'].iloc[xindex])+['lof/missense']], columns=list(x.columns)+['functional_category'])])
                # if both variants are missense or inframe
                if ('missense/inframe' in x_mum['canonical_vep_annotation_category'].tolist()) and ('missense/inframe' in x_dad['canonical_vep_annotation_category'].tolist()):
                    for xindex in range(x_mum[x_mum['canonical_vep_annotation_category']=='missense/inframe'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_mum[x_mum['canonical_vep_annotation_category']=='missense/inframe'].iloc[xindex])+['missense/missense']], columns=list(x.columns)+['functional_category'])])
                    for xindex in range(x_dad[x_dad['canonical_vep_annotation_category']=='missense/inframe'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_dad[x_dad['canonical_vep_annotation_category']=='missense/inframe'].iloc[xindex])+['missense/missense']], columns=list(x.columns)+['functional_category'])])
                # If both are synonymous
                if ('synonymous_variant' in x_mum['canonical_vep_annotation_category'].tolist()) and ('synonymous_variant' in x_dad['canonical_vep_annotation_category'].tolist()):
                    for xindex in range(x_mum[x_mum['canonical_vep_annotation_category']=='synonymous_variant'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_mum[x_mum['canonical_vep_annotation_category']=='synonymous_variant'].iloc[xindex])+['synonymous/synonymous']], columns=list(x.columns)+['functional_category'])])
                    for xindex in range(x_dad[x_dad['canonical_vep_annotation_category']=='synonymous_variant'].shape[0]):
                        chets=pd.concat([chets,pd.DataFrame([list(x_dad[x_dad['canonical_vep_annotation_category']=='synonymous_variant'].iloc[xindex])+['synonymous/synonymous']], columns=list(x.columns)+['functional_category'])])
    return(chets)
