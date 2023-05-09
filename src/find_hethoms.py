import numpy as np

def is_hethom(samples, samples_pop, hethoms, gt_types):
    """
    Get list of samples with a het and hom-alt in variant
    Args:
        samples (list str): list of sample IDs
        samples_pop (dict): dictionary of samples in each population, keys = populations (str), values = sample IDs (list str)
        hethoms (list str): list of sample IDs with a het or hom-alt in another variant in the gene (no point double counting)
        gt_types (list int): list of genotypes 0,1,2 (3=missing), order is sample as samples list
    Returns:
        list (str): list of sample IDs with a het or hom-alt in this variant
    """
    hh = []
    # for each hets and hom alts
    for i in np.where((gt_types>0) & (gt_types<3))[0]:
        if np.any(samples[i] in samples_pop[pop] for pop in samples_pop):
            # Add individual to list
            if samples[i] not in hethoms:
                hh += [samples[i]]
    return hh

def find_hethoms(af_pop, gnomad_af, samples, samples_pop, gt_types):
    """
    Get list of samples with a het and hom-alt in variant
    Args:
        af_pop (dict): dictionary of populations (keys) and allele frequencies in DDD (values)
        gnomad_af (list): list of allele freuquencies in gnomAD v2 exomes
        # not using this atm .... gnomad_genomes_af (list): list of allele freuquencies in gnomAD v2 genomes
        samples (list str): list of sample IDs
        samples_pop (dict): dictionary of samples in each population, keys = populations (str), values = sample IDs (list str)
        gt_types (list int): list of genotypes 0,1,2 (3=missing), order is sample as samples list
    Returns:
        list (str): list of sample IDs with a het or hom-alt in this variant if it passes allele frequency thresholds
    """
    hethoms=[]

    if np.all([f <= 0.01 for f in list(af_pop.values())]):
        gnomad_af=[i for i in gnomad_af if i!='']
        #gnomad_genomes_af=[i for i in gnomad_genomes_af if i!='']

    # if AFs in gnomAD
        if len(gnomad_af)>0:# or len(gnomad_genomes_af)>0:
        # if all gnomAD AFs < 0.01
            if np.all([float(i)<=0.01 for i in gnomad_af]):# and np.all([float(i)<=0.01 for i in gnomad_genomes_af]):
        # if hom-alt or het
                tmp=is_hethom(samples, samples_pop, hethoms, gt_types)
                if len(tmp)!=0:
                    hethoms+=tmp
            
        # if AFs not in gnomad, just go by the parents AFs
        else:
            tmp=is_hethom(samples,samples_pop, hethoms, gt_types)
            if len(tmp)!=0:
                hethoms+=tmp   
    return hethoms
