import pandas as pd
import numpy as np
import scipy.stats

def calc_pab(AD, genotypes, an, depths):
    """
    Calculate the binomial p-value of the allelic depth for each SNV heterozygote to use as a filter. For the homozygotes return p=0.5 which will pass all the filters.
    Args:
        AD (np.array list int): array of lists with the allelic depth for each individual in order of sample list
        genotypes (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        an (int): the allele number
        depths (list int): list of depths for each individual. v.gt_depths in cyvcf2
    Returns:
        list float: list of p-values for each individual
    """
    return np.array([scipy.stats.binom_test(AD[i][an], depths[i], 0.5, alternative='less') if genotypes[i][0]!=genotypes[i][1] and an in genotypes[i] else 0.5 for i in range(len(depths))])


def calc_ab(AD, genotypes, an):
    """
    Calculate the allelic balance for each indel heterozygote to use as a filter. For the homozygotes return ab=0.5 which will pass all the filters.
    Args:
        AD (np.array list int): array of lists with the allelic depth for each individual in order of sample list
        genotypes (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        an (int): the allele number
    Returns:
        list float: list of allelic balances for each individual
    """
    return np.array([0.5 if genotypes[i][0]==genotypes[i][1] else AD[i][an]/(AD[i][an]+AD[i][0]) if (AD[i][an]+AD[i][0])!=0 else 0 for i in range(len(genotypes))])    


def apply_SNV_genotype_QC(v_g, quals, depths, GQ, DP, pAB, an, pab_obs):
    """
    Apply SNV genotype QC, masking genotypes which fail the QC
    Args:
        v_g (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        quals (list int): the allele number. v.gt_quals in cyvcf2
        depths (list int): list of depths for each individual. v.gt_depths in cyvcf2
        GQ (int): genotype quality threshold
        DP (int): depth threshold
        pAB (float): binomial of allelic depths threshold
        an (int): allele number
        pab_obs (list float): the observed binomial of allelic depths
    Returns:
        tuple: (variant, pass rate) list of genotypes (format list [a1,a2,phased]) masked if failing thresholds, and the pass rate of the QC
    """
    failed=0
    for i in np.where((quals<GQ)|(depths<DP)|(pab_obs<pAB))[0]:
        v_g[i]=[-1]*2 + [False]
        failed+=1
    return (v_g,(len(quals)-failed)/len(quals))


def apply_indel_genotype_QC(v_g, quals, depths, GQ, DP, AB, an, ab_obs):
    """
    Apply indel genotype QC, masking genotypes which fail the QC
    Args:
        v_g (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        quals (list int): the allele number. v.gt_quals in cyvcf2
        depths (list int): list of depths for each individual. v.gt_depths in cyvcf2
        GQ (int): genotype quality threshold
        DP (int): depth threshold
        AB (float): allelic balance threshold
        an (int): allele number
        ab_obs (list float): the observed allelic balances
    Returns:
        tuple: (variant, pass rate) list of genotypes (format list [a1,a2,phased]) masked if failing thresholds, and the pass rate of the QC
    """
    failed=0
    for i in np.where((quals<GQ)|(depths<DP)|(ab_obs<AB))[0]:
        v_g[i]=[-1]*2 + [False]
        failed+=1

    return (v_g,(len(quals)-failed)/len(quals))


def indel_QC(v_g, quals, depths, combination, an, ab_obs):
    """
    Run indel genotype QC if variant QC passes
    Args:
        v_g (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        quals (list int): the allele number. v.gt_quals in cyvcf2
        depths (list int): list of depths for each individual. v.gt_depths in cyvcf2
        combination (list): list of QC thresholds [GQ, DP, AB, VQSLOD]
        an (int): allele number
        ab_obs (list float): the observed allelic balances
    Returns:
        list: list of genotypes (format list [a1,a2,phased]) masked if failing thresholds or return False if variant fails variant QC or too many genotypes fail QC
    """
#    if float(v.INFO.get('VQSLOD'))<combination[3]:
#        return False
#    else:
    v_g,pass_rate=apply_indel_genotype_QC(v_g, quals, depths, combination[0], combination[1], combination[2], an, ab_obs)
    if pass_rate < combination[4]:
        return False
    return v_g


def SNV_QC(v_g, vqslod, quals, depths, combination, an, pab_obs):
    """
    Run SNV genotype QC if variant QC passes
    Args:
        v_g (list list [int, int, bool]): list of genotypes, for each ind format = [a1, a2, phased]
        vqslod (float): VQSLOD score for variant
        quals (list int): the allele number. v.gt_quals in cyvcf2
        depths (list int): list of depths for each individual. v.gt_depths in cyvcf2
        combination (list): list of QC thresholds [GQ, DP, pAB, VQSLOD]
        an (int): allele number
        pab_obs (list float): the observed binomial of allelic depths
    Returns:
        list: list of genotypes (format list [a1,a2,phased]) masked if failing thresholds or return False if variant fails variant QC or too many genotypes fail QC
    """
    if vqslod==None:
        return False
    elif float(vqslod)<combination[3]:
        return False
    else:
        v_g,pass_rate=apply_SNV_genotype_QC(v_g, quals, depths, combination[0], combination[1], combination[2], an, pab_obs)
        if pass_rate < combination[4]:
            return False
    return v_g

def run_QC(variant, alt, combination_snv, combination_indel, an):
    """
    Run QC
    Args:
        variant (cyvcf2 variant): variant in the cyvcf2 format
        alt (str): alternate allele
        combination_snv (list): SNV QC thresholds = [ GQ, DP, P(AB), VQSLOD, FPASS ]
        combination_indel (list): Indel QC thresholds = [ GQ, DP, P(AB), VQSLOD, FPASS ]
        an (int): allele num
    Returns:
        list: list of genotypes (format list [a1,a2,phased]) masked if failing thresholds or return False if variant fails variant QC or too many genotypes fail QC
    """
    if (variant.is_snp | ((len(variant.REF)==1) & (len(alt)==1)) | (variant.REF[1:] == alt[1:])):
        pab_obs=calc_pab(variant.format('AD'), variant.genotypes, an+1, variant.gt_depths)
        variant_genotypes = SNV_QC(variant.genotypes, variant.INFO.get('VQSLOD'), variant.gt_quals, variant.gt_depths, combination_snv, an+1, pab_obs)
    elif variant.is_indel:
        ab_obs=calc_ab(variant.format('AD'), variant.genotypes, an+1)
        variant_genotypes = indel_QC(variant.genotypes, variant.gt_quals, variant.gt_depths, combination_indel, an+1, ab_obs)
    else:
        return False
    return variant_genotypes
