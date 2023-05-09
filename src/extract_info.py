from src.utils import get_minimal_representation
import tabix
import numpy as np

def init_tabix(args):
    """
    Initialise tabix instances
    Args:
        args (class): easier to pass in all args than the 8 seperate paths
    Returns:
        tuple (tabix instances): All the tabix instances
    """
    cadd_tabix = tabix.open(args.cadd)
    cadd_indels_tabix = tabix.open(args.cadd_indels)
    mpc_tabix = tabix.open(args.mpc)
    gnomad_tabix = tabix.open(args.gmd)
    revel_tabix = tabix.open(args.revel)
    primateai_tabix = tabix.open(args.primateai)
    varity_tabix = tabix.open(args.varity)
    spliceai_tabix = tabix.open(args.spliceai)
    spliceai_indels_tabix = tabix.open(args.spliceai_indels)
    return (cadd_tabix, cadd_indels_tabix, mpc_tabix, gnomad_tabix, revel_tabix, primateai_tabix, varity_tabix, spliceai_tabix, spliceai_indels_tabix)
    
def get_CADD(chrom, pos, ref, alt, compliments, cadd_tabix, cadd_indels_tabix):
    """
    Pull out the CADD score for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        is_indel (bool): True if an indel
        compliments (dict): Dictionary of compliment alleles
        cadd_tabix (tabix instance): cadd tabix instance
        cadd_indels_tabix (tabix instance): cadd indels tabix instance
    Returns:
        tuple (float, float): CADD raw and phred scores
    """

    pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
    size=abs(len(ref1)-len(alt1))
    if size>0:
        tabix = cadd_indels_tabix
    else:
        tabix = cadd_tabix

    for var in tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[2], var[3])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return var[5]
    return ''

def get_spliceai(chrom, pos, ref, alt, compliments, spliceai_tabix, spliceai_indels_tabix):
    pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
    size=abs(len(ref1)-len(alt1))
    if size>0:
        tabix = spliceai_indels_tabix
    else:
        tabix = spliceai_tabix
    for var in tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[3], var[4])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return max([float(i) for i in var[7].split('|')[2:6]])
    return ''

def get_REVEL(chrom, pos, ref, alt, compliments, revel_tabix):
    """
    Pull out the REVEL score for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        compliments (dict): Dictionary of compliment alleles
        revel_tabix (tabix instance): revel tabix instance
    Returns:
        float: REVEL score
    """
    ## pull out the REVEL score for the relevant position
    for var in revel_tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[3], var[4])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return var[7]
    return ''

def get_PrimateAI(chrom, pos, ref, alt, compliments, primateai_tabix):
    """
    Pull out the PrimateAI score for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        compliments (dict): Dictionary of compliment alleles
        primateai_tabix (tabix instance): PrimateAI tabix instance
    Returns:
        float: PrimateAI score
    """
    for var in primateai_tabix.query('chr'+str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[2], var[3])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return var[10]
    return ''

def get_VARITY(chrom, pos, ref, alt, compliments, varity_tabix):
    """
    Pull out the VARITY score for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        compliments (dict): Dictionary of compliment alleles
        varity_tabix (tabix instance): VARITY tabix instance
    Returns:
        tuple (float, float): VARITY_LOO  and VARITYER_LOO scores
    """
    for var in varity_tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[2], var[3])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return [var[10],var[11]]
    return ['','']

def get_MPC(chrom, pos, ref, alt, compliments,mpc_tabix):
    """
    Pull out the MPC score for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        compliments (dict): Dictionary of compliment alleles
        mpc_tabix (tabix instance): MPC tabix instance
    Returns:
        float: MPC score
    """
    for var in mpc_tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[2], var[3])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                return var[18]
    return ''

def get_gnomad_nhomalt(chrom, pos, ref, alt, compliments, gnomad_tabix):
    """
    Pull out the number of homalt genotypes in gnomAD for the relevant position
    Args:
        chrom (int): Chromosome
        pos (int): Base pair position
        ref (str): Reference allele
        alt (str): Alternate allele
        compliments (dict): Dictionary of compliment alleles
        gnomad_tabix (tabix instance): gnomAD tabix instance
    Returns:
        int: number of gnomAD homalts
    """
    for var in gnomad_tabix.query(str(chrom).replace("chr",""), pos-1, pos):
        pos1, ref1, alt1 = get_minimal_representation(pos, ref, alt)
        pos2, ref2, alt2 = get_minimal_representation(var[1], var[3], var[4])
        if pos1==pos2:
            if ((alt2==alt1) & (ref2==ref1)) or ((ref2==alt1) & (alt2==ref1)):
                nhomalt = int(var[7].split(';')[int(np.where(np.isin([i.split('=')[0] for i in var[7].split(';')],'nhomalt'))[0])].split('=')[1])
                #ac = int(var[7].split(';')[int(np.where(np.isin([i.split('=')[0] for i in var[7].split(';')],'AC'))[0])].split('=')[1])
                #if ac>0:
                return nhomalt#/ac
#                else:
#                    return 0
    return ''
