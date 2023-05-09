def get_minimal_representation(pos, ref, alt):
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets,
    with different combinations of alternate alleles, can always be matched directly.
    Note that chromosome is ignored here - in xbrowse, we'll probably be dealing with 1D coordinates
    Args:
        pos (int): genomic position in a chromosome (1-based)
        ref (str): ref allele string
        alt (str): alt allele string
    Returns:
        tuple: (pos, ref, alt) of remapped coordinate

    Modified based on the original from Eric Minkel:
    https://github.com/ericminikel/minimal_representation/blob/master/minimal_representation.py
    """

    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1:
        return (pos, ref, alt)
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return (pos, ref, alt)

def extract_csq(x):
    for i in x.header_iter():
        if i.type=='INFO':
            if i['ID']=='CSQ':
                csq=i['Description'][1:-1].split(': ')[1].split('|')
    return(csq)