
# Note that this list of VEP annotations is current as of v93 with 2 included for backwards compatibility (VEP <= 75)
# From Konrad Karczewski
# Slight updates from Jack Kosmicki and Kaitlin Samocha
csq_order = ["transcript_ablation",
            "splice_donor_variant",
            "splice_acceptor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "initiator_codon_variant",  # deprecated
            "transcript_amplification",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "splice_region_variant",
            "incomplete_terminal_codon_variant",
            "start_retained_variant", # add on 2018-07-31
            "stop_retained_variant",
            "synonymous_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "non_coding_exon_variant",  # deprecated
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "nc_transcript_variant",  # deprecated
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "feature_elongation",
            "regulatory_region_variant",
            "feature_truncation",
            "intergenic_variant",
            ""]
csq_order_dict = dict(zip(csq_order, range(len(csq_order))))
rev_csq_order_dict = dict(zip(range(len(csq_order)), csq_order))

# list of coding consequences
coding_csq = ["transcript_ablation",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "initiator_codon_variant",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant"]

# categories for analysis
consequence_categories = {"lof":["transcript_ablation",
             "splice_donor_variant",
             "splice_acceptor_variant",
             "stop_gained",
             "frameshift_variant"],
    "missense/inframe":["stop_lost",
             "start_lost",
             "inframe_insertion",
             "inframe_deletion",
             "missense_variant"],
    "synonymous_variant":["synonymous_variant", "stop_retained_variant"]
}

def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst annotation (In this case, index of 'frameshift_variant', so 4)
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_index(csqs.split('&'))
    :param annnotation:
    :return most_severe_consequence_index:

    From Konrad Karczewski
    """
    return min([csq_order_dict[ann] for ann in csq_list])


def worst_csq_from_csq(csq):
    """
    Input possibly &-filled csq string (e.g. 'non_coding_exon_variant&nc_transcript_variant')
    Return the worst annotation (In this case, 'non_coding_exon_variant')
    :param consequence:
    :return most_severe_consequence:

    From Konrad Karczewski
    """
    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]

def compare_two_consequences(csq1, csq2):
    'From Konrad Karczewski'
    if csq_order_dict[worst_csq_from_csq(csq1)] < csq_order_dict[worst_csq_from_csq(csq2)]:
        return -1
    elif csq_order_dict[worst_csq_from_csq(csq1)] == csq_order_dict[worst_csq_from_csq(csq2)]:
        return 0
    return 1

def worst_csq_with_vep(annotation_list):
    """
    Takes list of VEP annotations [{'Consequence': 'frameshift', Feature: 'ENST'}, ...]
    Returns most severe annotation (as full VEP annotation [{'Consequence': 'frameshift', Feature: 'ENST'}])
    Also tacks on worst consequence for that annotation (i.e. worst_csq_from_csq)
    :param annotation_list:
    :return worst_annotation:

    From Konrad Karczewski

    Modified on 2019-12-11 by Kaitlin Samocha to prioritize protein coding genes
    """
    if len(annotation_list) == 0: return None
    worst = annotation_list[0]
    for annotation in annotation_list:
        result_comparison = compare_two_consequences(annotation['Consequence'], worst['Consequence'])

        if result_comparison < 0:
            worst = annotation
        elif result_comparison == 0 and annotation['CANONICAL'] == 'YES':
            worst = annotation
        #elif result_comparison == 0 and annotation['BIOTYPE'] == 'protein_coding':
        #    worst = annotation

    worst['major_consequence'] = worst_csq_from_csq(worst['Consequence'])
    return worst


def get_canonical_and_worst_csq(variant, csq, an, gene, consequence_categories=consequence_categories):
    annotations = [dict(zip(csq, i.split('|'))) for i in variant.INFO.get('CSQ').split(',') if len(csq) == len(i.split('|')) and (int(i.split('|')[18]) == an+1) and (i.split('|')[csq.index('BIOTYPE')] == 'protein_coding') and (i.split('|')[csq.index('Gene')] == gene)]
    
    worst_vep_csq=worst_csq_with_vep(annotations)
    vep_csq=worst_csq_with_vep([i for i in annotations if i['CANONICAL']=='YES'])
    
#    if vep_csq == []:
#        return ([], [], [], [])

    if vep_csq!=None:

        if len(vep_csq['Consequence'].split('&')) >1:
            wc=0
            for i in range(1, len(vep_csq['Consequence'].split('&'))):
                comp=compare_two_consequences(vep_csq['Consequence'].split('&')[wc], vep_csq['Consequence'].split('&')[i])
                if comp==1:
                    wc=i
                else:
                    continue
            vep_csq['Consequence']=vep_csq['Consequence'].split('&')[wc]

        vep_csq_cat=[i if vep_csq['Consequence'] in consequence_categories[i] else '' for i in consequence_categories][0]
    else:
        vep_csq_cat=''

    if worst_vep_csq!=None:
        worst_csq_cat=[i if worst_vep_csq['Consequence'] in consequence_categories[i] else '' for i in consequence_categories][0]
    else:
        worst_csq_cat=''

    return (vep_csq, vep_csq_cat, worst_vep_csq, worst_csq_cat)
