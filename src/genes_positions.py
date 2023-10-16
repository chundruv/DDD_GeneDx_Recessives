import pandas as pd
import numpy as np
import csv

def get_gene_pos(gen_path, chrom, g1, g2):
    """
    Get gene start and end positions
    Args:
        gen_path (str): path to gene positions file
        chrom (int): Chromosome
        g1 (int): gene start index
        g2 (int): gene end index
    Returns:
        dict: dictionary of chromosome (key), with value of dictionary with keys of gene symbol and values of gene postion (format = chrom:start_bp-end_bp)
    """
    gen=pd.read_csv(gen_path, sep='\t')

    gen=gen[(gen['Chromosome/scaffold name']==int(chrom))]
    
    gen.loc[:,'pos']=gen['Chromosome/scaffold name'].astype(str)+':'+gen['Gene start (bp)'].astype(str)+'-'+gen['Gene end (bp)'].astype(str)

    genes = gen['Gene stable ID'].unique()[int(g1):int(g2)+1]

    return [{i:gen[gen['Gene stable ID']==i]['pos'].tolist()[0] for i in genes}, {i:[gen[gen['Gene stable ID']==i]['Gene name'].tolist()[0], gen[gen['Gene stable ID']==i]['HGNC symbol'].tolist()[0]] for i in genes}]
