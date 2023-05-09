import pandas as pd
import numpy as np
import sys
import csv
import pyranges as pr
import argparse
from src.read_and_filter_data import read_file
from src.call_comphets import call_comphets
from src.genes_positions import get_gene_pos
from src.calc_expected import calc_expected
from src.calc_observed import calc_obs

parser = argparse.ArgumentParser( description='Run burden test' )
parseargs = parser.add_argument_group('General inputs')
parseargs.add_argument(
        '--chrom',
        type = str,
        required = True,
        help = 'Chromosome')
parseargs.add_argument(
        '--g1',
        type = int,
        required = False,
        default = 0,
        help = 'Segment start index')
parseargs.add_argument(
        '--g2',
        type = int,
        required = False,
        default = 17320,
        help = 'Segment start index')
parseargs.add_argument(
        '--input_dir',
        type = str,
        required = True,
        help = 'Input directory')
parseargs.add_argument(
        '--output_dir',
        type = str,
        required = True,
        default = "out/",
        help = 'Output directory')
parseargs.add_argument(
        '--genepos',
        type = str,
        required = True,
        default = "Data/genepos.txt",
        help = 'Gene positions file')
parseargs.add_argument(
        '--rohs_dir',
        type = str,
        required = True,
        help = 'Directory containing bcftools-roh calls')
parseargs.add_argument(
        '--idmap',
        type = str,
        required = False,
        help = 'ID map file')
parseargs.add_argument(
        '--popfile',
        type = str,
        required = True,
        help = 'Population assignment file')
parseargs.add_argument(
        '--pedfile',
        type = str,
        required = True,
        help = 'Pedigree file')
parseargs.add_argument(
        '--unrel_probands',
        type = str,
        required = True,
        help = 'File of unrelated proband IDs')
parseargs.add_argument(
        '--unrel_parents',
        type = str,
        required = True,
        help = 'File of unrelated parent IDs')
parseargs.add_argument(
        '--unaff_parents',
        type = str,
        required = True,
        help = 'File of unaffected parent IDs')
parseargs.add_argument(
        '--qcfail',
        type = str,
        required = True,
        help = 'File of individual IDs who fail QC')
parseargs.add_argument(
        '--skipGenes',
        type = str,
        required = False,
        help = 'List of genes to skip')

missense_filters = parser.add_argument_group('Change missense filtering thresholds')
missense_filters.add_argument(
        '--cadd_threshold',
        type = float,
        required = False,
	default = 24.18,
        help = 'CADD threshold')
missense_filters.add_argument(
        '--cadd_indel_threshold',
        type = float,
        required = False,
        default = 17.34,
        help = 'CADD indel threshold')
missense_filters.add_argument(
        '--revel_threshold',
        type = float,
        required = False,
        default = 0.36,
        help = 'REVEL threshold')
missense_filters.add_argument(
        '--varity_threshold',
        type = float,
        required = False,
        default = 0.25,
        help = 'VARITY-ER threshold')
missense_filters.add_argument(
        '--polyphen_threshold',
        type = float,
        required = False,
        default = 0.59,
        help = 'PolyPhen threshold')
missense_filters.add_argument(
        '--clinpred_threshold',
        type = float,
        required = False,
        default = 0.53,
        help = 'ClinPred threshold')
missense_filters.add_argument(
        '--moipred_threshold',
        type = float,
        required = False,
        default = 0.11,
        help = 'MOIpred recessive probability threshold')
missense_filters.add_argument(
        '--synsplice_threshold',
        type = float,
        required = False,
        default = 0.8,
        help = 'MOIpred recessive probability threshold')
args = parser.parse_args()

x, parents_populations, probands_populations, populations, N_haps, N_probands, unrel_parents, unrel_probands = read_file(args)
gene_pos, not_used = get_gene_pos(args.gen_path, args.chrom, args.g1, args.g2)

chets = call_comphets(x)

if args.skipGenes==None:
    genes=x['gene.stable.id'].unique()
else:
    skipgenes=pd.read_csv(args.skipGenes, header=None)
    genes=[g for g in x['gene.stable.id'].unique() if g not in skipgenes[0].tolist()]

classes=['lof', 'missense/inframe', 'synonymous_variant']
lds=[0.2, 0.4, 0.6, 0.8]

lmbda, lmbda1 = calc_expected(x, lds, genes, classes, parents_populations, populations, N_haps, a, consequence_classes, unrel_parents)

consequence_classes={'lof/lof':['lof'], 'lof/missense':['lof', 'missense/inframe'], 'missense/missense':['missense/inframe'], 'synonymous/synonymous':['synonymous_variant']}
OB, varIDs = calc_obs(x, unrel_parents, unrel_probands, populations, consequence_classes, genes, probands_populations, chets)

with open(args.output_dir+'/parts/chr'+str(args.chrom)+'_'+str(args.g1)+'_'+str(args.g2)+'.txt', 'w') as write_file:
    csvwrite=csv.writer(write_file, delimiter='\t')
    csvwrite.writerow(['Gene', 'Variant_class', 'Population', 'LD_thinning_r2', 'Observed_biallelic_genotypes', 'N_probands', 'N_parents', 'Expected_freq_biallelic_genotypes', 'corrected_Expected_freq_biallelic_genotypes'])
    for pop in populations:
        for c in consequence_classes:
            for gene in genes:
                for ld in lds:
                    csvwrite.writerow([gene, c, pop, ld, OB[c][gene][pop], N_probands[pop], N_haps[pop]/2, lmbda[c][gene][pop][ld], lmbda1[c][gene][pop][ld]])

with open(args.output_dir+'/vars/chr'+str(args.chrom)+'_'+str(args.g1)+'_'+str(args.g2)+'.txt', 'w') as write_file:
    csvwrite=csv.writer(write_file, delimiter='\t')
    csvwrite.writerow(['Gene', 'Variant_class', 'Population', 'Individual_ID', 'Variant_ID', 'Genotype'])
    for c in consequence_classes:
        for gene in genes:
            for pop in populations:
                for sample in varIDs[c][gene][pop]:
                    for varid in varIDs[c][gene][pop][sample]:
                        csvwrite.writerow([gene, c, pop, sample, varid, varIDs[c][gene][pop][sample][varid]])


