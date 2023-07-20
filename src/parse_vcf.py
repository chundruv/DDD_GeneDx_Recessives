from src.inheritance import get_inh
from src.process import *
from src.genes_positions import get_gene_pos
from src.find_hethoms import find_hethoms
from src.apply_qc import run_QC
from src.extract_info import *
from src.utils import extract_csq
from src.vep_funcs import get_canonical_and_worst_csq
from cyvcf2 import VCF
import pandas as pd
import numpy as np
from scipy import stats
import csv
import os

def variant_parser(args):
    """

    Parse through the VCF files using cyvcf2 (https://brentp.github.io/cyvcf2/docstrings.html)

    Input args:
        vcf =   The input vcf file. Either vcf.gz or bcf; index is required
        g2p =   DDG2P genes file downloaded from DECIPHER (https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz)
        pop =   Populations assignment file. Columns: "ID1 ID2 Pop" - For DDD ID1=stable ID, ID2=VCF ID. For GeneDx put same ID in both
        unrelpar    =   unrelated unaffected parents list. Columns "ID1 ID2" (ID1 and ID2 same as above) 
        gen =   gene position file downloaded from biomart; required columns: 'HGNC symbol', 'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)'
        ped =   family PED file; required columns: Family ID, Individual ID, Dad ID, Mum ID, Sex, Affected status
        cadd    =   SNV CADD file (https://cadd.gs.washington.edu/download)
        cadd-indels =   Indels CADD file (https://cadd.gs.washington.edu/download)
        revel   =   REVEL file bgzipped and tabix indexed (https://sites.google.com/site/revelgenomics/downloads?authuser=0)
        mpc =   Variant MPC file (ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/)
        primateai   =   PrimateAI file bgzipped and tabix indexed (https://basespace.illumina.com/s/cPgCSmecvhb4)
        varity  =   VARITY file bgzipped and tabix indexed (http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz)
        gmd =   gnomAD exomes v2 allele frequency file (https://gnomad.broadinstitute.org/downloads)
        gmd_genomes = gnomAD genomes v2 allele frequency file (https://gnomad.broadinstitute.org/downloads)
        mcg =   gene missense constraint file (not sure where this is from, I got it from Joanna Kaplanis)
        mcr =   missesne constraint regions file (https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/regional_missense_constraint.tsv)
        pli =   pLI scores (https://storage.googleapis.com/gcp-public-data--gnomad/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz)
        snv_qc = QC thresholds for SNVs format= GQ,DP,P(AB),VQSLOD,FPASS
        indel_qc = QC thresholds for indels format= GQ,DP,AB,VQSLOD,FPASS
        chr =   chromosome
        g1  =   start gene index
        g2  =   end gene index
        gmdversion  =   gnomAD AF prefix in VEP string
        outdir  =   output directory

    Writes output to file "${outdir}/chr${chrom}_${g1}_${g2}_recessive_candidate.txt"

    """

    chrom=args.chr
    g1=args.g1
    g2=args.g2
    
    x=VCF(args.vcf, gts012=True)
    
    samples=x.samples
    genes, gene_symbols = get_gene_pos(args.gen, chrom, g1, g2)
    ped, parents, unaffected, samples_pop, samples_pop_inds, unrelated_parents_pop_index, populations = process_ped(args.pop, args.unrelpar, args.ped, samples)
    csq = extract_csq(x)

    compliments = {"A":"T","C":"G","G":"C","T":"A"}

    gene_csqpos = csq.index('Gene')

    cadd_tabix, cadd_indels_tabix, mpc_tabix, gnomad_tabix, revel_tabix, primateai_tabix, varity_tabix, spliceai_tabix, spliceai_indels_tabix = init_tabix(args)

    # Fomart = [ GQ, DP, P(AB), VQSLOD, FPASS ]
    # SNV QC - DDD GQ>20, DP>7, P(AB)>1e-3, VQSLOD>-2, FPASS>0.5
    # SNV QC - GeneDx GQ>25, DP>10, P(AB)>1e-3, VQSLOD>-2.5, FPASS>0.7 
    # combination_snv = [20,7,1e-3,-2,0.5]
    combination_snv = args.snv_qc.split(',')
    # INDEL QC - DDD GQ>20, DP>7, AB>0.2, No VQSLOD, FPASS>0.5
    # INDEL QC - GeneDx GQ>30, DP>7, AB>0.3, No VQSLOD, FPASS>0.7
    # combination_indel = [20,7,0.2,-1000,0.5]
    combination_indel = args.indel_qc.split(',')

    csv_header=['stable_id','genotype','is_proband','child_id','child_genotype','dad_id','dad_genotype','mum_id','mum_genotype','child_inheritance','parent_inheritance','n_unaffected_homalt','variant_class',
            'variant_id','worst_vep_annotation_category','worst_vep_annotation','canonical_vep_annotation_category','canonical_vep_annotation','gene.stable.id', 'gene.name', 'HGNC.symbol','chrom','position','ref',
            'alt']+[pop+'_unrelated_parents_popents_af' for pop in populations]+['gnomAD2.1_AFR','gnomAD2.1_AMR','gnomAD2.1_EAS','gnomAD2.1_FIN','gnomAD2.1_NFE','gnomAD2.1_SAS','gnomAD2.1_ASJ',
            'gnomad2.1_nhomalt','population']+['n_unrelated_parents_'+pop for pop in populations]+['CADD_phred','LOFTEE','MPC','REVEL','PrimateAI','VARITYR_LOO','VARITYER_LOO','PolyPhen', 'spliceai']

    if os.path.exists(args.outdir)==False:
	    os.makedirs(args.outdir)
    
    with open(args.outdir+'/chr'+chrom+'_'+g1+'_'+g2+'_recessive_candidate.txt', 'w') as write_file:
        csvwrite=csv.writer(write_file, delimiter='\t')
        csvwrite.writerow(csv_header)
        
        for gene in genes.keys():

            for variant in x(genes[gene]):

                for alt in variant.ALT:

                    if alt !="*":

                        an = variant.ALT.index(alt)

                        if np.all(np.array([i.split('|')[gene_csqpos] for i in variant.INFO.get('CSQ').split(',')])!=gene):
                            continue

                        variant_genotypes = run_QC(variant, alt, combination_snv, combination_indel, an)

                        if variant_genotypes==False:
                            continue

                        (vcf_AF, gt_types) = store_gts(variant, variant_genotypes, an)
                
                        af_pop = calc_af_pop(populations, gt_types, unrelated_parents_pop_index)

                        vep_csq, vep_csq_cat, worst_vep_csq, worst_csq_cat = get_canonical_and_worst_csq(variant, csq, an, gene)

                        if worst_vep_csq==None:
                            continue
                        if len(worst_vep_csq)==0:
                            continue

                        if vep_csq==None:
                            vep_csq={'Consequence':'', 'LoF':'', 'PolyPhen':''}

                        gnomad_af=[worst_vep_csq[args.gmdversion+'_AF_afr'],
                        worst_vep_csq[args.gmdversion+'_AF_amr'],
                        worst_vep_csq[args.gmdversion+'_AF_eas'],
                        worst_vep_csq[args.gmdversion+'_AF_fin'],   
                        worst_vep_csq[args.gmdversion+'_AF_nfe'],
                        worst_vep_csq[args.gmdversion+'_AF_sas'],   
                        worst_vep_csq[args.gmdversion+'_AF_asj']]
                
                        if vcf_AF[an]>0.5:
                            gnomad_af = [1-float(i) if i!='' else '' for i in gnomad_af]
                
                        #### Here we parse through the individuals finding the heterozygote and homozygote alternative genotypes which pass the AF thresholds
                        hethoms = find_hethoms(af_pop, gnomad_af, samples, samples_pop, gt_types)
                                
                        if len(hethoms)>0:
                            cadd_phred = get_CADD(chrom, variant.POS, variant.REF, alt, compliments, cadd_tabix, cadd_indels_tabix)
                            spliceai = get_spliceai(chrom, variant.POS, variant.REF, alt, compliments, spliceai_tabix, spliceai_indels_tabix)
                            mpc = get_MPC(chrom, variant.POS, variant.REF, alt, compliments, mpc_tabix)
                            revel = get_REVEL(chrom, variant.POS, variant.REF, alt, compliments, revel_tabix)
                            primateai = get_PrimateAI(chrom, variant.POS, variant.REF, alt, compliments, primateai_tabix)
                            varity = get_VARITY(chrom, variant.POS, variant.REF, alt, compliments, varity_tabix)
                            gnomad_nhomalt = get_gnomad_nhomalt(chrom, variant.POS, variant.REF, alt, compliments, gnomad_tabix)

                            for i in hethoms:
                                pop, child_inh, parent_inh = get_inh(samples_pop_inds, parents, ped, samples, gt_types, variant, i)

                                if pop!='':            
                                    csvwrite.writerow(['NA' if x=='' else x for x in [i, [var if var!=3 else 'NA' for var in gt_types[np.where(np.isin(samples,i))]][0], i in ped,
                                        i if i in parents else 'NA',
                                        [var if var!=3 else 'NA' for var in gt_types[np.where(np.isin(samples,parents[i]))]][0] if i in parents else 'NA',
                                        i if i in ped else 'NA',[var if var!=3 else 'NA' for var in gt_types[np.where(np.isin(samples,ped[i][0]))]][0] if i in ped else 'NA',
                                        i if i in ped else 'NA',[var if var!=3 else 'NA' for var in gt_types[np.where(np.isin(samples,ped[i][1]))]][0] if i in ped else 'NA',
                                        child_inh, parent_inh, sum(gt_types[np.where(np.isin(samples,unaffected))[0]]==2),worst_vep_csq['VARIANT_CLASS'],
                                        str(variant.CHROM)+':'+str(variant.POS)+'_'+variant.REF+'_'+alt,worst_csq_cat, worst_vep_csq['Consequence'] ,vep_csq_cat,vep_csq['Consequence'], gene, gene_symbols[gene][0], gene_symbols[gene][1],variant.CHROM, 
                                        variant.POS,variant.REF, alt] + list(af_pop.values()) + gnomad_af +[gnomad_nhomalt, pop]+[sum(gt_types[unrelated_parents_pop_index[pop]]!=3) for pop in populations]+ 
                                        [cadd_phred, vep_csq['LoF'], mpc, revel, primateai,varity[0],varity[1],vep_csq['PolyPhen'], spliceai]])
