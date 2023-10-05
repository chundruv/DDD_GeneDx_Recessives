#!/bin/bash

python parse_vcfs.py --vcf chr1.bcf \
--pop test/population_table.txt \
--gen data/gene_pos.txt.gz \
--ped test/pedigree.ped \
--chr 1 --g1 1 --g2 10 --gmdversion gnomAD2.1 \
--unrelpar test/unrel_par.txt \
--snv-qc 20,7,1e-3,-2,0.5 \
--indel-qc 20,7,0.2,-1000,0.5 \
--clinpred ClinPred.txt.gz \
--cadd cadd_scores_hg19_v1.6/whole_genome_SNVs.tsv.gz \
--revel revel_scores.txt.gz \
--moipred MOI-Pred_hg19_liftover.txt.gz \
--cadd-indels cadd_scores_hg19_v1.6/InDels.tsv.gz \
--varity varity_all_predictions.txt.gz \
--spliceai spliceai_scores.masked.snv.hg19.vcf.gz \
--outdir test/output
