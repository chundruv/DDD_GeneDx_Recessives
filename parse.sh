#!/bin/bash

python parse.py --vcf test/test1.vcf.gz \
--pop test/population_table.txt \
--gen data/gene_pos.txt.gz \
--ped test/pedigree.ped \
--cadd /nfs/users/nfs_k/kc18/resources/cadd_scores/cadd_scores_hg19_v1.6/whole_genome_SNVs.tsv.gz \
--revel /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/revel_scores.txt.gz \
--chr 1 --g1 1 --g2 10 --gmdversion gnomAD2.1 \
--unrelpar unrel_par.txt \
--cadd-indels /nfs/users/nfs_k/kc18/resources/cadd_scores/cadd_scores_hg19_v1.6/InDels.tsv.gz \
--varity /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/varity/varity_all_predictions.txt.gz \
--spliceai /lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.masked.snv.hg19.vcf.gz \
--outdir test/output
