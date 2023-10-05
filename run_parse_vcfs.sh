#!/bin/bash

python parse_vcfs.py --vcf ~/vcfs_ddd/1.bcf \
--pop test/population_table.txt \
--gen data/gene_pos.txt.gz \
--ped test/pedigree.ped \
--cadd /nfs/users/nfs_k/kc18/resources/cadd_scores/cadd_scores_hg19_v1.6/whole_genome_SNVs.tsv.gz \
--revel /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/revel_scores.txt.gz \
--chr 1 --g1 1 --g2 10 --gmdversion gnomAD2.1 \
--unrelpar test/unrel_par.txt \
--snv-qc 20,7,1e-3,-2,0.5 \
--indel-qc 20,7,0.2,-1000,0.5 \
--clinpred /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/ClinPred.txt.gz \
--moipred /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/moipred/MOI-Pred_hg19_liftover.txt.gz \
--cadd-indels /nfs/users/nfs_k/kc18/resources/cadd_scores/cadd_scores_hg19_v1.6/InDels.tsv.gz \
--varity /nfs/users/nfs_k/kc18/ddd/recwest/data/missense_pathogenicity_scores/varity/varity_all_predictions.txt.gz \
--spliceai /lustre/scratch125/humgen/resources/SpliceAI_data_files/spliceai_scores.masked.snv.hg19.vcf.gz \
--outdir test/output
