<<<<<<< HEAD
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

=======
python parse.py --vcf 1.bcf 
--pop pops.txt 
--gen gene_pos_hg19_protein_coding_only_unique_EnsemblIDs_intersect_baits.txt 
--ped sanger_and_ega_ids.ped 
--cadd whole_genome_SNVs.tsv.gz 
--revel revel_scores.txt.gz 
--chr 1 --g1 1 --g2 10 --gmdversion gnomAD2.1 
--unrelpar unrelated_unaffected_parents_KINGgt0.04419417.txt 
--cadd-indels cadd_scores_hg19_v1.6/InDels.tsv.gz 
--varity varity_all_predictions.txt.gz 
--spliceai spliceai_scores.masked.snv.hg19.vcf.gz 
--outdir outdir
>>>>>>> f438d831a5e1f836d36d328b9fc0ae798fc6ecd4
