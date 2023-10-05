# Recessive burden testing

This repository contains code to perform recessive burden testing as seen in Chundru _et al_ 2023.

## Setup

The code requires a number of python packages. The easiest way to install all is to use the recessives_environment.yml:

`conda create -f recessives_environment.yml`

The other required files are:
* CADD SNV and indel files (https://cadd.gs.washington.edu/download)
* SpliceAI file (https://basespace.illumina.com/analyses/194103939/files?projectId=66029966 - NOTE: only need to download the one ~27G file not the whole >300GB folder - filename: spliceai_scores.masked.snv.hg19.vcf.gz)
* REVEL file (https://sites.google.com/site/revelgenomics/downloads?authuser=0)
* VARITY file (http://varity.varianteffect.org/downloads/varity_all_predictions.tar.gz)
* ClinPred file (https://sites.google.com/site/clinpred/download?authuser=0)
* MOIPred file (https://doi.org/10.5281/zenodo.5620519)

If you don't want to annotate with all of these then you can tweak the code to remove them, or if they are all annotated in the VEP string you can remove the annotation steo and get them from the VEP

## Running analysis
The code is set up to run in three steps.

### Step 1
The first step involves parsing the vcfs, extracting rare heterozygous and homozygous alternative genotypes in parents and probands, and annotating these genotypes with various metrics.

You will need:
* indexed VCF files (.vcf.gz or .bcf)
* pedigree file (tab delimited with header: "family_id	individual_id	dad_id	mum_id	sex	affected")
* population table (tab delimited no header, columns - ID, GIA-group, GIA-subgroup)
* list of unrelated, unaffected parents
* QC metrics to use (If genotypes not passing QC are already masked and variants removed, you can put 0,0,0,-1000,0 for both QC thresholds)
* all the metrics listed above (CADD, REVEL, etc.)

Edit the paths within and run `./parse_vcfs.sh`

The output should have the header (or similar):

`stable_id       genotype        is_proband      child_id        child_genotype  dad_id  dad_genotype    mum_id  mum_genotype    child_inheritance       parent_inheritance      n_unaffected_homalt     variant_class   variant_id      worst_vep_annotation_category   worst_vep_annotation    canonical_vep_annotation_category       canonical_vep_annotation        gene.stable.id  gene.name       HGNC.symbol     chrom   position        ref     alt     EUR1_unrelated_parents_popents_af       SAS1_unrelated_parents_popents_af       gnomAD2.1_AFR   gnomAD2.1_AMR   gnomAD2.1_EAS   gnomAD2.1_FIN   gnomAD2.1_NFE   gnomAD2.1_SAS   gnomAD2.1_ASJ   population      n_unrelated_parents_EUR1        n_unrelated_parents_SAS1        CADD_phred      LOFTEE  REVEL   VARITYR_LOO     VARITYER_LOO    PolyPhen        spliceai        clinpred        moipred_rp`

### Step 2

The second step applies filters the variants and calls compound heterozygous variants. Using these and the ROHs called in the sample we run a burden test for the following variant classes:
LoF/LoF - loss of function biallelic
LoF/Functional - loss of function - functional compound heterzygous variants
Functional/Functional - functional biallelic variants
Synonymous/Synonymous - synonymous biallelic variants

You will need:
* output from Step 1
* population table (same as from Step 1)
* pedigree 
* list of unrelated, unaffected parents
* list of unrelated probands
* ROHs output from bcftools-roh (No header, columns: "RG", Sample, Chromosome, Start, End, Length(bp), Number of markers, Quality)

To run, edit the paths within and run `./burden_test.sh`

The output will be in two folders, parts/ and vars/
parts/ will have the observed and expected for each gene in the chunk of the genome analysed, each population, LD R2 pruning for ROH calling, and variant class
vars/ has the observed biallelic variants

### Step 3

This aggregates the results from all the chunks and calculates the exome-wide burden and the per-gene enrichment

run `Rscript results.R $OUTDIR` where OUTDIR is the directory containing the output folders from step2

If you have any issues/questions please get in touch via the issues tab on this repository or v.chundru [at] exeter.ac.uk
