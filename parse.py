from src.parse_vcf import variant_parser
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser( description='Extract variant info from the VCF and other source files' )
    required = parser.add_argument_group('required arguments')
    required.add_argument(
            '--vcf',
            type = str,
            required = True,
            help = 'The input vcf file.')
    required.add_argument(
            '--pop',
            type = str,
            required = True,
            help = 'Populations assignment file')
    required.add_argument(
            '--unrelpar',
            type = str,
            required = True,
            help = 'unrelated parents list')
    required.add_argument(
            '--gen',
            type = str,
            required = True,
            help = 'gene position file')
    required.add_argument(
            '--ped',
            type = str,
            required = True,
            help = 'family PED file')
    required.add_argument(
            '--cadd',
            type = str,
            required = True,
            help = 'SNV CADD file')
    required.add_argument(
            '--cadd-indels',
            type = str,
            required = True,
            help = 'Indels CADD file')
    required.add_argument(
            '--spliceai',
            type = str,
            required = True,
            help = 'SNV spliceai file')
    required.add_argument(
            '--revel',
            type = str,
            required = True,
            help = 'REVEL file bgzipped and tabix indexed')
    required.add_argument(
            '--varity',
            type = str,
            required = True,
            help = 'VARITY file bgzipped and tabix indexed')
    required.add_argument(
            '--chr',
            type = str,
            required = True,
            help = 'chromosome number')
    required.add_argument(
            '--snv-qc',
            type = str,
            required = True,
            help = 'SNV QC thresholds, comma-seperated list - GQ,DP,P(AB),VQSLOD,FPASS ')
    required.add_argument(
            '--indel-qc',
            type = str,
            required = True,
            help = 'Indel QC thresholds, comma-seperated list - GQ,DP,AB,VQSLOD,FPASS ')
    required.add_argument(
            '--g1',
            type = str,
            required = True,
            help = 'starting gene number')
    required.add_argument(
            '--g2',
            type = str,
            required = True,
            help = 'ending gene number')
    required.add_argument(
            '--gmdversion',
            type = str,
            required = True,
            help = 'Version of gnomAD in VEP')
    required.add_argument(
            '--outdir',
            type = str,
            required = True,
            help = 'output directory')
    args = parser.parse_args()
    variant_parser(args)
