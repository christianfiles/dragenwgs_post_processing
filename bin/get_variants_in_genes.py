#!/usr/bin/env python

import argparse
import csv
import pandas as pd
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet

"""
Script for extracting variants in certain genes from a vep annotated VCF.

Must have the following annotations within the CSQ - SYMBOL, PICK, Existing_variation
Must have the following annotations within the INFO - gnomad_popmax_af

Usage:

python bin/get_variants_in_genes.py \
--vcf results/annotated_vcf/170714_D00501_0145_BHMM5GBCXY_syn8.norm.anno.vcf.gz \
--proband_id 16M12401 \
--ped results/ped/170714_D00501_0145_BHMM5GBCXY_syn8.ped \
--gene_list gene_list.txt \
--output test.csv \
--min_dp 10 \
--min_gq 10 \
--min_qual_snps 13 \
--min_qual_indels 13 \
--comment test

"""


parser = argparse.ArgumentParser(description='Find UPD events in NGS Trio Data')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='The path to the annotated VCF file. Must be bgzipped and tabixed.')
parser.add_argument('--proband_id', type=str, nargs=1, required=True,
				help='The Sample ID of the proband in the VCF.')
parser.add_argument('--ped', type=str, nargs=1, required=True,
				help='A ped file describing the family relationships.')
parser.add_argument('--gene_list', type=str, nargs=1, required=True,
				help='A text file with required gene names (symbols) e.g. BRCA1.')
parser.add_argument('--output', type=str, nargs=1, required=True, help='The output name prefix.')
parser.add_argument('--min_dp', type=int, nargs=1, required=True, help='The minimum genotype depth.')
parser.add_argument('--min_gq', type=int, nargs=1, required=True, help='The minimum genotype quality (GQ).')
parser.add_argument('--min_qual_snps', type=float, nargs=1, required=True, help='The minimum QUAL value for snps.')
parser.add_argument('--min_qual_indels', type=float, nargs=1, required=True, help='The minimum QUAL value for indels.')
parser.add_argument('--comment', type=str, nargs=1, required=True, help='Comment to add to report e.g. 20-1234')
args = parser.parse_args()


vcf = args.vcf[0]
proband_id = args.proband_id[0]
ped = args.ped[0]
min_dp = args.min_dp[0]
min_gq = args.min_gq[0]
min_qual_snps = args.min_qual_snps[0]
min_qual_indels = args.min_qual_indels[0]
gene_list = args.gene_list[0]
comment = args.comment[0]
output = args.output[0]

# define a filtering function

def passes_initial_filter(variant, proband_id, gene_dict, min_gq, min_qual_snps, min_qual_indels, min_dp):
	"""
	An initial filter to 

	We import if the variant passes quality filtering and is below 1% in gnomad exomes and gnomad genomes AND

	a) Is listed as pathogenic in clinvar OR
	b) Has a a relevant consequence
	OR is in the whitelist

	"""

	# If the proband has the variant and we pass the genotype and variant level filters
	if variant.has_alt(proband_id) and variant.passes_gt_filter(proband_id, min_gq=min_gq, min_dp=min_dp):

		# ignore mt variants for now
		if variant.chrom == 'MT':

			return False

		else:

			#get rid of poor quality variants
			if variant.is_snp():

				if variant.quality < min_qual_snps:

					return False

			else:

				if variant.quality < min_qual_indels:

					return False

		genes = variant.get_genes(feature_key='SYMBOL')

		in_panel = False

		for gene in genes:

			if gene in gene_dict:

				return True

	return False

def fix_genotype(df, column_key):
	"""
	Convert genotype in format G/A to HET etc
	"""

	genotype = df[column_key]
	variant  = df['variant_id']

	alt = variant.split('>')[1]

	if '/' in genotype:

		genotype = genotype.split('/')

	elif '|' in genotype:

		genotype = genotype.split('|')

	else:

		raise Exception('weird genotype')

	if genotype.count('.') == 2:

		return 'MISSING'

	elif genotype.count(alt) == 2:

		return 'HOM_ALT'

	elif genotype.count(alt) == 1:

		return 'HET'

	elif genotype.count(alt) == 0:

		return 'HOM_REF'

	else:

		return 'UNKNOWN'

# read ped
ped_df = pd.read_csv(ped, sep='\t', names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'affected'])

# get family id
filtered_ped = ped_df[ped_df['sample_id']==proband_id]
family_id = filtered_ped['family_id'].iloc[0]

# create family object
my_family = Family(family_id)	
my_family.read_from_ped_file(ped, family_id, proband_id)

# read gene file and convert to dict

gene_dict = {}

with open(gene_list) as csvfile:

	spamreader = csv.reader(csvfile, delimiter='\t')

	for row in spamreader:

		row_gene = row[0].strip()
		gene_dict[row_gene] = row_gene

print ('Extracting variants in the following genes.')
print (gene_dict)

# Create a new VariantSet object
my_variant_set = VariantSet()

# Associate the my_family object with my_variant_set
my_variant_set.add_family(my_family)

my_variant_set.read_variants_from_vcf(vcf,
									proband_variants_only=True,
									filter_func=passes_initial_filter,
									args=(proband_id, gene_dict,min_gq, min_qual_snps, min_qual_indels, min_dp))

print (f'{len(my_variant_set.variant_dict)} variants have been loaded into the variant set.')


variant_df = my_variant_set.to_df(min_parental_gq_dn= min_gq, min_parental_depth_dn=min_dp, min_parental_gq_upi=min_gq, min_parental_depth_upi= min_dp)

gt_fields= []

# convert genotype fields
for fm in my_family.get_all_family_member_ids():

	variant_df[f'{fm}_Genotype'] = variant_df.apply(fix_genotype, axis=1, args=( f'{fm}_GT',))

proband_id = my_family.get_proband_id()
variant_df['#SampleId'] = proband_id
variant_df['WorklistId'] = comment
variant_df['Variant'] = variant_df['variant_id']
variant_df['Genotype'] = variant_df[f'{proband_id}_Genotype']
variant_df['SYMBOL'] = variant_df['csq_SYMBOL']
variant_df['Feature'] = variant_df['csq_Feature']
variant_df['Consequence'] = variant_df['csq_Consequence']
variant_df['HGVSc'] = variant_df['csq_HGVSc']
variant_df['HGVSp'] = variant_df['csq_HGVSp']
variant_df['CLIN_SIG'] = variant_df['csq_CLIN_SIG']
variant_df['Existing_variation'] = variant_df['csq_Existing_variation']
variant_df['AutoPick'] = variant_df['csq_PICK']
variant_df['gnomad_popmax_af'] = variant_df['info_gnomad_popmax_af']

# add sample specific columns such as genotype and depth
for fm in my_family.get_all_family_member_ids():

	# ignore proband as we did it earlier
	if fm != proband_id:

		for field in ['_Genotype', '_DP', '_GQ', '_AD']:

			gt_fields.append(fm + field)


csv_fields = ['#SampleId', 'WorklistId', 'Variant', 'Genotype', f'{proband_id}_DP', f'{proband_id}_GQ',  f'{proband_id}_AD', 'SYMBOL', 'worst_consequence', 'Consequence', 'inheritance_models',  'HGVSc', 'HGVSp', 'CLIN_SIG', 'Existing_variation', 'AutoPick', 'gnomad_popmax_af'] + gt_fields

# save to file
variant_df[csv_fields].to_csv(output, index=False, sep='\t')