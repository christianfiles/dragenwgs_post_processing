#!/usr/bin/env python

import argparse
import pandas as pd
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet


parser = argparse.ArgumentParser(description='Create a text report with priority')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='The path to the VCF file. Must be bgzipped and tabixed.')
parser.add_argument('--proband_id', type=str, nargs=1, required=True,
				help='The Sample ID of the proband in the VCF.')
parser.add_argument('--ped', type=str, nargs=1, required=True,
				help='A ped file describing the family relationships.')
parser.add_argument('--max_mitomap_af', type=float, nargs=1, required=True,
				help='The maximum mitomap population frequency to filter on.')
parser.add_argument('--min_mt_af', type=float, nargs=1, required=True,
				help='The minimum AF to filter on.')
parser.add_argument('--whitelist', type=str, nargs=1, required=True,
				help='A list of variants to NOT be filtered out regardless of other criteria. Should be A CSV with the chromosome and position as first and second columns respectively.')
parser.add_argument('--output', type=str, nargs=1, required=True, help='The output name')
parser.add_argument('--worklist', type=str, nargs=1, required=True, help='Worklist ID for Audit purposes. Can be used as a comment.')

args = parser.parse_args()


vcf = args.vcf[0]
proband_id = args.proband_id[0]
ped = args.ped[0]
output_name = args.output[0]
worklist = args.worklist[0]
whitelist = args.whitelist[0]
min_mt_af = args.min_mt_af[0]
max_mitomap_af = args.max_mitomap_af[0]

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

def pathogenic_in_clinvar(clinvar_vep, clinvar_custom, clinvar_conflicting):
	"""
	Check multiple fields to see if pathogenic in clinvar.

	We check:

	1) CLIN_SIG - added by VEP
	2) CLIN_SIG - added by our own clinvar annnotation.
	3) CLINVAR_CONFLICTING - added by our own clinvar annotation.

	We check all three for presense of string pathogenic.


	"""
	for anno in clinvar_vep:


		if 'pathogenic' in anno.lower():

			return True

	for anno in clinvar_custom:

		# ignore this as we will pick it up via clinvar_conflicting
		if anno == 'Conflicting_interpretations_of_pathogenicity':

			pass

		elif 'pathogenic' in anno.lower():

			return True

	for anno in clinvar_conflicting:

		if 'pathogenic' in anno.lower():

			return True


	return False

def passes_initial_filter(variant, proband_id, max_mitomap, min_af):
	"""
	An initial filter to 

	We import if the variant passes quality filtering and is below 1% in gnomad exomes and gnomad genomes AND

	a) Is listed as pathogenic in clinvar OR
	b) Has a a relevant consequence

	OR is in the whitelist

	"""




	# If the proband has the variant and we pass the genotype and variant level filters
	if variant.has_alt(proband_id):

		# ignore mt variants for now
		if variant.chrom != 'MT':

			return False

		else:
			
			reads_ref = variant.get_ref_reads(proband_id) 
			reads_alt = variant.get_alt_reads(proband_id)

			# no reads
			if (reads_ref + reads_alt) == 0:

				return False
			
			af = reads_alt / (reads_ref + reads_alt)
			
			if af < min_af:
				
				return False

		mito_map = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='mitomap_AF',
																						  ad_het=max_mitomap,
																						  ad_hom_alt=max_mitomap,
																						  x_male =max_mitomap,
																						  x_female_het=max_mitomap,
																						  x_female_hom=max_mitomap,
																						  compound_het=max_mitomap,
																						  y=max_mitomap,
																						  mt=max_mitomap,
																				   )
		
	
		# If the variant is below x%
		if mito_map:

			return True

	return False

def get_af(df, sample_id):
	"""
	Convert AD to AF

	e.g. alt /ref+alt

	"""

	ref = df[f'{sample_id}_AD'].split(',')[0]
	alt = df[f'{sample_id}_AD'].split(',')[1]

	if float(ref)+float(alt) == 0:

		return 0

	return round(float(alt) / (float(ref)+float(alt)), 3)



args = parser.parse_args()

# read ped into df
ped_df = pd.read_csv(ped, sep='\t', names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'affected'])

print (ped_df)

# read white list
white_df = pd.read_csv(whitelist)
white_dict = {}

for row in white_df.itertuples():
	
	white_dict[f'{row.chrom}:{row.pos}'] = f'{row.chrom}:{row.pos}'

print (proband_id)
# have we got a singleton or a family?
filtered_ped = ped_df[ped_df['sample_id']==proband_id]

family_id = filtered_ped['family_id'].iloc[0]

affected = filtered_ped['affected'].iloc[0]

sex = filtered_ped['sex'].iloc[0]


# if any family members don't have sex set in ped
family_ped = ped_df[ped_df['family_id'] == family_id]


if sex =='0' or sex == 0:

	print('Sex cannot be zero - not running program creating emprty file.')
	f = open(output_name, 'w')
	f.write(f'No sex for sample {proband_id}. Program not run.')
	f.close()
	exit()

has_family= True

if family_id == '0' or family_id == 0:
	
	has_family = False
	
is_affected = False
	
if affected == 2:
	
	is_affected = True

# if we have a family and sample is affected
if has_family == True and is_affected == True:
	
	my_family = Family(family_id)
	
	my_family.read_from_ped_file(ped, family_id, proband_id)

elif has_family == False and is_affected == True:
	
	proband = FamilyMember(proband_id, proband_id, sex, is_affected)

	my_family = Family(proband_id)

	my_family.add_family_member(proband)

	my_family.set_proband(proband.get_id())

else:
	print('not running as not affected.')
	f = open(output_name, 'w')
	f.write(f'Sample {proband_id} is not affected. Program not run.')
	f.close()
	exit()


# Create a new VariantSet object
my_variant_set = VariantSet()

# Associate the my_family object with my_variant_set
my_variant_set.add_family(my_family)

# read from vcf
my_variant_set.read_variants_from_vcf(vcf,
									proband_variants_only=True,
									filter_func=passes_initial_filter,
									args=(proband_id, max_mitomap_af, min_mt_af))


# convert to pandas dataframe
variant_df = my_variant_set.to_df()

# catch NTC
if variant_df.shape[0] == 0:

	f = open(output_name, 'w')
	f.write(f'No variants for sample {proband_id}. Program not run.')
	f.close()
	exit()

# convert genotype fields
for fm in my_family.get_all_family_member_ids():

	variant_df[f'{fm}_Genotype'] = variant_df.apply(fix_genotype, axis=1, args=( f'{fm}_GT',))

gt_fields = []

proband_id = my_family.get_proband_id()
variant_df['#SampleId'] = proband_id
variant_df['WorklistId'] = worklist
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
variant_df['mitomap_AF'] = variant_df['csq_mitomap_AF']

# add sample specific columns such as genotype and depth
for fm in my_family.get_all_family_member_ids():

	# ignore proband as we did it earlier
	if fm != proband_id:

		for field in ['_Genotype', '_DP', '_AF']:

			gt_fields.append(fm + field)


for fm in my_family.get_all_family_member_ids():

	variant_df[f'{fm}_AF'] = variant_df.apply(get_af,axis=1,args=(fm,))



csv_fields = ['#SampleId', 'WorklistId', 'Variant', 'Genotype', f'{proband_id}_DP', f'{proband_id}_AF', 'SYMBOL', 'Feature', 'worst_consequence', 'Consequence', 'inheritance_models',  'HGVSc', 'HGVSp', 'CLIN_SIG', 'Existing_variation', 'csq_PICK', 'mitomap_AF' ] + gt_fields

# save to file
variant_df[csv_fields].to_csv(output_name, index=False, sep='\t')