import argparse
import pandas as pd
from pyvariantfilter.family import Family
from pyvariantfilter.family_member import FamilyMember
from pyvariantfilter.variant_set import VariantSet

parser = argparse.ArgumentParser(description='Create a text report with priority')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='vcf path.')
parser.add_argument('--proband_id', type=str, nargs=1, required=True,
				help='proband id.')
parser.add_argument('--ped', type=str, nargs=1, required=True,
				help='ped file containing family relationships.')
parser.add_argument('--gene_list', type=str, nargs=1, required=True,
				help='list of genes for virtual panel.')
parser.add_argument('--gnomad_ad', type=float, nargs=1, required=True,
				help='max gnomad freq to filter on for dominant variants.')
parser.add_argument('--gnomad_r', type=float, nargs=1, required=True,
				help='max gnomad freq to filter on for reccesive variants.')
parser.add_argument('--minqual_snps', type=float, nargs=1, required=True,
				help='min quality for snps')
parser.add_argument('--minqual_indels', type=float, nargs=1, required=True,
				help='min quality for indels')
parser.add_argument('--min_dp', type=int, nargs=1, required=True,
				help='min depth for variants')
parser.add_argument('--min_gq', type=int, nargs=1, required=True,
				help='min gq for variants')
parser.add_argument('--min_af_mt', type=float, nargs=1, required=True,
				help='min af for MT variants')
parser.add_argument('--output', type=str, nargs=1, required=True, help='output name')

args = parser.parse_args()

vcf = args.vcf[0]
proband_id = args.proband_id[0]
ped = args.ped[0]
gene_list = args.gene_list[0]
gnomad_ad = args.gnomad_ad[0]
gnomad_r = args.gnomad_r[0]
minqual_snps = args.minqual_snps[0]
minqual_indels = args.minqual_indels[0]
min_dp = args.min_dp[0]
min_gq = args.min_gq[0]
min_af_mt = args.min_af_mt[0]
output_name = args.output[0]

initial_af = max(gnomad_ad, gnomad_r)


# define a few functions to help

def passes_initial_filter(variant, proband_id, gene_dict):
	"""
	Filter variants from the VCF.
	
	We import if the variant passes quality filtering and is below 1% in gnomad exomes and gnomad genomes AND
	
	b) Is listed as pathogenic in clinvar OR
	c) Has a a relevant consequence
	
	"""




	# If the proband has the variant and we pass the genotype and variant level filters
	if variant.has_alt(proband_id) and variant.passes_gt_filter(proband_id, min_gq=min_gq):

		if variant.chrom == 'MT':

			af = variant.info_annotations

		else:

			if variant.is_snp():

				if variant.quality < minqual_snps:

					return False

			else:

				if variant.quality < minqual_indels:

					return False



		# The filter_on_numerical_transcript_annotation_lte() function allows us to filter on numerical values 
		# we can set different cutoffs for different variant types. For example ad_het is variants in which the 
		# proband is heterozygous on an autosome. In this case we get two boolean values describing whether the 
		# variant is below 1% in the gnomad genomes and gnomad exomes datasets.
		freq_filterg = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADg_AF_POPMAX',
																						  ad_het=initial_af,
																						  ad_hom_alt=initial_af,
																						  x_male =initial_af,
																						  x_female_het=initial_af,
																						  x_female_hom=initial_af,
																						  compound_het=initial_af,
																						  y=initial_af,
																						  mt=initial_af,
																						  )
		freq_filtere = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADe_AF_POPMAX',
																						  ad_het=initial_af,
																						  ad_hom_alt=initial_af,
																						  x_male =initial_af,
																						  x_female_het=initial_af,
																						  x_female_hom=initial_af,
																						  compound_het=initial_af,
																						  y=initial_af,
																						  mt=initial_af,
																						  )  
		
		genes = variant.get_genes(feature_key='SYMBOL')

		in_panel = False

		for gene in genes:

			if gene in gene_dict:

				in_panel = True
				break

		# Coopt the get_genes() function to get the clinvar annotation VEP field.
		clinvar = variant.get_genes(feature_key='CLIN_SIG')
		is_path_in_clinvar = False
		
		for anno in clinvar:
			
			if 'pathogenic' in anno.lower():
				is_path_in_clinvar = True
				break
				
		# If the variant is below 1% and pathogenic in clinvar then import
		if freq_filterg and freq_filtere and is_path_in_clinvar and in_panel:
			
			return True
		
		csq_filter = False
		
		if variant.get_worst_consequence() in {'transcript_ablation': None,
											   'splice_acceptor_variant': None,
											   'splice_donor_variant': None,
											   'stop_gained': None,
											   'frameshift_variant': None,
											   'stop_lost': None,
											   'splice_region_variant': None,
											   'start_lost': None,
											   'transcript_amplification': None,
											   'inframe_insertion': None,
											   'inframe_deletion': None,
											   'missense_variant': None,
											   'protein_altering_variant': None,
											   'incomplete_terminal_codon_variant': None,
											   'start_retained_variant': None,
											   'stop_retained_variant': None}:
		
			csq_filter = True
		
	   # If the variant is below 1% and has a relevant consequence then import
		if csq_filter and freq_filterg and freq_filtere and in_panel:
			
			return True
		
	return False

def passes_final_filter_trio(variant, compound_het_dict , inheritance):

	freq_filterg = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADg_AF_POPMAX',
																						  ad_het=gnomad_ad,
																						  ad_hom_alt=initial_af,
																						  x_male =initial_af,
																						  x_female_het=initial_af,
																						  x_female_hom=initial_af,
																						  compound_het=initial_af,
																						  y=initial_af,
																						  mt=initial_af,
																						  compound_het_dict=compound_het_dict
																						  )
	freq_filtere = variant.filter_on_numerical_transcript_annotation_lte(annotation_key='gnomADe_AF_POPMAX',
																						  ad_het=gnomad_ad,
																						  ad_hom_alt=initial_af,
																						  x_male =initial_af,
																						  x_female_het=initial_af,
																						  x_female_hom=initial_af,
																						  compound_het=initial_af,
																						  y=initial_af,
																						  mt=initial_af,
																						  compound_het_dict=compound_het_dict
																						  )  



	# Coopt the get_genes() function to get the clinvar annotation VEP field.
	clinvar = variant.get_genes(feature_key='CLIN_SIG')
	is_path_in_clinvar = False
	
	for anno in clinvar:
		
		if 'pathogenic' in anno.lower():
			is_path_in_clinvar = True
			break
			
	# If the variant is below 1% and pathogenic in clinvar then import
	if freq_filterg and freq_filtere and is_path_in_clinvar:
		
		return True

	
	# Get variants which match certain inheritance models
	if variant.matches_inheritance_model(inheritance,
										   compound_het_dict,
										   min_parental_gq_dn = min_gq,
										   min_parental_depth_dn = min_dp,
										   min_parental_gq_upi = min_gq,
										   min_parental_depth_upi = min_dp
										   ) and freq_filterg and freq_filtere:
		
		return True
		
	return False




# read ped into df
ped_df = pd.read_csv(ped, sep='\t', names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'affected'])


# read gene list to dict
gene_df = pd.read_csv(gene_list, sep=',')

gene_dict = {}

for row in gene_df.itertuples():
	
	gene_dict[row[1]] = row[1]


# have we got a singleton or a family?
filtered_ped = ped_df[ped_df['sample_id']==proband_id]

family_id = filtered_ped['family_id'].iloc[0]

affected = filtered_ped['affected'].iloc[0]

sex = filtered_ped['sex'].iloc[0]

has_family= False

if family_id != '0':
	
	has_family = True
	
is_affected = False
	
if affected == 2:
	
	is_affected = True


# create family object
if has_family == True and is_affected == True:
	
	my_family = Family(family_id)
	
	my_family.read_from_ped_file(ped, family_id, proband_id)

elif has_family == False and is_affected == True:
	
	proband = FamilyMember(proband_id, proband_id, sex, is_affected)

	my_family = Family(proband_id)

	my_family.add_family_member(proband)

	my_family.set_proband(proband.get_id())


# Create a new VariantSet object
my_variant_set = VariantSet()

# Associate the my_family object with my_variant_set
my_variant_set.add_family(my_family)


my_variant_set.read_variants_from_vcf(vcf,
									proband_variants_only=True,
									filter_func=passes_initial_filter,
									args=(proband_id, gene_dict))

# see whether we can phase comp hets by inheritance

if my_variant_set.family.proband_has_both_parents() == True:
	
	# Create an attribute my_variant_set.candidate_compound_het_dict where each transcript is a key the variants 
	# within that transcript are the values
	my_variant_set.get_candidate_compound_hets()

	# As we have both parents we can filter the compound hets.
	my_variant_set.filter_compound_hets()

	# Flatten the filtered compound hets (my_variant_set.filtered_compound_het_dict) into a a dictionary with each 
	# genuine compound het as a key.
	my_variant_set.get_filtered_compound_hets_as_dict()


	inheritance = ['uniparental_isodisomy','autosomal_reccessive','x_reccessive','x_dominant','de_novo', 'compound_het']
	
else:
	
	
	# Create an attribute my_variant_set.candidate_compound_het_dict where each transcript is a key the variants 
	# within that transcript are the values
	my_variant_set.get_candidate_compound_hets()

	# As we are pretending we do not have any parents we cannot phase the compound hets
	my_variant_set.get_unfiltered_compound_hets_as_dict()

	inheritance = ['uniparental_isodisomy', 'autosomal_dominant', 'autosomal_reccessive','x_reccessive','x_dominant','de_novo', 'compound_het']

my_variant_set.filter_variants(passes_final_filter_trio, args=(my_variant_set.final_compound_hets, inheritance ))


variant_df = my_variant_set.to_df()



gt_fields = []

for fm in my_family.get_all_family_member_ids():

	for field in ['_GT', '_DP', '_GQ']:


		gt_fields.append(fm + field)


csv_fields = ['variant_id', 'csq_SYMBOL', 'csq_Feature', 'csq_HGVSc', 'csq_HGVSp', 'csq_CLIN_SIG', 'csq_Existing_variation', 'worst_consequence', 'inheritance_models', 'csq_PICK', 'filter_status', 'csq_gnomADg_AF_POPMAX', 'csq_gnomADe_AF_POPMAX'] + gt_fields


variant_df[csv_fields].to_csv(output_name, index=False)










