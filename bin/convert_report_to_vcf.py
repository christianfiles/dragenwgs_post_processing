#!/usr/bin/env python

"""
Scrip to convert CSV file to VCF

1) Loop through CSV and get variants and family member ids

2) For each family member id loop through vcf and print out just those variants.

3) Then we have a vcf with those variants in csv in.

It's a bit of a hack but should work!

No quality filters applied as assuming we have already applied those in CSV report

"""

from pysam import VariantFile
import argparse
import csv

"""
Script to apply filters to WGS data - requires a VCF with only one sample. This is to be done before upload to Qiagen

"""

parser = argparse.ArgumentParser(description='filter a vcf based on variants in a csv.')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='vcf file location')
parser.add_argument('--csv', type=str, nargs=1, required=True,
				help='csv report')
parser.add_argument('--sample_id', type=str, nargs=1, required=True,
				help='sample id')

args = parser.parse_args()

vcf = args.vcf[0]
csv_file = args.csv[0]
sample_id = args.sample_id[0]

# read csv report and get family members and variants
variant_dict = {}

with open(csv_file, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',')

	count = 0

	for row in spamreader:

		if count == 0:

			pass

		else:

			variant_dict[row[2]] = row[2]

		count = count +1

# now loop through vcf

myvcf = VariantFile(vcf, "r")
print(myvcf.header, end='')

for variant in myvcf:

	chrom = variant.chrom
	pos = variant.pos
	ref = variant.ref
	alt = variant.alts[0]

	variant_id  = f'{chrom}:{pos}{ref}>{alt}'

	if variant_id in variant_dict:
		print (variant, end='')
