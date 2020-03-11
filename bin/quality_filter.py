#!/usr/bin/env python

from pysam import VariantFile
import argparse



parser = argparse.ArgumentParser(description='Perform quality filtering on single sample vcfs before we upload to Qiagen')
parser.add_argument('--vcf', type=str, nargs=1, required=True,
				help='vcf file location')
parser.add_argument('--sampleid', type=str, nargs=1, required=True,
				help='sample id')
parser.add_argument('--min_mt_af', type=float, nargs=1, required=True,
				help='The minimum allele frequency a MT variant must have to not fail e.g. 0.05')
parser.add_argument('--snp_qual', type=float, nargs=1, required=True,
				help='The minimum QUAL for a (non MT) SNP')
parser.add_argument('--indel_qual', type=float, nargs=1, required=True,
				help='The minimum QUAL for a (non MT) Indel')
parser.add_argument('--min_dp', type=int, nargs=1, required=True,
				help='The minimum sample DP (non MT) variant')
parser.add_argument('--min_gq', type=int, nargs=1, required=True,
				help='The minimum sample GQ (non MT) variant')


args = parser.parse_args()

vcf = args.vcf[0]
sampleid = args.sampleid[0]
min_mt_af = args.min_mt_af[0]
snp_qual = args.snp_qual[0]
indel_qual = args.indel_qual[0]
min_dp = args.min_dp[0]
min_gq = args.min_gq[0]


myvcf = VariantFile(vcf, "r")

myvcf.header.filters.add('LowAFMT', None, None, f'An MT variant with AF below {min_mt_af}')
myvcf.header.filters.add('SNPQual', None, None, f'An SNP variant with QUAL below {snp_qual}')
myvcf.header.filters.add('IndelQual', None, None, f'An Indel variant with QUAL below {indel_qual}')
myvcf.header.filters.add('LowDP', None, None, f'A non MT variant with DP below {min_dp}')
myvcf.header.filters.add('LowGQ', None, None, f'A non MT variant with GQ below {min_gq}')
myvcf.header.filters.add('lod_fstar', None, None, f'MT LOD is below default')
myvcf.header.filters.add('DRAGENHardQUAL', None, None, f'Default Hard Filter Applied to Genotypes by Dragen')

print(myvcf.header, end='')

for variant in myvcf:

	chrom = variant.chrom
	ref = variant.ref
	alt = variant.alts
	qual = variant.qual

	if chrom == 'MT':

		filters_to_add = []

		variant.filter.clear()

		ft = variant.samples[sampleid]['FT']
		af  = variant.samples[sampleid]['AF']

		if len(af) >1:

			print(variant)
			continue

		af = af[0]

		if af < min_mt_af:

			filters_to_add.append('LowAFMT')

		if ft != 'PASS':

			filters_to_add.append(ft)

		if len(filters_to_add) == 0 and ft == 'PASS':

			variant.filter.add('PASS')

		else:

			for filt in filters_to_add:

				variant.filter.add(filt)

	else:

		variant.filter.clear()

		filters_to_add = []

		try:

			gq = variant.samples[sampleid]['GQ']

		except:

			print (variant)
			continue

		if gq == None:

			gq = 0

		try:

			dp = variant.samples[sampleid]['DP']

		except:

			print (variant)
			continue

		if dp == None:

			dp = 0


		is_snp = True

		if len(ref) > 1:

			is_snp = False

		# filter snps
		for allele in alt:

			if len(alt) > 1:

				is_snp = False

		if is_snp == True:

			if qual < snp_qual:

				filters_to_add.append('SNPQual')

		else:

			if qual < indel_qual:

				filters_to_add.append('IndelQual')

		# filter dp

		if dp < min_dp:

			filters_to_add.append('LowDP')

		if gq < min_gq:

			filters_to_add.append('LowGQ')

		if len(filters_to_add) == 0:

			variant.filter.add('PASS')

		else:

			for filt in filters_to_add:

				variant.filter.add(filt)

		
	print (variant, end='')














