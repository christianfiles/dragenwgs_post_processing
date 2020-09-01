#!/usr/bin/env python

import argparse
import csv

"""
Read in the variant report from the test data and check all the variants are there
"""

parser = argparse.ArgumentParser(description='Create a text report with priority')
parser.add_argument('--variant_report', type=str, nargs=1, required=True,
				help='The path to the variant report.')
parser.add_argument('--mt_variant_report', type=str, nargs=1, required=True,
				help='The path to the MT variant report.')

args = parser.parse_args()


variant_report = args.variant_report[0]
mt_variant_report = args.mt_variant_report[0]

expected_variants = {
	'2:167060869C>A': 'compound_het',
	'2:167144946C>A': 'compound_het',
	'15:23890978G>A': '',
	'22:37538526A>G': 'autosomal_reccessive',
	'15:89865989C>G': 'compound_het',
	'15:89870432C>T': 'compound_het',
	'4:56878134G>GA': 'compound_het',
	'4:56885717A>T': 'compound_het',
	'13:111277601C>T': 'autosomal_reccessive',
	'19:1398963C>T': 'autosomal_reccessive|uniparental_isodisomy',
	'19:1399792C>T': '',
	'2:122288506G>A': 'compound_het',
	'2:122288566G>T': 'compound_het',
	'17:56283913GGCATGCCATTGGGACAGCCTCAGGTTTCT>G': 'compound_het',
	'17:56293449C>T': 'compound_het',
	'12:48370516T>A': 'autosomal_dominant|de_novo',
	'5:125882068C>G': 'compound_het',
	'5:125903988C>G': 'compound_het',
	'X:153296882G>A': 'x_reccessive'




	}

expected_variants_mt ={
	'MT:8356T>C': 'mitochrondrial',
	'MT:10010T>C': 'mitochrondrial'
}


found_variants = {}
found_variants_mt = {}

with open(variant_report) as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:

		found_variants[row[2]] = row[10]


with open(mt_variant_report) as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:

		found_variants_mt[row[2]] = row[10]


print ('Testing Variant Reports')
for variant in expected_variants:

	if variant not in found_variants:

		raise Exception(f'Did not find variant {variant}')

	if found_variants[variant] != expected_variants[variant]:

		raise Exception(f'Wrong inheritance for {variant}')


print ('Testing MT Variant Reports')
for variant in expected_variants_mt:

	if variant not in found_variants_mt:

		raise Exception(f'Did not find variant {variant}')
		
	if found_variants_mt[variant] != expected_variants_mt[variant]:

		raise Exception(f'Wrong inheritance for {variant}')

print ('All Tests Pass')