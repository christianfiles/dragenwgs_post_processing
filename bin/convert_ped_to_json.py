#!/usr/bin/env python

import json
import argparse
import csv

"""
Script to convert tsv ped file into json verison for qiagen upload

"""

parser = argparse.ArgumentParser(description='Create a json pedigree for upload to Qiagen file from text ped file.')
parser.add_argument('--pedfile', type=str, nargs=1, required=True,
				help='ped file location')
parser.add_argument('--seqid', type=str, nargs=1, required=True,
				help='sequencing id')
parser.add_argument('--output', type=str, nargs=1, required=True,
				help='output name')

args = parser.parse_args()

pedfile = args.pedfile[0]
seqid = args.seqid[0]
output = args.output[0]

ped_dict = {}

# create a dictionary with ped information
with open(args.pedfile[0], 'r') as csvfile:

	my_reader = csv.reader(csvfile, delimiter='\t')

	for row in my_reader:

		family_id = row[0]
		sample_id = row[1]

		if row[2] != '0' and row[3] != '0':

			is_proband = True
			dad = row[2]
			mum = row[3] 

		else:

			is_proband = False

		if is_proband == True: 

			row_dict = {'id': f'{seqid}_{row[1]}', 'proband': is_proband, 'father': f'{seqid}_{dad}', 'mother': f'{seqid}_{mum}', 'externalId': f'{seqid}_{row[1]}'}

		else:

			row_dict = {'id': f'{seqid}_{row[1]}', 'externalId': f'{seqid}_{row[1]}'}

		if family_id == '0':

			ped_dict[sample_id] = [row_dict]

		else:

			if family_id not in ped_dict:

				ped_dict[family_id] = [row_dict]

			else:

				ped_dict[family_id].append(row_dict)

# loop through ped dict
for family_id in ped_dict:

	proband_id = None

	# work out which is proband
	for member in ped_dict[family_id]:

		if 'proband' in member:

			proband_id = member['id']

	# decide what to call file
	if proband_id != None:

		file_name = f'{output}/{proband_id}.json'

	else:

		name = ped_dict[family_id][0]['id']

		file_name = f'{output}/{name}.json'


	# write file
	with open(file_name, 'w') as outfile:


		# try and order it
		new_ped_dict_family = []

		# loop through and put proband first
		for member in ped_dict[family_id]:

			if 'proband' in member:

				new_ped_dict_family.append(member)

		# now add other family members
		for member in ped_dict[family_id]:

			if 'proband' not in member:

				new_ped_dict_family.append(member)

		# write json
		json.dump(new_ped_dict_family, outfile)



















