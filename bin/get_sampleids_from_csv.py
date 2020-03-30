#!/usr/bin/env python

import argparse
import csv

parser = argparse.ArgumentParser(description='Get sample ids from the csv report')
parser.add_argument('--csv', type=str, nargs=1, required=True,
				help='csv report')

args = parser.parse_args()

csv_file = args.csv[0]

sample_ids = []


with open(csv_file, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',')

	count = 0

	for row in spamreader:

		if count == 0:

			# header row
			# sop get family ids
			for key in row:

				if key[-3:] == '_AD':

					sample_ids.append(key.split('_')[0])

		count = count +1



print (' '.join(sample_ids) )

