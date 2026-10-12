#!/usr/bin/env python3

import optparse
from collections import defaultdict

################################# Command line options

desc='Select the best summit from narrowPeak file of MACS3 using the lowest p-value'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 01-06-2024 - Author: FCicconardi')

parser.add_option('-b', '--BED', dest='bed', help='BED4+6 narrowPeak file from MACS3.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['bed']
for m in mandatories:
	if not opts.__dict__[m]:
		print("\nWARNING! One or more options not specified\n")
		parser.print_help()
		exit(-1)

#####################################################


peakDCT=defaultdict(list)

with open(opts.bed) as bed:
	bed=bed.readlines()
	for peak in bed:
		peak=peak.strip().split()
		coord=f"{peak[0]}:{peak[1]}-{peak[2]}"
		peakDCT[coord].append('\t'.join(peak))



for peak in peakDCT:
	besSummit=defaultdict(list)
	for summit in peakDCT[peak]:
		score=int(summit.split()[4])
		besSummit[score].append(summit)
	print(besSummit[max(besSummit.keys())][0])
		
