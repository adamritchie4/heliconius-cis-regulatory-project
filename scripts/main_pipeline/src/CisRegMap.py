#!/usr/bin/env python3

import optparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binomtest
from collections import defaultdict

################################# Command line options

desc='''This script processes TSS + gene body features to define regulatory zones (promoter, proximal, and gene body) from an annotation file 
and associates peaks with genes and regulatory zones. It generates two output tables: one associating peaks with genes 
and another associating peaks with regulatory zones. Additionally, it provides a summary of how the regulatory zones 
were defined.'''

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 09-10-2024 - Author: FCicconardi')

parser.add_option('-t', '--TSS-BED', dest='tbed', help='Path to the transcription start sites (TSS) BED file generated with get_TSSandGeneBody.py.', action='store', metavar='FILE')
parser.add_option('-e', '--exon-BED', dest='ebed', help='Path to the exon BED file.', action='store', metavar='FILE')
parser.add_option('-i', '--intron-BED', dest='ibed', help='Path to the intron BED file.', action='store', metavar ='FILE')
parser.add_option('-p', '--Peak-BED', dest='pbed', help='Path to the peaks BED file.', action='store', metavar='FILE')
parser.add_option('-r', '--RandomisedPeaks-BED', dest='rbed', help='Randomised set of peaks in BED3 file format.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['tbed','pbed','ebed','ibed']
for m in mandatories:
	if not opts.__dict__[m]:
		print("\nWARNING! One or more options not specified\n")
		parser.print_help()
		exit(-1)

############################## Reading files and parameters from sys import argv


###################################### Functions

def ReadBedLine(line):
	el=line.strip().split('\t')
	scf=el[0]
	TSS=int(el[1])
	gene=el[3]
	strand=el[5]
	if strand == '+':
		Gbody=[int(el[6]),int(el[7])]
	else: Gbody=[int(el[7]),int(el[6])]

	return scf,TSS,gene,strand,Gbody

def ReadRegionBedLine(line):
    el=line.strip().split('\t')
    scf=el[0]
    start=int(el[1])
    end=int(el[2])

    return scf,start,end



def DefPromoter(TSS,strand,Gbody):
	if strand == '+':
		promoter_st = TSS-1000
		promoter_en = TSS+1000
		if promoter_en > Gbody[1]: promoter_en = Gbody[1]
		if promoter_st < 0: promoter_st = 0
	elif strand == '-':
		promoter_st = TSS+1000
		promoter_en = TSS-1000
		if promoter_st > Gbody[1]: promoter_st = Gbody[1]
		if promoter_en < 0: promoter_en = 0

	region = sorted([promoter_st,promoter_en])
	return region
		
def DefProximal(TSS,strand,Prom_region):
	if strand == '+':
		prox_st = Prom_region[0]-4000
		if prox_st < 0: prox_st = 0
		region = [prox_st,Prom_region[0]]
	elif strand == '-':
		prox_st = Prom_region[1]+4000
		region = [Prom_region[1],prox_st]

	return region

def DefGeneBody(Gbody,Prom_region,strand):
	if strand == '+':
		GB_st = Prom_region[1]
		GB_end = Gbody[1]
	elif strand == '-':
		GB_st = Prom_region[0]
		GB_end = Gbody[0]

	region = sorted([GB_st,GB_end])
	return region

def CheckOverlapRegion(Regions,scf,pst,pend):
	Reg_gene='-'

	for region in Regions[scf]:
		gene=region[0]
		reg_st=region[1][0]
		reg_end=region[1][1]

		start1, end1 = pst,pend
		start2, end2 = region[1]

		#print(start1, end1, 'VS', start2, end2, gene)

		# Check if the regions overlap
		if start1 <= end2 and end1 >= start2:
			#print("Overlap!!!")
			Reg_gene=gene
			break

	return Reg_gene

def CheckOverlapBed(Regions,scf,pst,pend):

    if scf not in Regions:
        return False

    for start, end in Regions[scf]:
        if pst < end and pend > start:
            return True

    return False

###################################### Functions end ####################

###################################### Processing data --->
classes=['Promoter','Proximal','Exon','Intron','Distal']

BedDCT=defaultdict(list)
PeakDCT=defaultdict(list)

Promoter=defaultdict(list)
Proximal=defaultdict(list)
GeneBody=defaultdict(list)

Exon=defaultdict(list)
Intron=defaultdict(list)

PeakClassification=defaultdict(list)
RandomPeakDist=defaultdict(list)
RandomPeakDist_prop=defaultdict(list)

species=opts.tbed.split('.')[0]

print('Reading the TSS file:', opts.tbed, '...')

with open(opts.tbed) as tss_file:
	tss_lines=tss_file.readlines()
	for line in tss_lines:
		scf,TSS,gene,strand,Gbody=ReadBedLine(line)
		BedDCT[gene].append((scf,TSS,strand,Gbody[0],Gbody[1]))

		# 1. define Promoter region
		Prom_region = DefPromoter(TSS,strand,Gbody)
		Promoter[scf].append((gene,Prom_region,strand))

		# 2. define Proximal region
		Prox_region = DefProximal(TSS,strand,Prom_region)
		Proximal[scf].append((gene,Prox_region,strand))

		# 3. define Gene body region
		GB_region = DefGeneBody(Gbody,Prom_region,strand)
		GeneBody[scf].append((gene,GB_region,strand))

print('Reading the exon file:', opts.ebed, '...')

with open(opts.ebed) as exon_file:
    for line in exon_file:
        scf,start,end=ReadRegionBedLine(line)

        Exon[scf].append((start,end))

print('Reading the intron file:', opts.ibed, '...')

with open(opts.ibed) as intron_file:
    for line in intron_file:
        scf,start,end=ReadRegionBedLine(line)

        Intron[scf].append((start,end))

class_out=open('Class.tsv', 'w')
gene_out=open('AssGene.tsv', 'w')
#class_out=open(opts.pbed.replace('.tsv','Class.tsv'), 'w')
#gene_out=open(opts.pbed.replace('.tsv','AssGene.tsv'), 'w')


print('Reading and classifing peaks:', opts.pbed, '...')

with open(opts.pbed) as peak_file:
    peak_lines = peak_file.readlines()
    for line in peak_lines:
        if line.startswith('chrom'):
            continue
        el = line.strip().split('\t')
        scf, pst, pend = el[0], int(el[1]), int(el[2])
        pname = f"{scf}:{pst}-{pend}"
        prom_gene = CheckOverlapRegion(Promoter, scf, pst, pend)
        if prom_gene != '-':
            PeakClassification['Promoter'].append((pname, prom_gene))
            prox_gene, gb = '-', '-'
            prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy = 'Promoter', '-', '-', '-', '-'
            PeakDCT[pname].append((scf, pst, pend, prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy, prom_gene, prox_gene, gb))
        else:
            prox_gene = CheckOverlapRegion(Proximal, scf, pst, pend)
            if prox_gene != '-':
                PeakClassification['Proximal'].append((pname, prox_gene))
                gb = '-'
                prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy = '-', 'Proximal', '-', '-', '-'
                PeakDCT[pname].append((scf, pst, pend, prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy, prom_gene, prox_gene, gb))
            else:
                exon = CheckOverlapBed(Exon, scf, pst, pend)
                if exon:
                    gb = CheckOverlapRegion(GeneBody, scf, pst, pend)
                    PeakClassification['Exon'].append((pname, gb))
                    prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy = '-', '-', 'Exon', '-', '-'
                    PeakDCT[pname].append((scf, pst, pend, prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy, prom_gene, prox_gene, gb))
                else:
                    intron = CheckOverlapBed(Intron, scf, pst, pend)
                    if intron:
                        gb = CheckOverlapRegion(GeneBody, scf, pst, pend)
                        PeakClassification['Intron'].append((pname, gb))
                        prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy = '-', '-', '-', 'Intron', '-'
                        PeakDCT[pname].append((scf, pst, pend, prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy, prom_gene, prox_gene, gb))
                    else:
                        gb = '-'
                        PeakClassification['Distal'].append((pname, 'Distal'))
                        prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy = '-', '-', '-', '-', 'Distal'
                        PeakDCT[pname].append((scf, pst, pend, prom_clasfy, prox_clasfy, exon_clasfy, intron_clasfy, dist_clasfy, prom_gene, prox_gene, gb))

with open(opts.rbed) as random_peak_files:
    random_peak_files = random_peak_files.readlines()

rep = 0
for random_peak_file in random_peak_files:
    rep += 1
    with open(random_peak_file.strip()) as random_peak:
        if rep % 500 == 0:
            print('Reading and classifing random set of peaks:', random_peak_file.strip(), '...')
        RandomPeakClassification = defaultdict(list)
        nrp = 0
        random_peak = random_peak.readlines()
        for line in random_peak:
            if line.startswith('chrom'):
                continue
            nrp += 1
            el = line.strip().split('\t')
            scf, pst, pend = el[0], int(el[1]), int(el[2])
            pname = f"{scf}:{pst}-{pend}"
            prom_gene = CheckOverlapRegion(Promoter, scf, pst, pend)
            if prom_gene != '-':
                RandomPeakClassification['Promoter'].append((pname, prom_gene))
            else:
                prox_gene = CheckOverlapRegion(Proximal, scf, pst, pend)
                if prox_gene != '-':
                    RandomPeakClassification['Proximal'].append((pname, prox_gene))
                else:
                    exon = CheckOverlapBed(Exon, scf, pst, pend)
                    if exon:
                        RandomPeakClassification['Exon'].append((pname, 'Exon'))
                    else:
                        intron = CheckOverlapBed(Intron, scf, pst, pend)
                        if intron:
                            RandomPeakClassification['Intron'].append((pname, 'Intron'))
                        else:
                            RandomPeakClassification['Distal'].append((pname, 'Distal'))
        for cls in classes:
            RandomPeakDist[cls].append(len(RandomPeakClassification[cls]))
            RandomPeakDist_prop[cls].append(len(RandomPeakClassification[cls]) / nrp)

print('\n\n\n')

total=len(PeakDCT)

df={'Species': ['RandomSet',species]}
df_prop={'Species': ['RandomSet',species]}

for cls in classes:
        df[cls] = []
        df_prop[cls] = []

for cls in classes:
        expected_prop=np.average(RandomPeakDist_prop[cls])
        observed=len(PeakClassification[cls])

        print(cls, 'observed:', observed, 'Expected proportion:', expected_prop)

        ZoneEnrich=binomtest(observed, n=total, p=expected_prop, alternative='two-sided')
        print(ZoneEnrich)

        df[cls].append(np.average(RandomPeakDist[cls]))
        df[cls].append(observed)

        df_prop[cls].append(expected_prop)
        df_prop[cls].append(observed/total)


class_out.write("PeakName\tLocus\tPeakZone\n")
gene_out.write("PeakName\tLocus\tPromoterGene\tProximalGene\tGeneBodyGene\n")

for peak in PeakDCT:
        scf,pst,pend,prom_clasfy,prox_clasfy,exon_clasfy,intron_clasfy,dist_clasfy,prom_gene,prox_gene,gb=PeakDCT[peak][0]

        for cls in [prom_clasfy,prox_clasfy,exon_clasfy,intron_clasfy,dist_clasfy]:
                if cls != '-':
                        class_out.write(f"{peak}\t{scf}:{pst}-{pend}\t{cls}\n")

        gene_out.write(f"{peak}\t{scf}:{pst}-{pend}\t{prom_gene}\t{prox_gene}\t{gb}\n")

class_out.close()
gene_out.close()


df = pd.DataFrame(df)
df_prop = pd.DataFrame(df_prop)

# Set 'Species' as the index to facilitate easier plotting
df.set_index('Species', inplace=True)
df_prop.set_index('Species', inplace=True)

print('\n\n\n')
print(df)
print('\n')
print(df_prop)

print('\nProportion sums:')
print(df_prop.sum(axis=1))

#for n, x in enumerate([*df.index.values]):
#    print(f'n: {n}, x: {x}')


#for n, x in enumerate([*df.index.values]):
#    for proportion in df_prop.loc[x]:
#        print(f'n: {n}, x: {x}, proportion: {proportion}')



df_prop.plot(kind='bar', 
                    stacked=True, 
                    colormap='tab10', 
                    figsize=(10, 6))

# Add labels and legend
plt.legend(loc="upper right", ncol=2)
plt.xlabel("Species")
plt.ylabel("Proportion of peaks")
plt.title("Proportion and Absolute Number of Peaks by Region")

for n, x in enumerate([*df.index.values]):
    for (proportion, count, y_loc) in zip(df_prop.loc[x],
                                          df.loc[x],
                                          df_prop.loc[x].cumsum()):
                
        plt.text(x=n - 0.17,
                 y=(y_loc - proportion) + (proportion / 2),
                 s=f'{count}\n({np.round(proportion * 100, 1)}%)', 
                 color="black",
                 fontsize=12,
                 fontweight="bold")

# Set x-axis labels horizontally
plt.xticks(rotation=0)

# Save plot to a file
plt.savefig('stacked_peak_distribution.pdf')  # Change the filename and extension as needed
#plt.savefig(opts.pbed.replace('.tsv','.stacked_peak_distribution.pdf'))  # Change the filename and extension as needed


