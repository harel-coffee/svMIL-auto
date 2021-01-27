"""
	Collect all read counts from the GTEx pipeline output and combine them into 1 file

"""

import sys
import numpy as np
import glob
import re
from os import listdir

#folder with read counts
readCountDir = sys.argv[1]

#read the files & keep read counts in memory

readCountSubdirs = [f for f in listdir(readCountDir)]

geneCounts = dict()
genes = [] #keep gene names
samples = [] #sample names
fileCount = 0
for readCountSubdir in readCountSubdirs:

	if len(glob.glob(readCountDir + '/' + readCountSubdir + '/*gene_reads.gct')) < 1: #tmp for failed/not-yet-completed runs
		continue
	readCountFile = glob.glob(readCountDir + '/' + readCountSubdir + '/*gene_reads.gct')[0]
	
	sampleName = re.search('^.*\/(.+)\.A.+$', readCountFile).group(1)
	shortSampleName = sampleName.split('_')[1]

	samples.append(shortSampleName)

	with open(readCountFile, 'r') as inF:

		lineCount = 0
		for line in inF:
			if lineCount < 3:

				if lineCount == 1:
					splitLine = line.split('\t')
					geneCount = splitLine[0]
					print(geneCount)
					if geneCount != '56206':
						print('different file length')
				
				lineCount += 1
				continue
			line = line.strip()
			splitLine = line.split('\t')
			
			geneName = splitLine[1]
			geneId = splitLine[0]
			count = splitLine[2]
			
			if fileCount < 1:
				genes.append([geneId, geneName])
			
			if geneId not in geneCounts:
				geneCounts[geneId] = dict()
			if shortSampleName not in geneCounts[geneId]:
				geneCounts[geneId][shortSampleName] = count

	
	fileCount += 1

genes = np.array(genes, dtype='object')

#write to a combined file
outFile = sys.argv[2]
with open(outFile, 'w') as outF:
	samples = ['Name', 'Description'] + samples
	header = '\t'.join(samples)
	outF.write(header)
	outF.write('\n')

	for gene in genes:
		line = gene[0] + '\t' + gene[1] + '\t'
		for sample in samples:
			if sample in ['Name', 'Description']:
				continue
			geneCount = geneCounts[gene[0]][sample]
			line += geneCount
			line += '\t'
		line += '\n'
		outF.write(line)
		
#the data also need to be merged with the GTEx pipeline output

#load the GTEx counts



