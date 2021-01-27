### read the TCGA segmentation file

import sys
import numpy as np

path = sys.argv[1]
sys.path.insert(1, path) #path to the settings file

import settings


#for each sample, match it with the WGS name.

#get this from the settings
cnvFile = '../data/cnvs/brca/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'
cnvFile = '../data/cnvs/luad/LUAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'
cnvFile = '../data/cnvs/coad/COAD.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'
cnvFile = '../data/cnvs/ov/OV.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'
cnvFile = '../data/cnvs/lihc/LIHC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'

#read the icgc metadata file to get the WGS name for this sample
#if we have TCGA expression data to link with PCAWG, the samples need to be mapped to those IDs.
#read the metadata file from icgc
def getMetadataICGC(metadataFile):
	#get the metadata file to extract the mapping from wgs to rna-seq identifiers

	nameMap = dict()
	with open(metadataFile, 'r') as inF:

		header = dict()
		lineCount = 0
		for line in inF:
			line = line.strip()
			splitLine = line.split('\t')
			if lineCount < 1:

				for colInd in range(0, len(splitLine)):
					header[splitLine[colInd]] = colInd

				lineCount += 1
				continue

			#tcga_sample_uuid is the name we have in the SV directory
			#sample_id is the one at TCGA in the expression data.
			if header['tcga_sample_uuid'] < len(splitLine):
				wgsName = splitLine[header['tcga_sample_uuid']]
				expressionName = splitLine[header['sample_id']]

				expressionName = expressionName.split('-')
				newExpressionName = '-'.join(expressionName[0:4])

				nameMap[newExpressionName] = wgsName

	return nameMap

nameMap = getMetadataICGC(settings.files['metadataICGC'])

#then make a new file with the WGS sample names, and the actual CN rather than the log.
convertedCNVs = []
with open(cnvFile, 'r') as inF:

	lineCounter = 0
	for line in inF:

		if lineCounter < 1: #skip header
			lineCounter += 1
			continue

		line = line.strip()
		splitLine = line.split('\t')
		sample = splitLine[0]

		sample = sample.split('-')
		newSampleName = '-'.join(sample[0:4])

		if newSampleName not in nameMap:
			continue

		convertedSampleName = nameMap[newSampleName]
		
		segMean = float(splitLine[len(splitLine)-1])
		
		cn = (2 ** segMean)*2
		
		if cn > 2.3 or cn < 1.7:
			#new sample name, chr, start, end, probes, cn
			convertedCNVs.append([convertedSampleName, 'chr' + splitLine[1], splitLine[2],
								  splitLine[3], splitLine[4], str(cn)])

outFile = sys.argv[2]
with open(outFile, 'w') as outF:

	outF.write('Sample\tChromosome\tStart\tEnd\tNum_Probes\tCN\n')

	for line in convertedCNVs:
		outF.write('\t'.join(line))
		outF.write('\n')

