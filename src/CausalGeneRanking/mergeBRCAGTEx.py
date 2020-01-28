"""
	Merge the BRCA and GTEx expression data into 1 big table.
	Make sure that the genes are matched. 

"""

import sys
import numpy as np


#get the gene expression
def getBrcaExpression(expressionFile):
	
	#also use a map for the gene names
	geneNameConversionMap = dict()
	geneNameConversionFile = sys.argv[3]
	with open(geneNameConversionFile, 'r') as inF:
		
		lineCount = 0
		for line in inF:
			
			if lineCount < 1:
				lineCount += 1
				continue
			line = line.strip()
			splitLine = line.split("\t")
			ensgId = splitLine[3]
			splitEnsgId = ensgId.split('.') #we only keep everything before the dot
			geneName = splitLine[4]
			geneNameConversionMap[splitEnsgId[0]] = geneName
		
	brcaExpression = dict()
	samples = []
	expressionGeneBased = dict()
	addedSamples = dict()
	with open(expressionFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip()
			if lineCount == 0:
				samples += line.split("\t")
	
				lineCount += 1
				continue
			splitLine = line.split("\t")
			
			fullGeneName = splitLine[0]
			if fullGeneName not in geneNameConversionMap:
				continue
			geneName = geneNameConversionMap[fullGeneName] #get the gene name rather than the ENSG ID
	
			if geneName not in expressionGeneBased:
				expressionGeneBased[geneName] = dict()
	
			for col in range(0, len(splitLine)):
				
				if samples[col] == 'gene':
					continue

				if samples[col] not in expressionGeneBased[geneName]:
					expressionGeneBased[geneName][samples[col]] = float(splitLine[col])
					addedSamples[samples[col]] = 0
			
	return expressionGeneBased, list(addedSamples.keys())


def getGtexExpressionData(expressionFile, metadata):
	#Filter for breast tissue samples
	#check which column in the metadata == breast
	tissueId = 'Breast - Mammary Tissue'
	sampleIds = dict()
	with open(metadata, 'r') as inF:
		lineCount = 0
		for line in inF:
			
			if lineCount < 1:
				lineCount += 1
				continue
	
			splitLine = line.split('\t')
			tissue = splitLine[6]

			if tissue == tissueId:
				sampleIds[splitLine[0]] = 0
		
	#Read the GTEx expression data
	
	expressionData = dict()
	samples = []
	addedSamples = dict()
	with open(expressionFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			
			if lineCount < 3:
				lineCount += 1
				
				if lineCount == 3:
					
					splitLine = line.split('\t')
					for col in range(0, len(splitLine)):
						samples.append(splitLine[col])
				
				continue
			
			splitLine = line.split('\t')
			geneName = splitLine[1]
			
			if geneName not in expressionData:
				expressionData[geneName] = dict()
			
			for col in range(0, len(splitLine)):
				
				if samples[col] not in sampleIds:
					continue
	
				if samples[col] == 'Name' or samples[col] == 'Description':
					continue

				if samples[col] not in expressionData[geneName]:
					expressionData[geneName][samples[col]] = float(splitLine[col])
					addedSamples[samples[col]] = 0
	
	return expressionData, list(addedSamples.keys())

brcaExpression, brcaSamples = getBrcaExpression(sys.argv[4])
gtexExpression, gtexSamples = getGtexExpressionData(sys.argv[1], sys.argv[2])

#merge these together
allSamples = ['']

for sample in gtexSamples:
	allSamples.append(sample)
	
for sample in brcaSamples:
	allSamples.append(sample)

print(len(gtexSamples))
print(len(brcaSamples))
exit()
expressionArray = [allSamples]

for gene in gtexExpression:
	
	geneExpr = [gene]
	
	if gene not in brcaExpression:
		continue
	
	for sample in gtexSamples:
		geneExpr.append(gtexExpression[gene][sample])
	for sample in brcaSamples:
		
		if sample == 'gene':
			continue

		geneExpr.append(brcaExpression[gene][sample])
	
	expressionArray.append(geneExpr)
	
expressionArray = np.array(expressionArray, dtype='object')	
	


np.savetxt('/hpc/cog_bioinf/ridder/users/mnieboer/data/pipeline/read_counts/brcaGtexMerged.txt',
		   expressionArray, delimiter='\t', fmt='%s')


