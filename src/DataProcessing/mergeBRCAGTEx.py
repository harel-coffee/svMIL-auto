"""
	Merge the BRCA and GTEx expression data into 1 big table.
	Make sure that the genes are matched. 

"""

import sys
import numpy as np


#get the gene expression
def getBrcaExpression(expressionFile):
		
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
			
			geneId = splitLine[1]

			if geneId not in expressionGeneBased:
				expressionGeneBased[geneId] = dict()
	
			for col in range(0, len(splitLine)):
				
				if samples[col] in ['Name', 'Description']:
					continue

				if samples[col] not in expressionGeneBased[geneId]:
					expressionGeneBased[geneId][samples[col]] = float(splitLine[col])
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

brcaExpression, brcaSamples = getBrcaExpression(sys.argv[3])
gtexExpression, gtexSamples = getGtexExpressionData(sys.argv[1], sys.argv[2])

#merge these together
allSamples = ['']

for sample in gtexSamples:
	allSamples.append(sample)
	
for sample in brcaSamples:
	allSamples.append(sample)

print(len(gtexSamples))
print(len(brcaSamples))

expressionArray = [allSamples]

for gene in gtexExpression: #gtex has fewer genes (about 4) so this works best
	
	geneExpr = [gene]
	
	if gene not in brcaExpression:
		continue
	
	for sample in gtexSamples:
		geneExpr.append(gtexExpression[gene][sample])
	for sample in brcaSamples:
		
		if sample in ['Name', 'Description']:
			continue

		geneExpr.append(brcaExpression[gene][sample])
	
	expressionArray.append(geneExpr)
	
expressionArray = np.array(expressionArray, dtype='object')	
	
print(expressionArray.shape)

np.savetxt('../data/expression/brcaGtex_merged.txt',
		   expressionArray, delimiter='\t', fmt='%s')


