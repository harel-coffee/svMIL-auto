"""
	First look in the GTEx data. How is the expression distributed?
	What are the mean and std of all genes in this dataset? Is it approximately the same between the patients, or should it be normalized?
	
	Then, convert the expression data to RPKM to get the same as for the BRCA patients
	Make the same visualizations for the BRCA dataset. Can we easily compare to GTEx?

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

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
	
	samples = np.array(list(addedSamples.keys()), dtype='object')
	#convert back to array
	expressionArray = []
	for geneName in expressionData:
		geneData = [geneName]
		for sample in samples: # make sure to have all samples in the same order
			
			geneData.append(expressionData[geneName][sample])
			
		expressionArray.append(geneData)
	
	expressionArray = np.array(expressionArray, dtype='object')
	#save it to proper numpy readable format, because parsing the file is very slow
	#np.save('../../data/gtex/gtex_breast.npy', expressionData)
	
	dataHeader = '\t'.join([''] + list(samples))
	print(dataHeader)
	np.savetxt('../../data/gtex/gtex_breast_raw.txt', expressionArray, header = dataHeader, fmt='%s', delimiter='\t')
	
	#
	
	

#save the expression data output to disk nd quickly re-load it next run. 
#getGtexExpressionData(sys.argv[1], sys.argv[2])

def plotRawExpression(gtexExpression):
	#make boxplots of a handful of patients
	plotData = []
	samples = 50
	sampleInd = -1
	for sample in gtexExpression:
		sampleInd += 1
		if sampleInd > samples:
			continue
		
		exprValues = np.array(list(gtexExpression[sample].values()))
		exprValues = exprValues[exprValues < np.percentile(exprValues, 75)]
		plotData.append(exprValues)
	
	plt.boxplot(plotData)
	plt.show()
	#plt.savefig('brcaExpression_normalized.svg')

def getGtexExpression(expressionFile):
	#normalize this data. Use R for this. 
		brcaExpression = dict()
		expressionGeneBased = dict()
		samples = ['']
		with open(expressionFile, 'r') as inF:
			lineCount = 0
			for line in inF:
				line = line.strip()
				if lineCount == 0:
					samples += line.split("\t")
		
					lineCount += 1
					continue
				splitLine = line.split("\t")
				geneName = splitLine[0]
				
				if geneName not in expressionGeneBased:
					expressionGeneBased[geneName] = dict()
				
				for col in range(0, len(splitLine)):
					
					if samples[col] == '':
						continue
					
					if samples[col] not in brcaExpression:
						brcaExpression[samples[col]] = dict()
					if geneName not in brcaExpression[samples[col]]:
						brcaExpression[samples[col]][geneName] = float(splitLine[col])
						
					if samples[col] not in expressionGeneBased[geneName]:
						expressionGeneBased[geneName][samples[col]] = float(splitLine[col])
				
		return brcaExpression, expressionGeneBased

gtexExpression, gtexExpressionGeneBased = getGtexExpression('../../data/gtex/gtex_normalized.txt')

#plotRawExpression(gtexExpression)

#check the std of genes in the gtex data. That should not be huge.
# 
# for gene in gtexExpressionGeneBased:
# 	
# 	values = list(gtexExpressionGeneBased[gene].values())
# 	
# 	if np.std(values) > 100:
# 		print(gene)
# 		plt.boxplot(values)
# 		plt.show()
# 		exit()
# 	
# 	print(np.std(values) / np.mean(values))
# 
# exit()

#For each patient and gene, find the TAD that the gene is in. If the TAD is in the disrupted patients list, add the gene expression to the disrupted set.
#Otherwise, add the gene expression to the non-disrutped set.

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
	samples = ['']
	expressionGeneBased = dict()
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
				
				if samples[col] == '':
					continue
				
				if samples[col] not in brcaExpression:
					brcaExpression[samples[col]] = dict()
				if geneName not in brcaExpression[samples[col]]:
					brcaExpression[samples[col]][geneName] = float(splitLine[col])
				if samples[col] not in expressionGeneBased[geneName]:
					expressionGeneBased[geneName][samples[col]] = float(splitLine[col])
				
			
	return brcaExpression, expressionGeneBased	

brcaExpression, expressionGeneBased = getBrcaExpression(sys.argv[4])

for gene in gtexExpressionGeneBased:
	
	values = list(gtexExpressionGeneBased[gene].values())
	
	#if np.std(values) > 100:
	#	print(gene)
	#	plt.boxplot(values)
	#	plt.show()
	#	exit()
	
	print(np.std(values))

exit()

#plot this data
plotRawExpression(brcaExpression)







