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
			
			for col in range(0, len(splitLine)):
				
				if samples[col] not in sampleIds:
					continue
	
				if samples[col] == 'Name' or samples[col] == 'Description':
					continue
	
				if samples[col] not in expressionData:
					expressionData[samples[col]] = dict()
				
				if geneName not in expressionData[samples[col]]:
					expressionData[samples[col]][geneName] = float(splitLine[col])
	
	
	#save it to proper numpy readable format, because parsing the file is very slow
	np.save('../../data/gtex/gtex_breast.npy', expressionData)

#save the expression data output to disk nd quickly re-load it next run. 
#getGtexExpressionData(sys.argv[1], sys.argv[2])
gtexExpression = np.load('../../data/gtex/gtex_breast.npy', allow_pickle=True, encoding='latin1').item()

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

#plotRawExpression(gtexExpression)

#Convert to RPKM to compare to the BRCA dataset

#per sample, get the total number of reads. Divide this by 1 million
#get the read counts and divide it by the above number
#divide this value by the length of the gene in kb


#make the same plots for the BRCA dataset




