# read in the metadata and check the cancer types
import numpy as np
import glob


metadataFile = '/hpc/cuppen/shared_resources/HMF_data/DR-104/metadata/metadata.tsv'
cancerTypes = dict()
samples = []
with open(metadataFile, 'rb') as inF:
	
	lineCount = 0
	for line in inF:
		line = line.decode('ISO-8859-1')
		if lineCount < 1:
			lineCount += 1
			continue
		
		splitLine = line.split('\t')
		cancerType = splitLine[6]
		
		if cancerType not in cancerTypes:
			cancerTypes[cancerType] = 0
		cancerTypes[cancerType] += 1
		
		samples.append([splitLine[1], cancerType])
		
print(cancerTypes)
print(len(cancerTypes))

samples = np.array(samples, dtype='object')


#check for each sample if there is expression data
matchedExpressionSamples = dict()
rnaFolder = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/isofix/'
 #get all folfders in this main folder
 
folders = glob.glob(rnaFolder + '/*')
for sample in samples:

	if rnaFolder + sample[0] in folders:
		if sample[1] not in matchedExpressionSamples:
			matchedExpressionSamples[sample[1]] = 0

		matchedExpressionSamples[sample[1]] += 1
		
print(matchedExpressionSamples)
print(len(matchedExpressionSamples))
	