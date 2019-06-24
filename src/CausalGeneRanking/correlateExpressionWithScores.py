
import sys
import re
import numpy as np
from scipy import stats


#1. Compute the differential expression matrix per sample
#For every gene with at least 1 sample, make a null distribution of the expression values of the samples that do not have an SV affecting this gene
#Compute a p-value for that sample

# 
# geneScoreFile = sys.argv[1]
# geneScores = np.loadtxt(geneScoreFile, dtype="object")
# 
# expressionFile = sys.argv[2]
# 
# expressionData = []
# samples = []
# with open(expressionFile, 'r') as inF:
# 	lineCount = 0
# 	for line in inF:
# 		line = line.strip()
# 		if lineCount == 0:
# 			samples = line.split("\t")
# 			lineCount += 1
# 			continue
# 		if lineCount < 2:
# 			lineCount += 1
# 			continue
# 		splitLine = line.split("\t")
# 		fullGeneName = splitLine[0]
# 		geneName = fullGeneName.split("|")[0]
# 
# 		data = splitLine[1:len(splitLine)-1]
# 		fixedData = [geneName]
# 		fixedData += data
# 		expressionData.append(fixedData)
# 
# expressionData = np.array(expressionData, dtype="object")
# 
# featureMatrixFile = sys.argv[3]
# featureMatrix = np.loadtxt(featureMatrixFile, dtype="object")
# 
# pValueMatrix = np.zeros(featureMatrix.shape, dtype="object")
# 
# for gene in geneScores:
# 	
# 	#Get the index of the gene in the feature matrix
# 	geneInd = np.where(featureMatrix[:,0] == gene[0])[0]
# 	featureMatrix[geneInd,1] = gene[0]
# 	
# 	#1. get the samples affecting this gene
# 	if gene[0] not in expressionData[:,0]:
# 		continue
# 	geneExpression = expressionData[expressionData[:,0] == gene[0]][0]
# 	
# 	geneSamples = gene[31].split(",")
# 	if geneSamples[0] == "None":
# 		pValueMatrix[geneInd,1:] = 100
# 		continue
# 	matchedFullSampleNames = []
# 	#2. get the expression of the negative distribution
# 	for geneSample in geneSamples:
# 
# 		shortSampleName = geneSample.split("brca")[1]
# 		
# 		#match the sample name with the expression sample name
# 		for sampleInd in range(0, len(samples)):
# 			sample = samples[sampleInd]
# 			if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
# 				matchedFullSampleNames.append(sample) #keep this to check later for the negative set
# 				
# 	#Get 5 random samples that are not affecting this gene
# 	unmatchedSamples = np.setdiff1d(samples[1:len(samples)-1], matchedFullSampleNames) #exclude hybrid ref
# 	negativeSamples = []
# 	for sample in unmatchedSamples: #sample tumor samples, exclude normals
# 		splitSample = sample.split("-")
# 		code = int(splitSample[len(splitSample)-1])
# 		
# 		if code < 10: 
# 			negativeSamples.append(sample)
# 	
# 	#Get the expression of these samples
# 	negativeSampleExpressionValues = []
# 	for sample in negativeSamples:
# 		sampleInd = samples.index(sample)				
# 		negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
# 	
# 	negativeDistributionMean = np.mean(negativeSampleExpressionValues)
# 	
# 	#3. For each sample, compute the expression 
# 	for geneSample in geneSamples:
# 		shortSampleName = geneSample.split("brca")[1]
# 		
# 		
# 		#match the sample name with the expression sample name
# 		for sampleInd in range(0, len(samples)):
# 			sample = samples[sampleInd]
# 			if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
# 				
# 				#Get the last 2 numbers
# 				splitSample = sample.split("-")
# 				code = int(splitSample[len(splitSample)-1])
# 				
# 				if code < 10: #above 9 are the normal samples, which we do not want to include here
# 					sampleInd = samples.index(sample)
# 					
# 					sampleExpression = float(geneExpression[sampleInd])
# 		
# 					#Compute p-value
# 					#print sampleExpression
# 					#print negativeDistributionMean
# 					pValue = stats.ttest_1samp(negativeSampleExpressionValues, sampleExpression)[1]
# 					
# 					matrixSampleInd = np.where(featureMatrix[0,:] == geneSample)[0]
# 					pValueMatrix[0,matrixSampleInd] = geneSample
# 					pValueMatrix[geneInd,matrixSampleInd] = pValue
# 		
# np.savetxt("pValues.txt", pValueMatrix, delimiter='\t', fmt='%s')		

#Correlate the pvalue matrix with the negative feature matrix

#This needs to be a separate script

featureMatrixFile = sys.argv[1]
featureMatrix = np.loadtxt(featureMatrixFile, dtype="object")

pValueFile = sys.argv[2]
pValueMatrix = np.loadtxt(pValueFile, dtype="object")

#filter out genes without svs
filteredFeatureMatrix = []
filteredPValueMatrix = []
for row in range(1, featureMatrix.shape[0]):

	if np.sum(featureMatrix[row,:]) < 1:
		continue
	filteredFeatureMatrix.append(featureMatrix[row,1:])
	pValues = []
	for col in range(1,featureMatrix.shape[1]):
		if pValueMatrix[row,col] < 0.05:
			pValues.append(1)
		else:
			pValues.append(-1)
		#filteredPValueMatrix.append(pValues[row,1:])
	filteredPValueMatrix.append(pValues)

filteredFeatureMatrix = np.array(filteredFeatureMatrix, dtype="float") + 0.0001
filteredPValueMatrix = np.array(filteredPValueMatrix,dtype="float") + 0.0001

print filteredFeatureMatrix
print filteredPValueMatrix

for row in range(0, filteredFeatureMatrix.shape[0]):
	
	#Check for this gene how many samples have a positive correlation
	
	print filteredFeatureMatrix[row,:]
	print np.corrcoef(filteredFeatureMatrix[row,:], filteredPValueMatrix[row,:])[0,1]
	



#out = np.corrcoef(filteredFeatureMatrix,filteredPValueMatrix)
out = np.dot(filteredFeatureMatrix,filteredPValueMatrix.T)
print out



# 
# def correlation_matrix(df):
#     from matplotlib import pyplot as plt
#     from matplotlib import cm as cm
# 
#     fig = plt.figure()
#     ax1 = fig.add_subplot(111)
#     cmap = cm.get_cmap('jet', 30)
#     cax = ax1.imshow(df, interpolation="nearest", cmap=cmap)
#     ax1.grid(True)
#     #plt.title('Abalone Feature Correlation')
#     #labels=['Sex','Length','Diam','Height','Whole','Shucked','Viscera','Shell','Rings',]
#     #ax1.set_xticklabels(labels,fontsize=6)
#     #ax1.set_yticklabels(labels,fontsize=6)
#     # Add colorbar, make sure to specify tick locations to match desired ticklabels
#     #fig.colorbar(cax, ticks=[.75,.8,.85,.90,.95,1])
#     plt.show()
# 
# correlation_matrix(out)

