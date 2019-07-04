"""	
	Input pairs of non-coding SVs and genes and coding SVs and genes.
	For each SV, compute how many genes it affects in the non-coding way, and how many in the coding way. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import stats

nonCodingPairs = np.loadtxt(sys.argv[1], dtype="object")
codingPairs = np.loadtxt(sys.argv[2], dtype="object")

#For every SV in the non-coding set, count how many genes are affected by this SV (frequency in dataset)
nonCodingSVCounts = dict()
for pair in nonCodingPairs:
	
	splitPair = pair[0].split("_")
	svEntries = splitPair[1:]
	sv = "_".join(svEntries)
	
	if sv not in nonCodingSVCounts:
		nonCodingSVCounts[sv] = 0
	nonCodingSVCounts[sv] += 1
	
print nonCodingSVCounts	

#Repeat for the coding pairs
codingSVCounts = dict()
for pair in codingPairs:
	
	splitPair = pair.split("_")
	svEntries = splitPair[1:]
	sv = "_".join(svEntries)
	
	if sv not in codingSVCounts:
		codingSVCounts[sv] = 0
	codingSVCounts[sv] += 1
	
print codingSVCounts

#Plot counts as scatter
#Make x and y axis with values in the same entries.
#Append additional SVs that are unique to coding/noncoding

totalSVs = np.union1d(nonCodingSVCounts.keys(), codingSVCounts.keys())
print len(totalSVs)

svEffects = np.empty([len(totalSVs), 4], dtype="object")
svEffects[:,1] = 0
svEffects[:,2] = 0
ind = 0
for sv in nonCodingSVCounts:
	svEffects[ind,0] = sv
	svEffects[ind,1] = nonCodingSVCounts[sv]
	ind += 1
	
	if nonCodingSVCounts[sv] > 80:
		print "non-coding guy: ", sv
	
for sv in codingSVCounts:
	
	if sv in svEffects[:,0]:
		rowInd = svEffects[:,0] == sv
		svEffects[rowInd,2] = codingSVCounts[sv]
	else:
		svEffects[ind,0] = sv
		svEffects[ind,2] = codingSVCounts[sv]
		ind += 1
	if codingSVCounts[sv] > 1500:
		print "wow I'm huge:", sv


print svEffects


#Get the expression z-scores for every SV. 

expressionFile = sys.argv[3]

expressionData = []
samples = []
with open(expressionFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		if lineCount == 0:
			samples = line.split("\t")
			lineCount += 1
			continue
		if lineCount < 2:
			lineCount += 1
			continue
		splitLine = line.split("\t")
		fullGeneName = splitLine[0]
		geneName = fullGeneName.split("|")[0]

		data = splitLine[1:len(splitLine)-1] 
		fixedData = [geneName]
		fixedData += data
		expressionData.append(fixedData)

expressionData = np.array(expressionData, dtype="object")	
print expressionData

#Get the z-scores for every pair

#For every gene, get a list of all samples in which this gene is affected to exclude these and make a null distribution
geneSampleRef = dict()
for pair in nonCodingPairs[:,0]:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	
	if gene not in geneSampleRef:
		geneSampleRef[gene] = []
	geneSampleRef[gene].append(sample)

for pair in codingPairs:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[len(splitPair)-1]
	
	if gene not in geneSampleRef:
		geneSampleRef[gene] = []
	geneSampleRef[gene].append(sample)

#Set for every gene the expression values in all possible samples for lookup
geneSampleExpr = dict()
for gene in geneSampleRef:
	
	if gene not in expressionData[:,0]:
		continue
	
	geneSamples = geneSampleRef[gene]
	geneSampleExpr[gene] = dict()
	geneExpression = expressionData[expressionData[:,0] == gene][0]
	for geneSample in geneSamples:
		
		shortSampleName = geneSample.split("brca")[1]
		
		#match the sample name with the expression sample name
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
				
				splitSample = sample.split("-")
				code = int(splitSample[len(splitSample)-1])
				
				if code < 10: #above 9 are the normal samples, which we do not want to include here
					sampleInd = samples.index(sample)
					
					geneSampleExpr[gene][geneSample] = float(geneExpression[sampleInd])
print "done getting expr for samples"

#Also set the negative set for every gene consisting of the expression of all samples wthout any SV
negativeExpr = dict()
for gene in geneSampleExpr:
	matchedFullSampleNames = geneSampleExpr[gene].keys()
	
	#Get all the samples without an SV for this gene
	unmatchedSamples = np.setdiff1d(samples[1:len(samples)-1], matchedFullSampleNames) #exclude hybrid ref
	negativeSamples = []
	for sample in unmatchedSamples: #sample tumor samples, exclude normals
		splitSample = sample.split("-")
		code = int(splitSample[len(splitSample)-1])
		
		if code < 10: 
			negativeSamples.append(sample)
		
	#Get the expression of these samples
	negativeSampleExpressionValues = []
	for sample in negativeSamples:
		sampleInd = samples.index(sample)				
		negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
	
	negativeExpr[gene] = negativeSampleExpressionValues
print "negative expr done"

def getDEPairs(pairs, geneSampleRef, epressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr):
									
	for pair in pairs:
		splitPair = pair.split("_")
		gene = splitPair[0]
		pairSample = splitPair[len(splitPair)-1]
		shortPairSampleName = pairSample.split("brca")[1]
		sv = "_".join(splitPair[1:])
		if gene not in expressionData[:,0]:
			continue
		
		
		sampleExpressionValue = geneSampleExpr[gene][pairSample] #expression values of this gene in all samples
		matchedFullSampleNames = geneSampleExpr[gene].keys()
					
		
		negativeSampleExpressionValues = negativeExpr[gene]
		
		#Get the expression z-score for this pair
		if np.std(negativeSampleExpressionValues) == 0:
			continue
	
		z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
		pValue = stats.norm.sf(abs(z))*2
	
		perPairDifferentialExpression[pair] = pValue
		
	return perPairDifferentialExpression

#Get the p-value for each pair in coding & non-coding
perPairDifferentialExpression = getDEPairs(nonCodingPairs[:,0], geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print "done"
perPairDifferentialExpression = getDEPairs(codingPairs, geneSampleRef, expressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr)
print "coding done"
#Do multiple testing correction
#print perPairDifferentialExpression

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = perPairDifferentialExpression.keys()
perPairDifferentialExpressionArray[:,1] = perPairDifferentialExpression.values()


from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.save('codingNonCodingPairDEGs_shuffled.npy', perPairDifferentialExpressionArrayFiltered)
# exit()

#perPairDifferentialExpressionArrayFiltered = np.load('codingNonCodingPairDEGs.npy')
print perPairDifferentialExpressionArrayFiltered.shape

#For every pair, assign a +1 to the SV if it has a DEG gene
svEffects[:,3] = 0
for pairInd in range(0, svEffects.shape[0]):
	pair = svEffects[pairInd,0]
	for degPair in perPairDifferentialExpressionArrayFiltered[:,0]:
		splitDegPair = degPair.split("_")
		sv = "_".join(splitDegPair[1:])
		
		if sv == pair:
			svEffects[pairInd,3] += 1
		
	#Find how often this SV is linked to a DEG gene
	
#Plot
plt.scatter(svEffects[:,2], svEffects[:,1], c=svEffects[:,3]) #c=colors
plt.colorbar()
plt.show()


