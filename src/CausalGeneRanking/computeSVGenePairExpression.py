
import sys
import numpy as np
import re
from scipy import stats

# Get the expression z-scores for every SV. 

nonCodingPairs = np.loadtxt(sys.argv[1], dtype="object")
codingPairs = np.loadtxt(sys.argv[2], dtype="object")

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

#np.save('codingNonCodingPairDEGs.npy', perPairDifferentialExpressionArrayFiltered)

#instead of writing the whole array to a file, write per SV which DEGs these are linked to
svsAndDegs = dict()
for pair in perPairDifferentialExpressionArrayFiltered[:,0]:
	splitPair = pair.split("_")
	svEntries = splitPair[1:]
	sv = "_".join(svEntries)
	
	if sv not in svsAndDegs:
		svsAndDegs[sv] = 0
	svsAndDegs[sv] += 1

pairs = np.empty([len(svsAndDegs), 2], dtype="object")	
for svInd in range(0, len(svsAndDegs)):
	pairs[svInd, 0] = svsAndDegs.keys()[svInd]
	pairs[svInd,1] = svsAndDegs[svsAndDegs.keys()[svInd]]

np.savetxt(sys.argv[1] + "_degPairs.txt", pairs, delimiter="\t", fmt="%s")

# # Output DEG pairs for non-coding only
# perPairDifferentialExpression = getDEPairs(nonCodingPairs[:,0], geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
# print "done"
# 
# perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
# perPairDifferentialExpressionArray[:,0] = perPairDifferentialExpression.keys()
# perPairDifferentialExpressionArray[:,1] = perPairDifferentialExpression.values()
# 
# from statsmodels.sandbox.stats.multicomp import multipletests
# reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')
# 
# perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]
# 
# np.save('nonCodingPairDEGs.npy', perPairDifferentialExpressionArrayFiltered)
# 
