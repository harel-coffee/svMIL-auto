"""
	For every TAD boundary disruption by a non-coding SV, check which genes are in the same TAD that are potentially affected.
	For each affected gene, also record if there is a coding SV disrupting this gene to get to a list of non-coding effects only. How does this differ from the set of nc-only genes we find with our tool?
	How many of the genes found with this simplistic approach are DEG? 

"""

from __future__ import absolute_import
from __future__ import print_function
from inputParser import InputParser
import settings
import sys
import numpy as np
import re
from scipy import stats
from six.moves import range

#1. Get the TADs that are disrupted by non-coding SVs. 

somaticSVs = InputParser().getSVsFromFile(sys.argv[1], "all")
tads = InputParser().getTADsFromFile(sys.argv[2])
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#Remove translocations, these do not disrupt boundaries
filteredSomaticSVs = []
for sv in somaticSVs:
	if sv[0] == sv[3]:
		filteredSomaticSVs.append(sv)
filteredSomaticSVs = np.array(filteredSomaticSVs, dtype="object")


#For each TAD, determine which SVs start or end within the TAD, and cover the boundary
affectedGenes = []
nonCodingSamples = dict()
for tad in tads:
	
	#match on chromosome
	svSubset = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
	geneChrSubset = genes[genes[:,0] == tad[0]]
	
	#if the start of the SV is within the TAD, and the end is outside of the TAD
	startingSVs = svSubset[(svSubset[:,1] >= tad[1]) * (svSubset[:,1] <= tad[2]) * (svSubset[:,5] > tad[2])]
	
	#if the end of the SV is within the TAD, and the start is outside of the TAD
	endingSVs = svSubset[(svSubset[:,5] >= tad[1]) * (svSubset[:,5] <= tad[2]) * (svSubset[:,1] < tad[1])]
	
	#Which genes are in these TADs? For the starting, take all genes on the left
	# for sv in startingSVs:
	# 	
	# 	matchingStartGenes = (geneChrSubset[:,1] >= tad[1]) * (geneChrSubset[:,1] <= sv[1])
	# 	matchingEndGenes = (geneChrSubset[:,2] >= tad[1]) * (geneChrSubset[:,2] <= sv[1])
	# 
	# 	matchingGenes = geneChrSubset[matchingStartGenes * matchingEndGenes]
	# 	for gene in matchingGenes:
	# 		
	# 		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
	# 		affectedGenes.append(gene[3].name + "_" + svStr)
	# 		if gene[3].name not in nonCodingSamples:
	# 			nonCodingSamples[gene[3].name] = []
	# 		nonCodingSamples[gene[3].name].append(sv[7])
	# 	
	# #For the ending, take all genes on the right
	# for sv in endingSVs:
	# 	
	# 	matchingStartGenes = (geneChrSubset[:,1] >= sv[5]) * (geneChrSubset[:,1] <= tad[2])
	# 	matchingEndGenes = (geneChrSubset[:,2] >= sv[5]) * (geneChrSubset[:,2] <= tad[2])
	# 
	# 	matchingGenes = geneChrSubset[matchingStartGenes * matchingEndGenes]
	# 	for gene in matchingGenes:
	# 		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
	# 		affectedGenes.append(gene[3].name + "_" + svStr)
	# 		if gene[3].name not in nonCodingSamples:
	# 			nonCodingSamples[gene[3].name] = []
	# 		nonCodingSamples[gene[3].name].append(sv[7])

	if len(startingSVs) > 0 or len(endingSVs) > 0:
		matchingStartGenes = (geneChrSubset[:,1] >= tad[1]) * (geneChrSubset[:,1] <= tad[2])
		matchingEndGenes = (geneChrSubset[:,2] >= tad[1]) * (geneChrSubset[:,2] <= tad[2])
	
		matchingGenes = geneChrSubset[matchingStartGenes * matchingEndGenes]
		for gene in matchingGenes:
			svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
			affectedGenes.append(gene[3].name + "_" + svStr)
			if gene[3].name not in nonCodingSamples:
				nonCodingSamples[gene[3].name] = []
			nonCodingSamples[gene[3].name].append(sv[7])

print(affectedGenes)
	
#Find out which one of these is DEG
#Compare to samples without ANY SV at all

#First determine all samples that have a coding SV
codingSamples = dict()
for gene in genes:
	
	svSubset = filteredSomaticSVs[filteredSomaticSVs[:,0] == gene[0]]
	
	#Find SVs that overlap the gene
	matchingSVs = svSubset[(svSubset[:,5] > gene[1]) * (svSubset[:,1] < gene[2])]
	for sv in matchingSVs:
		
		if gene[3].name not in codingSamples:
			codingSamples[gene[3].name] = []
		codingSamples[gene[3].name].append(sv[7])


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
print(expressionData)
	

geneSampleRef = dict()
for gene in nonCodingSamples:
	geneSampleRef[gene] = nonCodingSamples[gene]
for gene in codingSamples:
	if gene not in geneSampleRef:
		geneSampleRef[gene] = codingSamples[gene]
	else:
		geneSampleRef[gene] += codingSamples[gene]


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
print("done getting expr for samples")

#Also set the negative set for every gene consisting of the expression of all samples wthout any SV
negativeExpr = dict()
for gene in geneSampleExpr:
	matchedFullSampleNames = list(geneSampleExpr[gene].keys())
	
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
print("negative expr done")


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
		matchedFullSampleNames = list(geneSampleExpr[gene].keys())
					
		
		negativeSampleExpressionValues = negativeExpr[gene]
		
		#Get the expression z-score for this pair
		if np.std(negativeSampleExpressionValues) == 0:
			continue
	
		z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
		pValue = stats.norm.sf(abs(z))*2
	
		perPairDifferentialExpression[pair] = pValue
		
	return perPairDifferentialExpression

#Get the p-value for each pair in coding & non-coding
perPairDifferentialExpression = getDEPairs(affectedGenes, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)	
print(perPairDifferentialExpression)

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.savetxt('naiveTadDisr_nonCodingDEGs_all.txt', perPairDifferentialExpressionArrayFiltered, delimiter='\t', fmt='%s')





