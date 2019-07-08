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
# 
# #Get the properties of SVs beforehand to match these later
# svProperties = []
# 
# for sv in svEffects[:,0]:
# 	
# 	splitSV = sv.split("_")
# 	chrom = splitSV[0]
# 	size = int(splitSV[5]) - int(splitSV[1])
# 	sample = splitSV[6]
# 	
# 	
# 	
# 	svProperties.append([sv, chrom + "_" + str(size) + "_" + sample])
# 
# svProperties = np.array(svProperties, dtype="object")
# print svProperties
# 
# #Go through the shuffled files coding & non-coding and compute how often the real noncoding/coding is higher for that SV
# 
# import glob
# 
# #First determine non-coding potential
# shuffledNCFiles = glob.glob(sys.argv[3] + "/geneSV*")
# 
# #For each file, 
# #Find which SV it is in SV effects
# #See if the score of that SV is higher or not
# #Keep the total number of times that it is higher in the shuffled case
# #Compute a p-value per SV and mark the original SVs in the plot
# 
# nonCodingEnrichment = dict() #For every SV, set if the nc/coding is higher in the random case than in the true case
# nonCodingEnrichmentTrueSV = dict()
# fileCount = 0
# for shuffledFile in shuffledNCFiles:
# 	print "shuffled file: ", fileCount
# 
# 	if shuffledFile == ".DS_Store":
# 		continue
# 	
# 	nonCodingPairs = np.loadtxt(shuffledFile, dtype="object")
# 	
# 	nonCodingSVCounts = dict()
# 	for pair in nonCodingPairs:
# 		
# 		splitPair = pair[0].split("_")
# 		svEntries = splitPair[1:]
# 		sv = "_".join(svEntries)
# 		
# 		if sv not in nonCodingSVCounts:
# 			nonCodingSVCounts[sv] = 0
# 		nonCodingSVCounts[sv] += 1
# 
# 	#Get the corresponding coding file
# 	splitFileName = shuffledFile.split("_")
# 	permutationId = splitFileName[len(splitFileName)-1]
# 	
# 	codingPairsFile = glob.glob(sys.argv[3] + "/geneCoding*_" + permutationId)[0]
# 	
# 	codingPairs = np.loadtxt(codingPairsFile, dtype="object")
# 	
# 	codingSVCounts = dict()
# 	for pair in codingPairs:
# 		
# 		splitPair = pair.split("_")
# 		svEntries = splitPair[1:]
# 		sv = "_".join(svEntries)
# 		
# 		if sv not in codingSVCounts:
# 			codingSVCounts[sv] = 0
# 		codingSVCounts[sv] += 1
# 		
# 	allSVs = np.union1d(nonCodingSVCounts.keys(), codingSVCounts.keys())
# 	for sv in allSVs:
# 		
# 		ncPotential = 0
# 		if sv in nonCodingSVCounts and sv in codingSVCounts:
# 			
# 			ncPotential = nonCodingSVCounts[sv] / float(codingSVCounts[sv])
# 		
# 		#Match to the true SV
# 		splitSV = sv.split("_")
# 		svSize = int(splitSV[5]) - int(splitSV[1])
# 		svMatchCode = splitSV[0] + "_" + str(svSize) + "_" + splitSV[6]
# 		
# 		trueSV = svProperties[svProperties[:,1] == svMatchCode]
# 		
# 		if len(trueSV) > 0: #SV match found
# 			
# 			#Is the nc/coding higher in the true case than in this random case?
# 			trueSVEffects = svEffects[svEffects[:,0] == trueSV[0][0]][0]
# 			
# 			ncPotentialTrueSV = 0
# 			if trueSVEffects[2] > 0:
# 				ncPotentialTrueSV = trueSVEffects[1] / float(trueSVEffects[2])
# 			enrichment = False
# 			if ncPotential > ncPotentialTrueSV:
# 				enrichment = True
# 			
# 			if trueSVEffects[0] not in nonCodingEnrichment:
# 				nonCodingEnrichment[trueSVEffects[0]] = []
# 			nonCodingEnrichment[trueSVEffects[0]].append(ncPotential)
# 			
# 			if trueSVEffects[0] not in nonCodingEnrichmentTrueSV:
# 				nonCodingEnrichmentTrueSV[trueSVEffects[0]] = ncPotentialTrueSV
# 	
# 	fileCount += 1		
# 		
# #for every SV, compute the p-value
# svSignificance = []
# for sv in nonCodingEnrichment:
# 	#trueCount = [int(i) for i in nonCodingEnrichment[sv]]
# 	#trueCount = sum(trueCount)
# 	
# 	if np.std(nonCodingEnrichment[sv]) == 0:
# 		continue
# 	
# 	z = (nonCodingEnrichmentTrueSV[sv] - np.mean(nonCodingEnrichment[sv])) / float(np.std(nonCodingEnrichment[sv]))	
# 	pValue = stats.norm.sf(abs(z))*2
# 	
# 	#proportion = (trueCount + 1) / float(len(nonCodingEnrichment[sv]) + 1)
# 	svSignificance.append([sv, pValue])
# 	
# svSignificance = np.array(svSignificance, dtype="object")
# 
# np.savetxt('Output/significantNCProportion.txt', svSignificance, fmt='%s', delimiter='\t')
# 
# #Do multiple testing correction
# 
# from statsmodels.sandbox.stats.multicomp import multipletests
# reject, pAdjusted, _, _ = multipletests(svSignificance[:,1], method='bonferroni')
# 
# svSignificanceCorrected = []
# for svInd in range(0, svSignificance.shape[0]):
# 	
# 	if reject[svInd] == True:
# 		svSignificanceCorrected.append([svSignificance[svInd,0], pAdjusted[svInd]])
# 	
# 
# svSignificanceCorrected = np.array(svSignificanceCorrected, dtype="object")	
# 	
# np.savetxt('Output/significantNCProportion_multipleTestCorrected.txt', svSignificanceCorrected, fmt='%s', delimiter='\t')
# 
# perPairDifferentialExpressionArrayFiltered = np.load('codingNonCodingPairDEGs.npy')
# print perPairDifferentialExpressionArrayFiltered.shape
# 
# #For every pair, assign a +1 to the SV if it has a DEG gene
# svEffects[:,3] = 0
# for pairInd in range(0, svEffects.shape[0]):
# 	pair = svEffects[pairInd,0]
# 	for degPair in perPairDifferentialExpressionArrayFiltered[:,0]:
# 		splitDegPair = degPair.split("_")
# 		sv = "_".join(splitDegPair[1:])
# 		
# 		if sv == pair:
# 			svEffects[pairInd,3] += 1
# plt.scatter(svEffects[:,2], svEffects[:,1], c=svEffects[:,3]) #c=colors
# 
#Do overlay because the colormap is not working separately
svSignificanceCorrected = np.loadtxt('Output/significantNCProportion_multipleTestCorrected.txt', dtype="object")
print "Number of significant gene SV pairs: ", svSignificanceCorrected.shape

degPairs = np.load(sys.argv[4])

#Add in the plot which SVs are significant
degCount = 0
signNCDegPairs = dict()
degGenes = []
degPairsSign = []
degPairsFull = []
for svInd in range(0, svEffects.shape[0]):
	sv = svEffects[svInd,0]
	
	if sv in svSignificanceCorrected[:,0]:
		if sv not in signNCDegPairs:
			signNCDegPairs[sv] = []
		
		#This SV is linked to more DEGs than by random chance
		#Find which genes it is linked to
		# genes = []
		# for pair in codingPairs:
		# 	splitPair = pair.split("_")
		# 	svEntries = splitPair[1:]
		# 	if "_".join(svEntries) == sv:
		# 		genes.append(splitPair[0])
		
		genes = []
		for pair in nonCodingPairs:
			splitPair = pair[0].split("_")
			svEntries = splitPair[1:]
			if "_".join(svEntries) == sv:
				genes.append(splitPair[0])
		
		#print "sv: ", sv
		#print "genes: ", genes
		
		splitSV = sv.split("_")
		#Find which genes are DEG in that sample
		for gene in genes:
			pair = gene + "_" + sv
			if pair in degPairs[:,0]:
				print "deg pair: ", pair
				degCount += 1
				signNCDegPairs[sv].append(pair)
				degGenes.append(gene)
				degPairsSign.append(gene + "_" + splitSV[len(splitSV)-1])
				degPairsFull.append(gene + "_" + sv)
				
				# if gene == "CDH1":
				# 	print "CDH1: ", sv
				# if gene == "ERBB2":
				# 	print "ERBB2: ", sv
				# if gene == "COL1A1":
				# 	print "COL1A1: ", sv
				# if gene == "SPOP":
				# 	print "SPOP: ", sv
				# if gene == "H3F3B":
				# 	print "H3F3B: ", sv
				
				#coding
				if gene == "CDH1":
					print "CDH1: ", sv
				if gene == "RNF43":
					print "RNF43: ", sv
				if gene == "CLTC":
					print "CLTC: ", sv
				if gene == "DDX5":
					print "DDX5: ", sv
				if gene == "PRKAR1A":
					print "PRKAR1A: ", sv
				if gene == "CDK12":
					print "CDK12: ", sv
				
				
print "Number of significant nc potential that are also linked to DEG genes: ", degCount

#np.savetxt("Output/ncPotentialDEGGenesCoding.txt", np.array(degGenes, dtype="object"), delimiter="\t", fmt="%s")

#Check for the DEG pairs which is also found in the 'naive' method

#Show some stats on how many of the nc potential SVs are also linked to DEG genes, in how many different patients etc
samples = []
svWithGenes = 0
geneCounts = []
for signNCDegPair in signNCDegPairs:
	splitPair = signNCDegPair.split("_")
	if splitPair[len(splitPair)-1] not in samples:
		samples.append(splitPair[len(splitPair)-1])

	if len(signNCDegPairs[signNCDegPair]) > 0:
		svWithGenes += 1
	geneCounts.append(len(signNCDegPairs[signNCDegPair]))
	
	if len(signNCDegPairs[signNCDegPair]) > 15:
		print signNCDegPair

print "no of samples: ", len(samples), samples
print "no of svs with at least 1 gene: ", svWithGenes
print "gene counts per sv: ", geneCounts

#Determine what their types are by searching through the original file
from inputParser import InputParser

sizes = []
somaticSVs = InputParser().getSVsFromFile(sys.argv[5], "all")
for signNCDegPair in signNCDegPairs:
	
	for sv in somaticSVs:
		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		if signNCDegPair == svStr:
			print sv
			print "type: ", sv[8].svType
			print sv[5] - sv[1]
			sizes.append(sv[5] - sv[1])

print sizes
print min(sizes)
print max(sizes)


#Which of the significant pairs are also found in the naive way? any that is not found in that way?
naiveDegPairs = np.loadtxt('naiveTadDisr_nonCodingDEGs.txt', dtype="object")

foundCount = 0
notFoundCount = 0
for degPair in degPairsFull:
	if degPair in naiveDegPairs[:,0]:
		print "pair found: ", degPair
		foundCount += 1
	else:
		print "pair not found: ", degPair
		notFoundCount += 1
		for sv in somaticSVs:
			svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
			splitPair = degPair.split("_")
			if svStr == "_".join(splitPair[1:]):
				print sv[8].svType
		
print foundCount
print notFoundCount



exit()

# 
# exit()
# 		
# 		#plt.scatter(svEffects[svInd,2], svEffects[svInd,1], marker="*", color='red') #plot the significant SVs
# 	
# #plt.show()
# 
# 
# exit()




# Get the expression z-scores for every SV. 

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

# #Get the p-value for each pair in coding & non-coding
# perPairDifferentialExpression = getDEPairs(nonCodingPairs[:,0], geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
# print "done"
# perPairDifferentialExpression = getDEPairs(codingPairs, geneSampleRef, expressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr)
# print "coding done"
# #Do multiple testing correction
# #print perPairDifferentialExpression
# 
# perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
# perPairDifferentialExpressionArray[:,0] = perPairDifferentialExpression.keys()
# perPairDifferentialExpressionArray[:,1] = perPairDifferentialExpression.values()
# 
# 
# from statsmodels.sandbox.stats.multicomp import multipletests
# reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')
# 
# perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]
# 
# np.save('codingNonCodingPairDEGs_shuffled.npy', perPairDifferentialExpressionArrayFiltered)

# Output DEG pairs for non-coding only
perPairDifferentialExpression = getDEPairs(nonCodingPairs[:,0], geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print "done"

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = perPairDifferentialExpression.keys()
perPairDifferentialExpressionArray[:,1] = perPairDifferentialExpression.values()

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.save('nonCodingPairDEGs.npy', perPairDifferentialExpressionArrayFiltered)
exit()


# # exit()
# 
# #perPairDifferentialExpressionArrayFiltered = np.load('codingNonCodingPairDEGs.npy')
# print perPairDifferentialExpressionArrayFiltered.shape
# 
# #For every pair, assign a +1 to the SV if it has a DEG gene
# svEffects[:,3] = 0
# for pairInd in range(0, svEffects.shape[0]):
# 	pair = svEffects[pairInd,0]
# 	for degPair in perPairDifferentialExpressionArrayFiltered[:,0]:
# 		splitDegPair = degPair.split("_")
# 		sv = "_".join(splitDegPair[1:])
# 		
# 		if sv == pair:
# 			svEffects[pairInd,3] += 1
# 		
# 	#Find how often this SV is linked to a DEG gene
# 	
# #Plot
# plt.scatter(svEffects[:,2], svEffects[:,1], c=svEffects[:,3]) #c=colors
# plt.colorbar()
# plt.show()


