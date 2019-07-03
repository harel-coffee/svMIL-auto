"""	
	Input pairs of non-coding SVs and genes and coding SVs and genes.
	For each SV, compute how many genes it affects in the non-coding way, and how many in the coding way. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import re

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

svEffects = np.empty([len(totalSVs), 3], dtype="object")
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

for pair in nonCodingPairs[:,0]:
	splitPair = pair.split("_")
	gene = splitPair[0]
	pairSample = splitPair[len(splitPair)-1]
	shortPairSampleName = pairSample.split("brca")[1]
	
	if gene not in expressionData[:,0]:
		continue
	
	geneExpression = expressionData[expressionData[:,0] == gene][0]
	sampleExpressionValue = 0 #expression values of this gene in all samples
	
	geneSamples = geneSampleRef[gene] #all samples in which this gene is affected by an SV, whether coding or non-coding
	
	matchedFullSampleNames = []
	for geneSample in geneSamples:

		shortSampleName = geneSample.split("brca")[1]
		
		#match the sample name with the expression sample name
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
				
				matchedFullSampleNames.append(sample) #All samples that should be excluded from the null distribution
			
			#Also find the expression for the current pair sample
			if re.search(shortPairSampleName, sample, re.IGNORECASE) is not None:
				#Get the last 2 numbers
				splitSample = sample.split("-")
				code = int(splitSample[len(splitSample)-1])
				
				if code < 10: #above 9 are the normal samples, which we do not want to include here
					sampleInd = samples.index(sample)
					
					sampleExpressionValue = float(geneExpression[sampleInd])
				
				
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
	
	#Get the expression z-score for this pair
	if np.std(negativeSampleExpressionValues) == 0:
		continue

	z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
	print pair
	print z

exit()
#### append the expression z-score as another column


#Plot
plt.scatter(svEffects[:,2], svEffects[:,1]) #c=colors
plt.show()


