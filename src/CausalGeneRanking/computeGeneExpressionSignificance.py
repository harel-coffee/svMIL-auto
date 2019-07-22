"""
	Compute the differential expression difference between genes with SVs and genes without SVs. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import stats

#Plot the number of samples that every gene has in a histogram
geneScoreFile = sys.argv[1]
geneScores = np.loadtxt(geneScoreFile, dtype="object")

sampleCounts = dict()

for gene in geneScores:
	
	samples = gene[31].split(",")
	if samples[0] != "None":
		print samples[0]
		if len(samples) not in sampleCounts:
			sampleCounts[len(samples)] = 0
		sampleCounts[len(samples)] += 1

plt.bar(sampleCounts.keys(), sampleCounts.values())
# plt.show()

#Make the gene subsets given a threshold of number of samples
threshold = 2
filteredGenes = []
for gene in geneScores:
	samples = gene[31].split(",")
	if len(samples) > threshold:
		
		filteredGenes.append(gene)
		
filteredGenes = np.array(filteredGenes, dtype="object")
print filteredGenes.shape

#Write to outfile to check cosmic overlap etc
header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\ttotal\tsamples"
				
#Write to numpy output file	
np.savetxt(geneScoreFile + "_filtered.txt", filteredGenes, delimiter='\t', fmt='%s', header=header)

#Get the expression values for the samples in the positive subset
expressionFile = sys.argv[2]

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

pValues = []
for gene in filteredGenes:
	if gene[0] not in expressionData[:,0]:
		continue
	geneExpression = expressionData[expressionData[:,0] == gene[0]][0]
	sampleExpressionValues = [] #expression values of this gene in all samples
	
	geneSamples = gene[31].split(",")
	matchedFullSampleNames = []
	for geneSample in geneSamples:

		shortSampleName = geneSample.split("brca")[1]
		
		#match the sample name with the expression sample name
		for sampleInd in range(0, len(samples)):
			sample = samples[sampleInd]
			if re.search(shortSampleName, sample, re.IGNORECASE) is not None:
				matchedFullSampleNames.append(sample) #keep this to check later for the negative set
				#Get the last 2 numbers
				splitSample = sample.split("-")
				code = int("".join(list(splitSampleName[3])[0:2]))
				
				if code < 10: #above 9 are the normal samples, which we do not want to include here
					sampleInd = samples.index(sample)
					
					sampleExpressionValues.append(float(geneExpression[sampleInd]))
			
				
	#Get 5 random samples that are not affecting this gene
	unmatchedSamples = np.setdiff1d(samples[1:len(samples)-1], matchedFullSampleNames) #exclude hybrid ref
	negativeSamples = []
	for sample in unmatchedSamples: #sample tumor samples, exclude normals
		splitSample = sample.split("-")
		code = int("".join(list(splitSampleName[3])[0:2]))
		
		if code < 10: 
			negativeSamples.append(sample)
		# if len(negativeSamples) == len(matchedFullSampleNames):
		# 	break
	
	#Get the expression of these samples
	negativeSampleExpressionValues = []
	for sample in negativeSamples:
		sampleInd = samples.index(sample)				
		negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
	
	
	#Do a t-test and compute the p-value for this gene
	z = (np.mean(sampleExpressionValues) - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
			
	pValue = stats.norm.sf(abs(z))*2
	
	# posMean = np.mean(sampleExpressionValues)
	# posStd = np.std(sampleExpressionValues)
	# negMean = np.mean(negativeSampleExpressionValues)
	# negStd = np.std(negativeSampleExpressionValues)
	# 
	# pValue = stats.ttest_ind_from_stats(posMean, posStd, len(sampleExpressionValues), negMean, negStd, len(negativeSampleExpressionValues))[1]
	pValues.append([gene[0], pValue])
	
	


pValues = np.array(pValues, dtype="object")

pValues = pValues[pValues[:,1].argsort()]
print pValues
from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(pValues[:,1], method='bonferroni')

filteredPValues = pValues[reject]

signGenes = []
signCount = 0
for pValue in filteredPValues:
	
	gene = filteredGenes[filteredGenes[:,0] == pValue[0]][0]
	signGenes.append(gene)
	
	print pValue
	signCount += 1
print "Number of significant genes: ", signCount

signGenes = np.array(signGenes)
#Write to outfile to check cosmic overlap etc
header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\ttotal\tsamples"
				
#Write to numpy output file	
np.savetxt(geneScoreFile + "_signgt3.txt", signGenes, delimiter='\t', fmt='%s', header=header)


exit()		


#Get the expression values for the sample sin the negative subset