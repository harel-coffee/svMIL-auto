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

#Get the number of DEGs linked to each SV in the non-shuffled case
nonCodingDegPairs = np.loadtxt(sys.argv[3], dtype="object")
codingDegPairs = np.loadtxt(sys.argv[4], dtype="object")


totalSVs = np.union1d(nonCodingDegPairs[:,0], codingDegPairs[:,0])
print len(totalSVs)

svEffects = np.empty([len(totalSVs), 4], dtype="object")
svEffects[:,1] = 0
svEffects[:,2] = 0
ind = 0
for svInd in range(0, nonCodingDegPairs.shape[0]):
	svEffects[ind,0] = nonCodingDegPairs[svInd,0]
	svEffects[ind,1] = int(nonCodingDegPairs[svInd,1])
	print int(nonCodingDegPairs[svInd,1])
	ind += 1

for svInd in range(0, codingDegPairs.shape[0]):
	sv = codingDegPairs[svInd,0]
	if sv in svEffects[:,0]:
		rowInd = svEffects[:,0] == sv
		svEffects[rowInd,2] = int(codingDegPairs[svInd,1])
	else:
		svEffects[ind,0] = sv
		svEffects[ind,2] = int(codingDegPairs[svInd,1])
		ind += 1
	
print svEffects

#Get the properties of SVs beforehand to match these later
svProperties = []

for sv in svEffects[:,0]:
	
	splitSV = sv.split("_")
	chrom = splitSV[0]
	size = int(splitSV[5]) - int(splitSV[1])
	sample = splitSV[6]

	svProperties.append([sv, chrom + "_" + str(size) + "_" + sample])

svProperties = np.array(svProperties, dtype="object")
print svProperties

#Go through the shuffled files coding & non-coding and compute how often the real noncoding/coding is higher for that SV
# 
# import glob
# 
# #First determine non-coding potential
# shuffledNCFiles = glob.glob(sys.argv[5] + "/*degPairsNonCoding.txt")
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
# 	
# 	for svInd in range(0, nonCodingPairs.shape[0]):
# 		sv = nonCodingPairs[svInd,0]
# 		splitSV = sv.split("_")
# 		
# 		#Find out which SV this shuffled one originally is so that we can match back to the DEG pairs
# 		svSize = int(splitSV[5]) - int(splitSV[1])
# 		svMatchCode = splitSV[0] + "_" + str(svSize) + "_" + splitSV[6]
# 		trueSV = svProperties[svProperties[:,1] == svMatchCode]
# 		
# 		#Make the plot based on the number of DEGs that this SV is linked to
# 		nonCodingSVCounts[sv] = int(nonCodingPairs[svInd,1])
# 		
# 	#Get the corresponding coding file
# 	splitFileName = shuffledFile.split("_")
# 	
# 	codingPairsFile = "_".join(splitFileName[0:len(splitFileName)-1]) + "_degPairsCoding.txt" 
# 	
# 	codingPairs = np.loadtxt(codingPairsFile, dtype="object")
# 	
# 	codingSVCounts = dict()
# 	for svInd in range(0, codingPairs.shape[0]):
# 		sv = codingPairs[svInd,0]
# 		
# 		#Make the plot based on the number of DEGs that this SV is linked to
# 		codingSVCounts[sv] = int(codingPairs[svInd,1])
# 		
# 	allSVs = np.union1d(nonCodingSVCounts.keys(), codingSVCounts.keys())
# 	for sv in allSVs:
# 		
# 		ncPotential = 0
# 		if sv in nonCodingSVCounts and sv in codingSVCounts:
# 
# 			ncPotential = nonCodingSVCounts[sv]
# 		
# 		#Match to the true SV
# 		splitSV = sv.split("_")
# 		svSize = int(splitSV[5]) - int(splitSV[1])
# 		svMatchCode = splitSV[0] + "_" + str(svSize) + "_" + splitSV[6]
# 		
# 		trueSV = svProperties[svProperties[:,1] == svMatchCode]
# 		
# 		if len(trueSV) > 0: #SV match found
# 			#Is the nc/coding higher in the true case than in this random case?
# 			
# 			
# 			trueSVEffects = svEffects[svEffects[:,0] == trueSV[0][0]][0]
# 
# 			ncPotentialTrueSV = trueSVEffects[1]
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
# 			
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
# np.savetxt('Output/significantNCProportion_DEG.txt', svSignificance, fmt='%s', delimiter='\t')
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
# np.savetxt('Output/significantNCProportion_multipleTestCorrected_DEG.txt', svSignificanceCorrected, fmt='%s', delimiter='\t')
# exit()
# print "plotting pairs"
# plt.scatter(svEffects[:,2], svEffects[:,1]) #c=colors
# 
# print "plotting significance: "
# #Do overlay because the colormap is not working separately
# svSignificanceCorrected = np.loadtxt('Output/significantNCProportion_multipleTestCorrected_DEG.txt', dtype="object")
# for svInd in range(0, svEffects.shape[0]):
# 	sv = svEffects[svInd,0]
# 	
# 	if sv in svSignificanceCorrected[:,0]:
# 		plt.scatter(svEffects[svInd,2], svEffects[svInd,1], marker="*", color='red') #plot the significant SVs
# plt.show()
# exit()

svSignificanceCorrected = np.loadtxt('Output/significantNCProportion_multipleTestCorrected_DEG.txt', dtype="object")


###1. Set threshold on high/low coding effects
codingThreshold = 20

#filter by threshold
filteredSVs = []
otherSVs = []
filteredSignificantSVs = []
for sv in svEffects:

	if sv[2] <= codingThreshold and sv[1] > 0:
		filteredSVs.append(sv)
		if sv[0] in svSignificanceCorrected[:,0]:
			filteredSignificantSVs.append(sv)
	elif sv[1] > 0:
		otherSVs.append(sv)

filteredSVs = np.array(filteredSVs, dtype='object')
otherSVs = np.array(otherSVs, dtype='object')
filteredSignificantSVs = np.array(filteredSignificantSVs, dtype='object')

### Check the DEG genes of these groups
degPairs = np.load(sys.argv[1] + "_nonCodingPairDEGs.npy")

degs = dict()
allDegGenes = []
allDegGenesHighCoding = []
allDegGenesSignificant = []
allDegGenesSignificantThreshold = []
signSVEffects = []
for pair in degPairs:
	splitPair = pair[0].split("_")
	sv = "_".join(splitPair[1:])
	
	if sv in filteredSVs[:,0]:
		
		if splitPair[0] not in allDegGenes:
			allDegGenes.append(splitPair[0])
	
	if sv in otherSVs[:,0]:
		
		if splitPair[0] not in allDegGenesHighCoding:
			allDegGenesHighCoding.append(splitPair[0])

	if sv in svSignificanceCorrected[:,0]:
		
		if splitPair[0] not in allDegGenesSignificant:
			allDegGenesSignificant.append(splitPair[0])
			
		if sv not in degs:
			degs[sv] = []
		degs[sv].append(splitPair[0])
		
		
		
	
	if sv in filteredSignificantSVs[:,0]:
		
		if splitPair[0] not in allDegGenesSignificantThreshold:
			allDegGenesSignificantThreshold.append(splitPair[0])


print len(allDegGenes)
print len(allDegGenesHighCoding)
print len(allDegGenesSignificant)
print len(allDegGenesSignificantThreshold)

np.savetxt('Output/allDegGenesLowCoding.txt', allDegGenes, delimiter='\t', fmt='%s')
np.savetxt('Output/allDegGenesHighCoding.txt', allDegGenesHighCoding, delimiter='\t', fmt='%s')
np.savetxt('Output/allDegGenesSignificant.txt', allDegGenesSignificant, delimiter='\t', fmt='%s')
np.savetxt('Output/allDegGenesSignificantThreshold.txt', allDegGenesSignificantThreshold, delimiter='\t', fmt='%s')

#Check for every DEG gene linked to the significant SVs if there are also other samples in which this gene is affected (in coding way)
codingDegPairs = np.load(sys.argv[1] + "_codingPairDEGs.npy")
svsWithOtherSampleEvidence = dict()
signGeneSVPairs = []
for pair in degPairs:
	splitPair = pair[0].split("_")
	sv = "_".join(splitPair[1:])
	gene = splitPair[0]

	if sv in filteredSVs[:,0]:
		signGeneSVPairs.append(gene + "_" + sv)
		#check if the gene linked to the SV is found more often
		sample = splitPair[len(splitPair)-1]
		for pair2 in codingDegPairs:
			splitPair = pair2[0].split("_")
			sv2 = "_".join(splitPair[1:])
			sample2 = splitPair[len(splitPair)-1]

			#we are looking for effects in other samples, so the SV should not be the same
			if sample != sample2: #if the samples are not the same, the SV is also never the same
				#check if the gene is the same
				gene2 = splitPair[0]
				if gene == gene2:
					if sv not in svsWithOtherSampleEvidence:
						svsWithOtherSampleEvidence[sv] = []
					svsWithOtherSampleEvidence[sv].append(gene2)

print len(svsWithOtherSampleEvidence)
multiSampleCount = 0
genes = []
for sv in svsWithOtherSampleEvidence:
	if len(svsWithOtherSampleEvidence[sv]) != len(np.unique(svsWithOtherSampleEvidence[sv])):
		print sv
		print np.unique(svsWithOtherSampleEvidence[sv])
		print "Number of genes also in coding: ", len(svsWithOtherSampleEvidence[sv])
		print "Number of unique genes: ", len(np.unique(svsWithOtherSampleEvidence[sv]))
		multiSampleCount += 1
		genes += list(np.unique(svsWithOtherSampleEvidence[sv]))
print multiSampleCount
np.savetxt('Output/leftSideGenes_multiCoding.txt', np.unique(genes), fmt='%s')
exit()
#Which of the significant pairs are also found in the naive way? any that is not found in that way?
naiveDegPairs = np.loadtxt('naiveTadDisr_nonCodingDEGs.txt', dtype="object")

foundCount = 0
notFoundCount = 0
for degPair in signGeneSVPairs:
	if degPair in naiveDegPairs[:,0]:
		print "pair found: ", degPair
		foundCount += 1
	else:
		print "pair not found: ", degPair
		notFoundCount += 1
		# for sv in somaticSVs:
		# 	svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		# 	splitPair = degPair.split("_")
		# 	if svStr == "_".join(splitPair[1:]):
		# 		print sv[8].svType
		
print foundCount
print notFoundCount
exit()

## Make a file with the scores of the filtered SVs vs the rest
#1. Get the scores of the filteredSVs
lowCodingFeatures = []
highCodingFeatures = []
for pair in nonCodingPairs:
	splitPair = pair[0].split("_")
	svEntries = splitPair[1:]
	sv = "_".join(svEntries)
	if sv in filteredSVs:
		lowCodingFeatures.append(pair)
	if sv in otherSVs:
		highCodingFeatures.append(pair)

lowCodingFeatures = np.array(lowCodingFeatures, dtype="object")
highCodingFeatures = np.array(highCodingFeatures, dtype='object')

leftFeatures = lowCodingFeatures[:,1:]
rightFeatures = highCodingFeatures[:,1:]

#Check for the significant SVs vs the rest with low coding
# svSignificanceCorrected = np.loadtxt('Output/significantNCProportion_multipleTestCorrected_DEG.txt', dtype="object")
# significantFeatures = []
# for pair in nonCodingPairs:
# 	splitPair = pair[0].split("_")
# 	svEntries = splitPair[1:]
# 	sv = "_".join(svEntries)
# 	if sv in filteredSVs and sv not in svSignificanceCorrected[:,0]:
# 		lowCodingFeatures.append(pair)
# 	if sv in svSignificanceCorrected[:,0]:
# 		significantFeatures.append(pair)
# 		
# lowCodingFeatures = np.array(lowCodingFeatures, dtype="object")
# significantFeatures = np.array(significantFeatures, dtype='object')
# 
# leftFeatures = lowCodingFeatures[:,1:]
# rightFeatures = significantFeatures[:,1:]


#Make a plot showing how frequently each feature is gained or lost in these sets
eQTLLosses = leftFeatures[:,0].astype(float)
enhancerLosses = leftFeatures[:,1].astype(float)
promoterLosses = leftFeatures[:,2].astype(float)
cpgLosses = leftFeatures[:,3].astype(float)
tfLosses = leftFeatures[:,4].astype(float)
hicLosses = leftFeatures[:,5].astype(float)
h3k9me3Losses = leftFeatures[:,6].astype(float)
h3k4me3Losses = leftFeatures[:,7].astype(float)
h3k27acLosses = leftFeatures[:,8].astype(float)
h3k27me3Losses = leftFeatures[:,9].astype(float)
h3k4me1Losses = leftFeatures[:,10].astype(float)
h3k36me3Losses = leftFeatures[:,11].astype(float)
dnaseLosses = leftFeatures[:,12].astype(float)	
	

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses)]
lossData = np.array(lossData)
#lossData = lossData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
#lossData = -np.log(lossData)
print lossData

width = 0.35

plt.bar(np.arange(len(lossData)), lossData, width, label='Low coding', color='blue')
# plt.xticks(range(0, len(lossData)),
	   # ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()

eQTLLosses = rightFeatures[:,0].astype(float)
enhancerLosses = rightFeatures[:,1].astype(float)
promoterLosses = rightFeatures[:,2].astype(float)
cpgLosses = rightFeatures[:,3].astype(float)
tfLosses = rightFeatures[:,4].astype(float)
hicLosses = rightFeatures[:,5].astype(float)
h3k9me3Losses = rightFeatures[:,6].astype(float)
h3k4me3Losses = rightFeatures[:,7].astype(float)
h3k27acLosses = rightFeatures[:,8].astype(float)
h3k27me3Losses = rightFeatures[:,9].astype(float)
h3k4me1Losses = rightFeatures[:,10].astype(float)
h3k36me3Losses = rightFeatures[:,11].astype(float)
dnaseLosses = rightFeatures[:,12].astype(float)

lossData = [np.sum(eQTLLosses), np.sum(enhancerLosses), np.sum(promoterLosses), np.sum(cpgLosses),
			np.sum(tfLosses), np.sum(hicLosses), np.sum(h3k9me3Losses), np.sum(h3k4me3Losses), np.sum(h3k27acLosses),
			np.sum(h3k27me3Losses), np.sum(h3k4me1Losses), np.sum(h3k36me3Losses), np.sum(dnaseLosses)]
lossData = np.array(lossData)
#lossData = lossData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
#lossData = -np.log(lossData)
print lossData

plt.bar(np.arange(len(lossData)) + width, lossData, width, label='High coding', color='red')
plt.xticks(np.arange(len(lossData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
#plt.ylim([0,10])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/leftRight_losses.svg')

#Gains

eQTLGains = leftFeatures[:,13].astype(float)
enhancerGains = leftFeatures[:,14].astype(float)
promoterGains = leftFeatures[:,15].astype(float)
cpgGains = leftFeatures[:,16].astype(float)
tfGains = leftFeatures[:,17].astype(float)
hicGains = leftFeatures[:,18].astype(float)
h3k9me3Gains = leftFeatures[:,19].astype(float)
h3k4me3Gains = leftFeatures[:,20].astype(float)
h3k27acGains = leftFeatures[:,21].astype(float)
h3k27me3Gains = leftFeatures[:,22].astype(float)
h3k4me1Gains = leftFeatures[:,23].astype(float)
h3k36me3Gains = leftFeatures[:,24].astype(float)
dnaseGains = leftFeatures[:,25].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains)]

gainData = np.array(gainData)
#gainData = gainData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
#gainData = -np.log(gainData)
print gainData

plt.bar(np.arange(len(gainData)), gainData, width, label='Low coding', color='blue')
# plt.xticks(np.arange(len(gainData)),
# 		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
# plt.show()
	
eQTLGains = rightFeatures[:,13].astype(float)
enhancerGains = rightFeatures[:,14].astype(float)
promoterGains = rightFeatures[:,15].astype(float)
cpgGains = rightFeatures[:,16].astype(float)
tfGains = rightFeatures[:,17].astype(float)
hicGains = rightFeatures[:,18].astype(float)
h3k9me3Gains = rightFeatures[:,19].astype(float)
h3k4me3Gains = rightFeatures[:,20].astype(float)
h3k27acGains = rightFeatures[:,21].astype(float)
h3k27me3Gains = rightFeatures[:,22].astype(float)
h3k4me1Gains = rightFeatures[:,23].astype(float)
h3k36me3Gains = rightFeatures[:,24].astype(float)
dnaseGains = rightFeatures[:,25].astype(float)

gainData = [np.sum(eQTLGains), np.sum(enhancerGains), np.sum(promoterGains), np.sum(cpgGains),
			np.sum(tfGains), np.sum(hicGains), np.sum(h3k9me3Gains), np.sum(h3k4me3Gains), np.sum(h3k27acGains),
			np.sum(h3k27me3Gains), np.sum(h3k4me1Gains), np.sum(h3k36me3Gains), np.sum(dnaseGains)]

gainData = np.array(gainData)
#gainData = gainData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
#gainData = -np.log(gainData)

print gainData

plt.bar(np.arange(len(gainData)) + width, gainData, width, label='High coding',color='red')
plt.xticks(np.arange(len(gainData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)

#plt.ylim([0,10])
plt.legend(loc='best')
plt.tight_layout()
plt.show()

exit()


svSignificanceCorrected = np.loadtxt('Output/significantNCProportion_multipleTestCorrected_DEG.txt', dtype="object")
print "Number of significant gene SV pairs: ", svSignificanceCorrected.shape

#Number of unique samples

samples = []
svWithGenes = 0
geneCounts = []
for sv in svSignificanceCorrected:
	splitSV = sv[0].split("_")
	if splitSV[len(splitSV)-1] not in samples:
		samples.append(splitSV[len(splitSV)-1])

print "no of samples: ", len(samples), samples

#Check which SVs have nc potential but no coding SVs
for sv in svSignificanceCorrected:
	
	svEffect = svEffects[svEffects[:,0] == sv[0]]
	print svEffect


#Find the DEG genes that these SVs are linked to
print nonCodingPairs


exit()

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
# # 
# # exit()
# # 		
# # 		#
# # 	
# # #plt.show()
# # 
# # 
# # exit()


