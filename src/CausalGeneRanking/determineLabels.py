from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
# 
# # Combining DEGs, recurrence and germline
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# germlineScores = np.loadtxt(sys.argv[2], dtype="object")
# degLabels = np.loadtxt(sys.argv[3], dtype="object")
# somaticRanks = np.loadtxt(sys.argv[4], dtype="object")
# germlineRanks = np.loadtxt(sys.argv[5], dtype="object")
# 
# #get the deg genes
# degs = []
# for pair in degLabels[:,0]:
# 	splitPair = pair.split("_")
# 	gene = splitPair[0]
# 	degs.append(gene)
# 
# #1. How many genes are linked to germline SVs? And how many of these are also linked to somatic? Which are exclusive? 
# 
# somaticGenes = []
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 
# 	if score == True:
# 		
# 		if splitScore[0] not in somaticGenes:
# 			somaticGenes.append(splitScore[0])
# 
# 
# germlineGenes = []
# for geneSVPair in germlineScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 
# 	if score == True:
# 		
# 		if splitScore[0] not in germlineGenes:
# 			germlineGenes.append(splitScore[0])
# 
# print len(somaticGenes)
# print len(germlineGenes)
# 
# print len(np.intersect1d(somaticGenes, germlineGenes))
# somaticUniqueGenes = np.setdiff1d(somaticGenes, germlineGenes)
# germlineUniqueGenes = np.setdiff1d(germlineGenes, somaticGenes)
# 
# #Of the somatic unique genes, which are recurrent and which are DEG?
# recurrenceThreshold = 3
# somaticRecurrentGenes = []
# for gene in somaticRanks:
# 	
# 	samples = gene[31]
# 	
# 	if samples == "None":
# 		continue
# 	
# 	splitSamples = samples.split(",")
# 	sampleNum = len(splitSamples)
# 	
# 	if sampleNum >= recurrenceThreshold and gene[0] in somaticUniqueGenes:
# 		
# 		somaticRecurrentGenes.append(gene[0])
# 
# print "Number of somatic recurrent genes: ", len(somaticRecurrentGenes)
# somaticUniqueDegs = np.intersect1d(somaticUniqueGenes, degs)
# print "Number of somatic unique genes that are DEG: ", len(somaticUniqueDegs)
# 
# #Of the germline unique genes, which are recurrent and which are DEG? 
# recurrenceThreshold = 5
# germlineRecurrentGenes = []
# for gene in germlineRanks:
# 	
# 	samples = gene[31]
# 	
# 	if samples == "None":
# 		continue
# 	
# 	splitSamples = samples.split(",")
# 	sampleNum = len(splitSamples)
# 	
# 	if sampleNum >= recurrenceThreshold and gene[0] in germlineUniqueGenes:
# 		
# 		germlineRecurrentGenes.append(gene[0])
# 
# print "Number of germline genes that are recurrent: ", len(germlineRecurrentGenes)
# print "Number of germline unique genes that are DEG: ", len(np.intersect1d(germlineUniqueGenes, degs))
# 
# #How many SVs do we get when we only use the ones that are linked to the unique germline set?
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in germlineScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# negativeBags = [] #based on genes unique to germline
# recurrentNegativeBags = [] #based on genes unique to germline that are also recurrent
# negativeFeatures = []
# genes = []
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent
# 	negativeStatus = True
# 	negativeRecurrentStatus = True
# 	for gene in genesPerBag[sv]:	
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		#Using all genes unique to germline
# 		if gene not in germlineUniqueGenes:
# 			negativeStatus = False
# 		#Using all genes that are germline + recurrent (VERY likely negative)
# 		if gene not in germlineRecurrentGenes:
# 			negativeRecurrentStatus = False
# 			
# 	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent	
# 	if negativeStatus == True:
# 		negativeBags.append(svBagContents[sv])
# 	if negativeRecurrentStatus == True:
# 		negativeFeatures.append(svBagContents[sv][0])
# 		recurrentNegativeBags.append(svBagContents[sv])
# 		genes.append(gene)
# 	
# print "Number of SVs in negative set using only genes unique to germline: ", len(negativeBags)
# print "Number of SVs in negative set using only genes unique to germline AND recurrent: ", len(recurrentNegativeBags)
# 
# #How many SVs do we get in the positive set?
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# positiveBags = [] #based on genes unique to somatic
# recurrentPositiveBags = [] #based on genes unique to somatic that are also recurrent
# degPositiveBags = [] #based on genes unique to somatic that are also DEGs
# recurrentDegPositiveBags = []
# #features = []
# genes = []
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	#An SV is only in the positive set if it does not affect ANY gene that is not unique to somatic/recurrent
# 	positiveStatus = True
# 	positiveRecurrentStatus = True
# 	positiveDegStatus = True
# 	positiveRecurrentDegStatus = True
# 	for gene in genesPerBag[sv]:	
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		#Using all genes unique to somatic
# 		if gene not in somaticUniqueGenes:
# 			positiveStatus = False
# 		#Using all genes that are somatic + recurrent
# 		if gene not in somaticRecurrentGenes:
# 			positiveRecurrentStatus = False
# 		if gene not in somaticUniqueDegs:
# 			positiveDegStatus = False
# 		if gene not in somaticUniqueDegs or gene not in somaticRecurrentGenes:
# 			positiveRecurrentDegStatus = False
# 		
# 			
# 	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent	
# 	if positiveStatus == True:
# 		positiveBags.append(svBagContents[sv])
# 		genes.append(gene)
# 	if positiveRecurrentStatus == True:
# 		recurrentPositiveBags.append(svBagContents[sv])
# 	if positiveDegStatus == True:
# 		degPositiveBags.append(svBagContents[sv])
# 		
# 	if positiveRecurrentDegStatus == True:
# 		
# 		recurrentDegPositiveBags.append(svBagContents[sv])
# 
# print "Number of SVs in positive set using only genes unique to somatic: ", len(positiveBags)
# print "Number of SVs in positive set using only genes unique to somatic AND recurrent: ", len(recurrentPositiveBags)
# print "Number of SVs in positive set using only genes unique to somatic AND DEGs: ", len(degPositiveBags)
# print "Number of SVs in positive set using only genes unique to somatic AND DEGs AND recurrent: ", len(recurrentDegPositiveBags)
# exit()
# #Load the pairs at the left & right side of the PCA
# leftPairs = np.load('Output/similarityMatrices/uniqueSomaticGermlineRecurrentLeftPairs.txt.npy')
# rightPairs = np.load('Output/similarityMatrices/uniqueSomaticGermlineRecurrentRightPairs.txt.npy')
# 
# leftGenes = []
# for pair in leftPairs:
# 	splitPair = pair.split("_")
# 	leftGenes.append(splitPair[0])
# leftGenes = np.unique(leftGenes)
# np.savetxt('Output/leftGenes.txt', leftGenes, fmt="%s")
# 
# rightGenes = []
# for pair in rightPairs:
# 	splitPair = pair.split("_")
# 	rightGenes.append(splitPair[0])
# rightGenes = np.unique(rightGenes)
# np.savetxt('Output/rightGenes.txt', rightGenes, fmt="%s")
# exit()
# 
# #Gather the features of these pairs specifically
# rightFeatures = []
# leftFeatures = []
# for pair in leftPairs:
# 	if pair in somaticScores[:,0]:
# 		pairScores = somaticScores[somaticScores[:,0] == pair,1:somaticScores.shape[1]-1]
# 		leftFeatures.append(pairScores[0])
# 	if pair in germlineScores[:,0]:
# 		pairScores = germlineScores[germlineScores[:,0] == pair,1:germlineScores.shape[1]-1]
# 		leftFeatures.append(pairScores[0])
# 
# for pair in rightPairs:
# 	if pair in somaticScores[:,0]:
# 		pairScores = somaticScores[somaticScores[:,0] == pair,1:somaticScores.shape[1]-1]
# 		rightFeatures.append(pairScores[0])
# 	if pair in germlineScores[:,0]:
# 		pairScores = germlineScores[germlineScores[:,0] == pair,1:germlineScores.shape[1]-1]
# 		rightFeatures.append(pairScores[0])


# leftFeatures = np.array(leftFeatures)
# rightFeatures = np.array(rightFeatures)

leftFeatures = np.loadtxt(sys.argv[1], dtype="object")
rightFeatures = np.loadtxt(sys.argv[2], dtype="object")

leftFeatures = leftFeatures[:,1:]
rightFeatures = rightFeatures[:,1:]


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
lossData = lossData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
lossData = -np.log(lossData)
print(lossData)

width = 0.35

plt.bar(np.arange(len(lossData)), lossData, width, label='Left', color='blue')
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
lossData = lossData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
lossData = -np.log(lossData)
print(lossData)

plt.bar(np.arange(len(lossData)) + width, lossData, width, label='Right pairs', color='yellow')
plt.xticks(np.arange(len(lossData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)
plt.ylim([0,10])
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
gainData = gainData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
gainData = -np.log(gainData)
print(gainData)

plt.bar(np.arange(len(gainData)), gainData, width, label='Left', color='blue')
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
gainData = gainData / float(leftFeatures.shape[0] + rightFeatures.shape[0])
gainData = -np.log(gainData)

print(gainData)

plt.bar(np.arange(len(gainData)) + width, gainData, width, label='Right pairs',color='yellow')
plt.xticks(np.arange(len(gainData) + width / 2),
		   ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'], rotation=90)

plt.ylim([0,10])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/leftRight_gains.svg')

exit()


# Show the features for the specific groups
# features = np.array(negativeFeatures)
# print features
# 
# #Make feature distribution plots
# #Do this here because the order of the features is different in the ranking compared to the scores file
# eQTLLosses = features[:,0].astype(float)
# enhancerLosses = features[:,1].astype(float)
# promoterLosses = features[:,2].astype(float)
# cpgLosses = features[:,3].astype(float)
# tfLosses = features[:,4].astype(float)
# hicLosses = features[:,5].astype(float)
# h3k9me3Losses = features[:,6].astype(float)
# h3k4me3Losses = features[:,7].astype(float)
# h3k27acLosses = features[:,8].astype(float)
# h3k27me3Losses = features[:,9].astype(float)
# h3k4me1Losses = features[:,10].astype(float)
# h3k36me3Losses = features[:,11].astype(float)
# dnaseLosses = features[:,12].astype(float)
# 
# lossData = [eQTLLosses, enhancerLosses, promoterLosses, cpgLosses, tfLosses, hicLosses, h3k9me3Losses, h3k4me3Losses, h3k27acLosses, h3k27me3Losses, h3k4me1Losses, h3k36me3Losses, dnaseLosses]
# 
# fig, ax = plt.subplots()
# bp = plt.boxplot(lossData)
# ax.set_xticklabels(['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'])
# ax.set_ylim(0,18)
# plt.xticks(rotation=70)
# plt.show()
# 
# 
# eQTLGains = features[:,13].astype(float)
# enhancerGains = features[:,14].astype(float)
# promoterGains = features[:,15].astype(float)
# cpgGains = features[:,16].astype(float)
# tfGains = features[:,17].astype(float)
# hicGains = features[:,18].astype(float)
# h3k9me3Gains = features[:,19].astype(float)
# h3k4me3Gains = features[:,20].astype(float)
# h3k27acGains = features[:,21].astype(float)
# h3k27me3Gains = features[:,22].astype(float)
# h3k4me1Gains = features[:,23].astype(float)
# h3k36me3Gains = features[:,24].astype(float)
# dnaseGains = features[:,25].astype(float)
# 
# 
# gainData = [eQTLGains, enhancerGains, promoterGains, cpgGains, tfGains, hicGains, h3k9me3Gains, h3k4me3Gains, h3k27acGains, h3k27me3Gains, h3k4me1Gains, h3k36me3Gains, dnaseGains]
# 
# fig, ax = plt.subplots()
# bp = plt.boxplot(gainData)
# ax.set_xticklabels(['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'DNAseI'])
# plt.xticks(rotation=70)
# ax.set_ylim(0,18)
# plt.show()
# 
# 
