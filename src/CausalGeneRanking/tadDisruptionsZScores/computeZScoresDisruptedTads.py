"""
	
	First make a plot of just the expression in TADs for which we have a non-coding SV disrupting the boundary

"""
import sys

path = sys.argv[4]
sys.path.insert(1, path)

#this code depends on the input parser from the linking part. This is quick and dirty, is there a better solution?
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser
import numpy as np
import os

import glob
import re
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from genomicShuffler import GenomicShuffler

###parameters
geneNameConversionFile = sys.argv[1]
expressionFile = sys.argv[2]
mutationDir = sys.argv[3]
outDir = sys.argv[5]

specificOutDir = outDir + '/tadDisruptionsZScores/'

if not os.path.exists(specificOutDir):
	os.makedirs(specificOutDir)

#For each TAD, determine which genes are there

#first get all genes and their positions
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
#allGenes = causalGenes

genes = []
for gene in allGenes:
	genes.append([gene[0], gene[1], gene[2], gene[3].name])

genes = np.array(genes, dtype='object')


#also use a map for the gene names
geneNameConversionMap = dict()
with open(geneNameConversionFile, 'r') as inF:
	
	lineCount = 0
	for line in inF:
		
		if lineCount < 1:
			lineCount += 1
			continue
		line = line.strip()
		splitLine = line.split("\t")
		ensgId = splitLine[3]
		splitEnsgId = ensgId.split('.') #we only keep everything before the dot
		geneName = splitLine[4]
		geneNameConversionMap[splitEnsgId[0]] = geneName

#For each patient and gene, find the TAD that the gene is in. If the TAD is in the disrupted patients list, add the gene expression to the disrupted set.
#Otherwise, add the gene expression to the non-disrutped set.

#get the gene expression

expressionData = []
samples = []
with open(expressionFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		if lineCount == 0:
			samples = ['']
			samples += line.split("\t")

			lineCount += 1
			continue
		splitLine = line.split("\t")
		fullGeneName = splitLine[0]
		if settings.general['source'] == 'HMF':
			if fullGeneName not in geneNameConversionMap:
				continue
			geneName = geneNameConversionMap[fullGeneName] #get the gene name rather than the ENSG ID
		else:
			geneName = fullGeneName.split("|")[0]
			if geneName == 'gene_id':
				continue #skip this line

		data = splitLine[1:len(splitLine)]

		fixedData = [geneName]
		fixedData += data

		expressionData.append(fixedData)

expressionData = np.array(expressionData, dtype="object")
print(expressionData)

#Get all SVs
if settings.general['source'] == 'HMF':
	svDir = settings.files['svDir']
	svData = InputParser().getSVsFromFile_hmf(svDir)
elif settings.general['source'] == 'TCGA':
	svData = InputParser().getSVsFromFile(settings.files['svFile'], '')

else:
	print('Other data sources not supported')
	exit(1)
#fix this
filteredSVs = svData

#Get the TADs.
#For each TAD, if there is an SV disrupting either of the boundaries, it goes into the disrupted class for that patient.
#All other patients that do not have SVs disrupting this TAD go into the negative group.
tadFile = settings.files['tadFile']

print("Getting TADs")
tadData = InputParser().getTADsFromFile(tadFile)

tadDisruptions = dict() #keep each TAD, and a list of which patients have an SV disrupting that TAD.
for tad in tadData:

	tadStr = tad[0] + '_' + str(tad[1]) + '_' + str(tad[2])
	tadDisruptions[tadStr] = []

	#Check if there is any SV overlapping the TAD on either side.

	#for intra-chromosomal SVs, we check if these overlap the boundary of the TAD.
	#for inter-chromosomal SVs, we check if these are inside a TAD.

	#get a subset of SVs on the right chromosome
	svChr1Subset = filteredSVs[filteredSVs[:,0] == tad[0]]
	intraSVSubset = svChr1Subset[svChr1Subset[:,0] == svChr1Subset[:,3]]

	#In the intra set, check if there are SVs that start before the TAD, but end within the TAD.

	#This TAD is disrupted if an SV starts or ends in it.
	#so an SV either starts before the end, and ends after the start
	#or the sv ends after the start, but starts before the end. 

	#SVs that start before the TAD end, and end after the TAD start. 
	startMatches = (intraSVSubset[:,1] >= tad[1]) * (intraSVSubset[:,1] <= tad[2]) * (intraSVSubset[:,5] >= tad[2])
	endMatches = (intraSVSubset[:,5] >= tad[1]) * (intraSVSubset[:,1] <= tad[1]) * (intraSVSubset[:,5] <= tad[2])
	
	#startMatches = intraSVSubset[:,1] <= tad[2]
	#endMatches = intraSVSubset[:,5] >= tad[1]
			
	#either of these must be true for the TAD to be disrupted by an intra-chromosomal SV.
	allMatches = intraSVSubset[startMatches + endMatches]
	
	for match in allMatches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]

		if match[7] not in tadDisruptions[tadStr]:
			
			tadDisruptions[tadStr].append([match[7], match])

	#then repeat for interchromosomal SVs.
	#check if there is any bp in this TAD to count it as disrupted.
	
	#get a subset of all SVs that are interchromosomal, then subset for the ones that are either on chr1 or chr2
	
	interSVSubset = filteredSVs[filteredSVs[:,0] != filteredSVs[:,3]]
	interSVsChr1 = interSVSubset[interSVSubset[:,0] == tad[0]]
	interSVsChr2 = interSVSubset[interSVSubset[:,3] == tad[0]]
	
	#which bps of the chr1 set are in the TAD
	chr1Matches = (interSVsChr1[:,1] >= tad[1]) * (interSVsChr1[:,1] <= tad[2])
	chr2Matches = (interSVsChr2[:,5] >= tad[1]) * (interSVsChr2[:,5] <= tad[2])
	
	allChr1Matches = interSVsChr1[chr1Matches]
	allChr2Matches = interSVsChr2[chr2Matches]
	
	
	for match in allChr1Matches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		
		if match[7] not in tadDisruptions[tadStr]:

			tadDisruptions[tadStr].append([match[7], match])
			
	for match in allChr2Matches:
		
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		
		if match[7] not in tadDisruptions[tadStr]:

			tadDisruptions[tadStr].append([match[7], match])
	
#check which TADs we have with disruptions, validate

disrCount = 0
nonDisrCount = 0
for tad in tadData:
	
	tadStr = tad[0] + '_' + str(tad[1]) + '_' + str(tad[2])
	
	if len(tadDisruptions[tadStr]) > 0:
		disrCount += 1
	else:
		nonDisrCount += 1
		
print('disrupted tads: ', disrCount)
print('non-disrupted tads: ', nonDisrCount)

allPatientsWithDisruptions = dict()
for tad in tadDisruptions:
	for match in tadDisruptions[tad]:
		allPatientsWithDisruptions[match[0]] = 0

#check how many unique SVs disrupt a TAD pair
disruptingSVs = dict()
typeDistribution = dict()
for tad in tadDisruptions:
	
	for matchList in tadDisruptions[tad]:
		
		match = matchList[1]
		svStr = match[0] + '_' + str(match[1]) + '_' + str(match[2]) + '_' + match[3] + '_' + str(match[4]) + '_' + str(match[5]) + '_' + match[7]
		svType = match[8].svType
		
		
		if svType not in typeDistribution:
			typeDistribution[svType] = 0
		
		if svStr not in disruptingSVs: 
			typeDistribution[svType] += 1

			disruptingSVs[svStr] = 0


print('disrupting SVs: ', len(disruptingSVs))
print('per type: ', typeDistribution)
#
# exit()
#
# #to shuffle across patients, first transpose, the shuffle, then transpose back.
#
# genes = expressionData[:,0]
# expression = expressionData[:,1:]
# expressionT = expression.T
# print(expressionT)
# print(expressionT.shape)
# np.random.shuffle(expressionT)
# print(expressionT)
# print(expressionT.shape)
# shuffledExpression = expressionT.T
# print(shuffledExpression)
# print(shuffledExpression.shape)
#
# shuffledExpressionData = np.empty(expressionData.shape, dtype='object')
# shuffledExpressionData[:,0] = genes
# shuffledExpressionData[:,1:] = shuffledExpression
#
# print(samples)
# print(shuffledExpressionData)
#
# expressionData = shuffledExpressionData

mutDir = outDir + '/patientGeneMutationPairs/'
snvPatients = np.load(mutDir + 'snvPatients.npy', allow_pickle=True, encoding='latin1').item()

svPatientsDel = np.load(mutDir + 'svPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svPatientsInv = np.load(mutDir + 'svPatientsInv.npy', allow_pickle=True, encoding='latin1').item()
svPatientsItx = np.load(mutDir + 'svPatientsItx.npy', allow_pickle=True, encoding='latin1').item()

cnvPatientsDel = np.load(mutDir + 'cnvPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()


disruptedTadExpression = []
nonDisruptedTadExpression = []

disruptedPairs = dict() #store for each patient a dict with genes, and in there, the expression of the gene. 
nonDisruptedPairs = dict()


#in case of tcga data, update the sample names that do not match the SVs.
if settings.general['source'] == 'TCGA':
	fixedSamples = []
	for sample in samples:

		if sample == '':
			fixedSamples.append('')
			continue
		if sample == 'Hybridization REF':
			continue

		splitPatientID = sample.split("-")
		shortPatientID = settings.general['cancerType'] + splitPatientID[2]
		
		if shortPatientID not in allPatientsWithDisruptions:
			continue #some patients never have svs, so no need to look at those. 
		
		fixedSamples.append(shortPatientID)
	samples = fixedSamples

#tcga data is a mess and many patients are missing from the mutations. So add them back here..
if settings.general['source'] == 'TCGA':
	for sample in samples:
		if sample not in snvPatients:
			snvPatients[sample] = []
		if sample not in svPatientsDel:
			svPatientsDel[sample] = []
		if sample not in svPatientsDup:
			svPatientsDup[sample] = []
		if sample not  in svPatientsInv:
			svPatientsInv[sample] = []
		if sample not in svPatientsItx:
			svPatientsItx[sample] = []
		if sample not in cnvPatientsAmp:
			cnvPatientsAmp[sample] = []
		if sample not in cnvPatientsDel:
			cnvPatientsDel[sample] = []

for patient in samples:
	
	if patient == '' or patient == 'Hyridization REF':
		continue
	
	if patient not in disruptedPairs:
		disruptedPairs[patient] = dict()

for gene in allGenes:
	if gene[3].name not in nonDisruptedPairs:
		nonDisruptedPairs[gene[3].name] = dict()

#for each gene, label the genes that overlap multiple TADs.
filteredGenes = []
for gene in allGenes:
	
	#check which tads is overlaps with
	tadChrSubset = tadData[tadData[:,0] == gene[0]]
	
	#gene starts before tad end, ends after start
	matches = tadChrSubset[(tadChrSubset[:,2] >= gene[1]) * (tadChrSubset[:,1] <= gene[2])]
	
	if len(matches) > 1:

		continue
	filteredGenes.append(gene)

filteredGenes = np.array(filteredGenes, dtype='object')


#for the non-disrupted, the patient does not matter, so keep all patient just in 1 array.

tadPositiveAndNegativeSet = []

tadInd = 0

#We define the positive/negative TADs based on ALL SVs. the negative set is only truly negative if there is also no other SV type disrupting it. 
for tad in tadDisruptions:

	if tadInd > 500:
		continue

	patientCount = dict()
		
	for sv in tadDisruptions[tad]:
		
		patient = sv[0]
		if patient not in patientCount:
			patientCount[patient] = []
		patientCount[patient].append(sv[1][8].svType)
	
	#get the patient names that have disrupted TADs
	patientsWithDisruptions = []
	for tadDisr in tadDisruptions[tad]:
		
		patientsWithDisruptions.append(tadDisr[0])

	#find the genes that are in this TAD
	
	splitTad = tad.split("_")
	geneChrMatches = filteredGenes[filteredGenes[:,0] == splitTad[0]]
	
	#the gene is inside the tad if the start is inside the tad, or if the end is inside the tad
	
	#is the start of the gene inside the tad?
	geneMatches = (geneChrMatches[:,1] <= int(splitTad[2])) * (geneChrMatches[:,2] >= int(splitTad[1]))
	
	allMatches = geneChrMatches[geneMatches]
	#add genes only when these are not affected by any mutation

	#go through all the genes and all the patients and their expression values.
	positivePatients = []
	negativePatients = []
	svTypes = []
	for gene in allMatches:
		
		#extract the row for this gene
		if gene[3].name not in expressionData[:,0]:
			continue

		geneExpr = expressionData[expressionData[:,0] == gene[3].name][0]
		
		#for each patient, append the expression to either the disrupted or non-disrupted based on the tad patient list
		
		for sv in tadDisruptions[tad]:

			svType = sv[1][8].svType

			patient = sv[0]
			
			if patient not in samples:
				continue
			patientInd = samples.index(patient)
			
			if patient not in positivePatients:
				positivePatients.append(patient)
				svTypes.append(sv[1][8].svType)
			
			#if this patient has multiple SVs disrupting this TAD, we do not know which one is causing an effect.
			#so we skip that TAD for this patient
			# if len(patientCount[patient]) > 1 and 'DUP' in patientCount[patient]:
			# 	continue
			
			#we filter genes in the TAD based on the SV type.
			if svType == 'DEL':
				#in case of DEL, we filter genes with any mutation. Deleted genes are not relevant

				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient] or \
				gene[3].name in svPatientsDup[patient]:

					continue

			elif svType == 'DUP':
				
				
				#in case of a DUP, we keep genes that are disrupted by the TAD disrupting DUP,
				#because those are the ones that see the effect.
				#because the CNV amp may overlap with the dup, ignore that one too. 
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in snvPatients[patient]:
					continue
				
			elif svType == 'INV':
				#only ignore genes that are in the INV.
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsDup[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient]:
					continue
				
			else:
				#if ITX, skip all mutations. 
				if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsItx[patient] or gene[3].name in cnvPatientsDel[patient] or \
				gene[3].name in cnvPatientsAmp[patient] or gene[3].name in snvPatients[patient] or \
				gene[3].name in svPatientsDup[patient]:
					continue
			
			disruptedPairs[patient][gene[3].name] = float(geneExpr[patientInd])
			

		for patientInd in range(1, len(samples)):
			
			#for the negative set, we consider all TADs that are not disrupted by
			#ANY SV.
			
			
			patient = samples[patientInd]
			
			
			#in the negative set, filter all genes that have ANY mutation. 
			if gene[3].name in svPatientsDel[patient] or gene[3].name in svPatientsInv[patient] or \
				gene[3].name in svPatientsDup[patient] or gene[3].name in svPatientsItx[patient] or \
				gene[3].name in cnvPatientsDel[patient] or gene[3].name in cnvPatientsAmp[patient] or \
				gene[3].name in snvPatients[patient]:
					continue
			
			#only if there is at least 1 gene without a mutation in this patient,
			#do we add the TAD as part of the negative set. 
			if patient not in patientsWithDisruptions:
				if patient not in negativePatients:
					negativePatients.append(patient)
	

			if patient not in patientsWithDisruptions:
				
				nonDisruptedPairs[gene[3].name][patient] = float(geneExpr[patientInd])
				
	#if len(otherSVs) > 0 and len(positivePatients) == 0:
	#	continue #so for this TAD, if it is not affected by e.g. a deletion, but it is by another SV type, then it should not be listed as specific for DEL. 
	tadPositiveAndNegativeSet.append([tad, positivePatients, negativePatients, svTypes])
	tadInd += 1

tadPositiveAndNegativeSet = np.array(tadPositiveAndNegativeSet, dtype='object')
np.savetxt(specificOutDir + '/tadPositiveAndNegativeSet.txt', tadPositiveAndNegativeSet, fmt='%s', delimiter='\t')

#For each gene in the disrupted group, compute the z-score of the gene compared to the expression of all patients in the negative group

#store the z-scores in an aggregated way, gene/patient pairs. Then we can get the overall plots. 
zScores = []
pValues = []
zScoresPerGene = dict()
scorePairs = [] #keep a list of the pairs that have a z-score originally, and which were set to 0 for the ranks. 
for patient in disruptedPairs:

	for gene in disruptedPairs[patient]:
		
		if gene not in zScoresPerGene:
			zScoresPerGene[gene] = dict()
			

		expr = disruptedPairs[patient][gene]

		#get the negative set
		if gene not in nonDisruptedPairs:
			continue
		negExprPatients = nonDisruptedPairs[gene]
		
		
		negExpr = []
		for negPatient in negExprPatients:
			
			negExpr.append(negExprPatients[negPatient])
		
		if np.std(negExpr) == 0:
			continue

		#compute the z-score
		z = (float(expr) - np.mean(negExpr)) / float(np.std(negExpr))
		
			
		posExpr = []
		for patientName in disruptedPairs:
			
			
			if gene in disruptedPairs[patientName]:
				posExpr.append(disruptedPairs[patientName][gene])
			
		zScoresPerGene[gene][patient] = z
		scorePairs.append(patient + '_' + gene)
		
#Go through the patients and compute the ranking of the z-scores inside the list.


for gene in zScoresPerGene:
	
	patients = list(zScoresPerGene[gene].keys())
	zScores = np.array(list(zScoresPerGene[gene].values()))

	#instead of ranks, we can also add a normalized/sigmoid value of the z-score. 

	for zScoreInd in range(0, len(zScores)):
		
		z = zScores[zScoreInd]
		patient = patients[zScoreInd]
		
		pValue = stats.norm.sf(abs(z))*2
		
		pValues.append([patient + '_' + gene, pValue, z])

pValues = np.array(pValues, dtype='object')

print(pValues.shape)

#Do MTC on the pValues
reject, pAdjusted, _, _ = multipletests(pValues[:,1], method='bonferroni') #fdr_bh or bonferroni
print(reject)
print(pAdjusted)
signPatients = []
for pValueInd in range(0, len(pValues[:,1])):

	signPatients.append([pValues[pValueInd][0], pValues[pValueInd][1], pAdjusted[pValueInd], reject[pValueInd], pValues[pValueInd][1], pValues[pValueInd][2]])

signPatients = np.array(signPatients, dtype='object')

print(signPatients.shape)
np.savetxt(specificOutDir + '/zScores.txt', signPatients, fmt='%s', delimiter='\t')




