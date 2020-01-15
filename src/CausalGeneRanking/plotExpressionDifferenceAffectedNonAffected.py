"""
	Get the expression of all genes in affected TADs, and in non-affected TADs.
	Show the difference between these groups in a boxplot. 
	
	
"""

import numpy as np
import matplotlib.pyplot as plt
import ast
from inputParser import InputParser
import sys
import settings
import random

#Load the affected and non-affected TAD data

tadPositiveAndNegativeSet = []
with open('tadPositiveAndNegativeSet.txt', 'r') as inF:
	
	for line in inF:
		
		splitLine = line.split('\t')
		tad = splitLine[0]
		positiveSet = ast.literal_eval(splitLine[1])
		negativeSet = ast.literal_eval(splitLine[2])
		svTypes = ast.literal_eval(splitLine[3])
		
		tadPositiveAndNegativeSet.append([tad, positiveSet, negativeSet, svTypes])

tadPositiveAndNegativeSet = np.array(tadPositiveAndNegativeSet, dtype='object')

#get the gene expression
expressionFile = sys.argv[1]

geneNameConversionMap = dict()
geneNameConversionFile = sys.argv[2]
with open(geneNameConversionFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		ensgId = splitLine[3]
		splitEnsgId = ensgId.split('.') #we only keep everything before the dot
		geneName = splitLine[4]
		geneNameConversionMap[splitEnsgId[0]] = geneName

expressionData = dict()
with open(expressionFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		if lineCount == 0:
			samples = line.split("\t")
			
			#get random indices to shuffle the patients
			shuffledInd = np.arange(0, len(samples))
			random.shuffle(shuffledInd)
			
			
			lineCount += 1
			continue
		splitLine = line.split("\t")
		fullGeneName = splitLine[0]
		if fullGeneName not in geneNameConversionMap:
			continue
		geneName = geneNameConversionMap[fullGeneName] #get the gene name rather than the ENSG ID

		data = np.array(splitLine)[1:]

		
		#shuffle expression.
		#data = data[shuffledInd]

		if geneName not in expressionData:
			expressionData[geneName] = dict()
		
		for sample in range(0, len(samples)):

			expressionData[geneName][samples[sample]] = data[sample]
			
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#exclude mutated patients
### change this path later!!!
svPatients = np.load('/hpc/cog_bioinf/ridder/users/mnieboer/data/dataProcessing/integrative-omics/src/CausalGeneRanking/svPatients.npy', allow_pickle=True, encoding='latin1').item()
snvPatients = np.load('/hpc/cog_bioinf/ridder/users/mnieboer/data/dataProcessing/integrative-omics/src/CausalGeneRanking/snvPatients.npy', allow_pickle=True, encoding='latin1').item()
cnvPatients = np.load('/hpc/cog_bioinf/ridder/users/mnieboer/data/dataProcessing/integrative-omics/src/CausalGeneRanking/cnvPatients.npy', allow_pickle=True, encoding='latin1').item()

#load the combinations that we found in rule-based, as an extra filter.
ruleBasedCombinations = np.loadtxt('Output/RankedGenes/0/BRCA/nonCoding_geneSVPairs.txt_', dtype='object')

	
###for plotting just expression
#go through each TAD, and divide the gene expression into the 2 respective groups.
# affectedExpression = []
# nonAffectedExpression = []
# svType = sys.argv[3]
# for tad in tadPositiveAndNegativeSet:
# 	
# 	#if svType not in tad[3]: #skip this that if no patient has this sv type at all
# 	#	continue
# 	
# 	splitTad = tad[0].split('_')
# 	
# 	#get the genes that are inside this TAD
# 	geneChrSubset = allGenes[allGenes[:,0] == splitTad[0]]
# 	
# 	genes = geneChrSubset[(geneChrSubset[:,2] >= int(splitTad[1])) * (geneChrSubset[:,1] <= int(splitTad[2]))]
# 	
# 	if len(tad[1]) > 0: #this TAD is affected
# 
# 		patientInd = -1
# 		for affectedPatient in tad[1]:
# 			patientInd += 1
# 			#check if the SV type of the patient matches the SV type restriction
# 			#if tad[3][patientInd] != svType:
# 			#	continue
# 			
# 			#add the gene expression values of this patient to the array. 
# 			for gene in genes:
# 				
# 				if gene[3].name not in expressionData:
# 					continue
# 				
# 				if gene[3].name in svPatients[affectedPatient] or gene[3].name in snvPatients[affectedPatient] or gene[3].name in cnvPatients[affectedPatient]:
# 					continue
# 					
# 				
# 				expression = expressionData[gene[3].name][affectedPatient]
# 				affectedExpression.append(float(expression))
# 	else: 
# 		
# 		for negativePatient in tad[2]:
# 			
# 			for gene in genes:
# 				
# 				if gene[3].name not in expressionData:
# 					continue
# 				
# 				if gene[3].name in svPatients[negativePatient] or gene[3].name in snvPatients[negativePatient] or gene[3].name in cnvPatients[negativePatient]:
# 					continue
# 				
# 				expression = expressionData[gene[3].name][negativePatient]
# 				nonAffectedExpression.append(float(expression))	
# 
# print(len(affectedExpression))
# print(len(nonAffectedExpression))
# 		
# data = [affectedExpression, nonAffectedExpression]
# 
# plt.boxplot(data)
# plt.show()
# 
# affectedExpression = np.array(affectedExpression)
# nonAffectedExpression = np.array(nonAffectedExpression)
# 
# #filter the data
# filteredAffected = affectedExpression[affectedExpression < np.percentile(affectedExpression,95)]
# filteredNonAffected = nonAffectedExpression[nonAffectedExpression < np.percentile(nonAffectedExpression,95)]
# 
# filteredData = [filteredAffected, filteredNonAffected]
# plt.boxplot(filteredData)
# plt.show()
# exit()
### for plotting z-scores compared to non-affected
#go through each TAD, and divide the gene expression into the 2 respective groups.
svType = sys.argv[3]
affectedZ = []
for tad in tadPositiveAndNegativeSet:
	
	if svType not in tad[3]: #skip this that if no patient has this sv type at all
		continue
	
	splitTad = tad[0].split('_')
	
	#get the genes that are inside this TAD
	geneChrSubset = allGenes[allGenes[:,0] == splitTad[0]]
	
	genes = geneChrSubset[(geneChrSubset[:,2] >= int(splitTad[1])) * (geneChrSubset[:,1] <= int(splitTad[2]))]
	
	if len(tad[1]) > 0: #this TAD is affected
		
		#first get the negative set for each gene.
		for gene in genes:
			
			if gene[3].name not in expressionData:
				continue
			
			negativeExpr = []
			for negativePatient in tad[2]:
				
				if gene[3].name in svPatients[negativePatient] or gene[3].name in snvPatients[negativePatient] or gene[3].name in cnvPatients[negativePatient]:
					continue
				
				expression = expressionData[gene[3].name][negativePatient]
				negativeExpr.append(float(expression))
			
			if len(negativeExpr) < 1: #stop if we have no patients in the negative set
				continue
			
			#go through each affected patient in this TAD, and compute the z-score to the negative set.
			patientInd = -1
			for affectedPatient in tad[1]:
				patientInd += 1
				
				#check if the SV type of the patient matches the SV type restriction
				if tad[3][patientInd] != svType:
					continue
				
				if gene[3].name in svPatients[affectedPatient] or gene[3].name in snvPatients[affectedPatient] or gene[3].name in cnvPatients[affectedPatient]:
					continue
				
				#if gene[3].name not in cnvPatients[affectedPatient]:
				#	continue
				
				expression = expressionData[gene[3].name][affectedPatient]
				
				if np.std(negativeExpr) == 0: #stop if we cannot compute a z-score
					continue
				
				z = (float(expression) - np.mean(negativeExpr)) / float(np.std(negativeExpr))
				affectedZ.append(z)
		
plt.boxplot(affectedZ)
plt.show()

affectedZ = np.array(affectedZ)

#filter the data
filteredZ = affectedZ[affectedZ < np.percentile(affectedZ,95)]

plt.boxplot(filteredZ)
plt.ylim([-2,7])
plt.show()
# 
# ##filtering out pairs not identified with rule-based. What do we see?
# ruleBasedPairs = []
# for combination in ruleBasedCombinations:
# 	
# 	splitPair = combination[0].split('_')
# 	
# 	ruleBasedPairs.append(splitPair[0] + '_' + splitPair[7])
# 	
# 
# affectedZ = []
# for tad in tadPositiveAndNegativeSet:
# 	
# 	splitTad = tad[0].split('_')
# 	
# 	#get the genes that are inside this TAD
# 	geneChrSubset = allGenes[allGenes[:,0] == splitTad[0]]
# 	
# 	genes = geneChrSubset[(geneChrSubset[:,2] >= int(splitTad[1])) * (geneChrSubset[:,1] <= int(splitTad[2]))]
# 	
# 	if len(tad[1]) > 0: #this TAD is affected
# 		
# 		#first get the negative set for each gene.
# 		for gene in genes:
# 			
# 			if gene[3].name not in expressionData:
# 				continue
# 			
# 			negativeExpr = []
# 			for negativePatient in tad[2]:
# 				
# 				if gene[3].name in svPatients[negativePatient] or gene[3].name in snvPatients[negativePatient] or gene[3].name in cnvPatients[negativePatient]:
# 					continue
# 				
# 				expression = expressionData[gene[3].name][negativePatient]
# 				negativeExpr.append(float(expression))
# 			
# 			if len(negativeExpr) < 1: #stop if we have no patients in the negative set
# 				continue
# 			
# 			#go through each affected patient in this TAD, and compute the z-score to the negative set. 
# 			for affectedPatient in tad[1]:
# 				
# 				if gene[3].name + '_' + affectedPatient not in ruleBasedPairs:
# 					continue
# 				
# 				if gene[3].name in svPatients[affectedPatient] or gene[3].name in snvPatients[affectedPatient] or gene[3].name in cnvPatients[affectedPatient]:
# 					continue
# 				
# 				expression = expressionData[gene[3].name][affectedPatient]
# 				
# 				if np.std(negativeExpr) == 0: #stop if we cannot compute a z-score
# 					continue
# 				
# 				z = (float(expression) - np.mean(negativeExpr)) / float(np.std(negativeExpr))
# 				affectedZ.append(z)
# 	
# 		
# plt.boxplot(affectedZ)
# plt.show()
# 
# affectedZ = np.array(affectedZ)
# 
# #filter the data
# filteredZ = affectedZ[affectedZ < np.percentile(affectedZ,95)]
# 
# plt.boxplot(filteredZ)
# plt.show()
