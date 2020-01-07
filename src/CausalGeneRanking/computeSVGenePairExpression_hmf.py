
from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import re
from scipy import stats
from os import listdir
from os.path import isfile, join
from six.moves import range
import glob
from scipy.stats import rankdata
from inputParser import InputParser
import settings

# Get the expression z-scores for every SV. 

nonCodingPairs = np.loadtxt(sys.argv[1], dtype="object")
codingPairs = np.loadtxt(sys.argv[2], dtype="object")

#Get the SNV data for these patients, make sure to map to the same identifiers
snvDir = sys.argv[3]

#use this conversion map to get the same gene names as output by the prioritizer
geneNameConversionMap = dict()
geneNameConversionFile = sys.argv[4]
with open(geneNameConversionFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		ensgId = splitLine[3]
		splitEnsgId = ensgId.split('.') #we only keep everything before the dot
		geneName = splitLine[4]
		geneNameConversionMap[splitEnsgId[0]] = geneName

import gzip
#search through the SNVs and link these to genes.
vcfs = glob.glob(snvDir + '/**/*.somatic.vcf.gz', recursive=True)

variantsList = []
addedVariants = [] #check if based on pairs no duplicates are added. 
mutations = []
for vcf in vcfs:
	print(vcf)
	
	#get the samplename from the vcf
	sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
	
	#open the .gz file
	with gzip.open(vcf, 'rb') as inF:
		
		for line in inF:
			line = line.strip().decode('utf-8')

			if re.search('^#', line): #skip header
				continue
			
			#skip the SV if it did not pass.
			splitLine = line.split("\t")
			filterInfo = splitLine[6]
			if filterInfo != 'PASS':
				continue
	
			#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
			infoField = splitLine[7]
			
			geneSearch = re.search('(ENSG\d+)', infoField)
			if geneSearch:
				geneMatch = re.search('(ENSG\d+)', infoField).group(1)
				if geneMatch not in geneNameConversionMap: #skip genes for which we do not know the name
					continue
				geneName = geneNameConversionMap[geneMatch]
				
				mutations.append([geneName, sampleName])

mutations = np.array(mutations, dtype="object")
print(mutations)

cnvs = []
tsvs = glob.glob(snvDir + '/**/*.gene.tsv', recursive=True)

for tsv in tsvs:
	
	#get the samplename from the vcf
	sampleName = re.search('.*\/([A-Z\d]+)\.', tsv).group(1)
	
	#open the .gz file
	with open(tsv, 'r') as inF:
		
		lineCount = 0
		for line in inF:

			if lineCount < 1: #skip header
				lineCount += 1
				continue

			splitLine = line.split("\t")
			
			gene = splitLine[3]
			
			if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
				continue

			cnvs.append([gene, sampleName])

cnvs = np.array(cnvs, dtype='object')

geneCheck = 'ACAD8'
patientCheck = 'CPCT02050143T'

expressionFile = sys.argv[5]

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
		if fullGeneName not in geneNameConversionMap:
			continue
		geneName = geneNameConversionMap[fullGeneName] #get the gene name rather than the ENSG ID

		
		data = splitLine[1:len(splitLine)] 
		fixedData = [geneName]
		fixedData += data
		expressionData.append(fixedData)

expressionData = np.array(expressionData, dtype="object")	

#Get the z-scores for every pair

#For every gene, get a list of all samples in which this gene is affected to exclude these and make a null distribution
geneSampleRef = dict()
for pair in nonCodingPairs[:,0]:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[7]
	
	if gene not in geneSampleRef:
		geneSampleRef[gene] = []
	geneSampleRef[gene].append(sample)

for pair in codingPairs:
	splitPair = pair.split("_")
	gene = splitPair[0]
	sample = splitPair[7]

	if gene not in geneSampleRef:
		geneSampleRef[gene] = []
	geneSampleRef[gene].append(sample)

for pair in mutations:

	if pair[0] not in geneSampleRef:
		geneSampleRef[pair[0]] = []
	geneSampleRef[pair[0]].append(pair[1])

for pair in cnvs:

	if pair[0] not in geneSampleRef:
		geneSampleRef[pair[0]] = []
	geneSampleRef[pair[0]].append(pair[1])


#Set for every gene the expression values in all possible samples for lookup
geneSampleExpr = dict()
for gene in geneSampleRef:
	
	
	if gene not in expressionData[:,0]:
		continue
	
	geneSamples = geneSampleRef[gene]

	geneSampleExpr[gene] = dict()
	geneExpression = expressionData[expressionData[:,0] == gene][0]
	for geneSample in geneSamples:
		
		#get the expression data for this sample
		if geneSample not in samples: #this sample has no RNA-seq data
			continue
		sampleInd = samples.index(geneSample)
		geneSampleExpr[gene][geneSample] = float(geneExpression[sampleInd])
		
			
				
print("done getting expr for samples")

#Also set the negative set for every gene consisting of the expression of all samples wthout any SV
negativeExpr = dict()
for gene in geneSampleExpr:
	matchedFullSampleNames = list(geneSampleExpr[gene].keys())
	
	#Get all the samples without an SV for this gene
	unmatchedSamples = np.setdiff1d(samples[1:len(samples)], matchedFullSampleNames) #exclude hybrid ref
	negativeSamples = []
	for sample in unmatchedSamples: #sample tumor samples, exclude normals
		
		negativeSamples.append(sample)
	
	geneExpression = expressionData[expressionData[:,0] == gene][0]	
	#Get the expression of these samples
	negativeSampleExpressionValues = []
	for sample in negativeSamples:
		if sample not in samples: #this sample has no RNA-seq data
			continue
		sampleInd = samples.index(sample)
		
		if gene == 'SEPT2':
			print('neg sample: ', sample)
			print(float(geneExpression[sampleInd]))
			
		negativeSampleExpressionValues.append(float(geneExpression[sampleInd]))
	
	negativeExpr[gene] = negativeSampleExpressionValues
exit()	
	
print("negative expr done")

def getDEPairs(pairs, geneSampleRef, epressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr):
	
	zScores = dict()							 	
	for pair in pairs:
		splitPair = pair.split("_")
		gene = splitPair[0]
		pairSample = splitPair[7]

		sv = "_".join(splitPair[1:])
		if gene not in expressionData[:,0]:
			continue
		
		if pairSample not in geneSampleExpr[gene]:
			continue
		
		sampleExpressionValue = geneSampleExpr[gene][pairSample] #expression values of this gene in all samples
		matchedFullSampleNames = list(geneSampleExpr[gene].keys())
					
		
		negativeSampleExpressionValues = negativeExpr[gene]
		
		#Get the expression z-score for this pair
		if np.std(negativeSampleExpressionValues) == 0:
			continue
	
		z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
		pValue = stats.norm.sf(abs(z))*2
		
		if gene == 'SEPT2' and pairSample == 'CPCT02030271T':
			print('expr: ', sampleExpressionValue)
			print('z:', z)
			print('ean neg: ', np.mean(negativeSampleExpressionValues))
			print('std neg: ', np.std(negativeSampleExpressionValues))
			import matplotlib.pyplot as plt
			plt.boxplot(negativeSampleExpressionValues)
			plt.show()
			exit()
		
	
		perPairDifferentialExpression[pair] = pValue
		zScores[pair] = z
	return perPairDifferentialExpression, zScores

def getDEPairsSNVs(pairs, geneSampleRef, epressionData, perPairDifferentialExpression, geneSampleExpr, negativeExpr):
									
	for pair in pairs:
		
		gene = pair[0]
		pairSample = pair[1]
		
		if gene not in expressionData[:,0]:
			continue
		
		if pairSample not in geneSampleExpr[gene]: #sometimes there is no expr data for that sample
			continue 
		
		sampleExpressionValue = geneSampleExpr[gene][pairSample] #expression values of this gene in all samples
		matchedFullSampleNames = list(geneSampleExpr[gene].keys())
					
		negativeSampleExpressionValues = negativeExpr[gene]
		
		#Get the expression z-score for this pair
		if np.std(negativeSampleExpressionValues) == 0:
			continue
	
		z = (sampleExpressionValue - np.mean(negativeSampleExpressionValues)) / float(np.std(negativeSampleExpressionValues))
		pValue = stats.norm.sf(abs(z))*2
		
		
		
	
		perPairDifferentialExpression[pair[0] + "_" + pair[1]] = pValue
		
	return perPairDifferentialExpression

# Output DEG pairs for non-coding only
from statsmodels.sandbox.stats.multicomp import multipletests
perPairDifferentialExpression, zScores = getDEPairs(nonCodingPairs[:,0], geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print("done")

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 5], dtype="object")
perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())

reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArray[:,2] = reject
perPairDifferentialExpressionArray[:,3] = list(zScores.values())

np.savetxt(sys.argv[1] + '_nonCodingPairsZScores.txt', perPairDifferentialExpressionArray, delimiter='\t', fmt='%s')

#Then convert the z-scores to ranks.

#The rank should be within the patient.
#so first gather all patients, then all genes.
genes = []
zScoresPerGene = dict()
scorePairs = []
for pair in perPairDifferentialExpressionArray:
	splitPair = pair[0].split('_')
	
	gene = splitPair[0]
	patient = splitPair[7]
	genes.append(gene)
	
	if gene not in zScoresPerGene:
		zScoresPerGene[gene] = dict()
	zScoresPerGene[gene][patient] = pair[3]
	scorePairs.append(patient + '_' + gene)

#add remaining genes that we did not link to SVs.
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

for gene in allGenes:
	if gene[3].name not in zScoresPerGene:
		zScoresPerGene[gene[3].name] = dict()

for gene in zScoresPerGene:
	for patient in samples:
		if patient == '':
			continue

		if patient not in zScoresPerGene[gene]:
			zScoresPerGene[gene][patient] = 0
		

#Get the ranks
ranks = []
for gene in zScoresPerGene:
	
	patients = list(zScoresPerGene[gene].keys())
	zScores = np.array(list(zScoresPerGene[gene].values()))

	rankedZScores = rankdata(zScores, method='min')
	
	#normalize rans between -1 and 1
	minVal = np.min(rankedZScores)
	maxVal = np.max(rankedZScores)
	
	normalizedRanks = 2 * ((rankedZScores - minVal) / (maxVal - minVal)) - 1
	
	for zScoreInd in range(0, len(zScores)):
		omitted = True
		if patients[zScoreInd] + '_' + gene in scorePairs:
			omitted = False
		ranks.append([patients[zScoreInd] + '_' + gene, zScores[zScoreInd], normalizedRanks[zScoreInd], omitted])


#np.save(sys.argv[1] + '_nonCodingPairDEGs.npy', perPairDifferentialExpressionArray)
np.savetxt(sys.argv[1] + '_nonCodingPairsRanks.txt', ranks, delimiter='\t', fmt='%s')
exit()

perPairDifferentialExpression = getDEPairs(codingPairs, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print("done")

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.save(sys.argv[1] + '_codingPairDEGs.npy', perPairDifferentialExpressionArrayFiltered)

#finally repeat for SNVs as well

perPairDifferentialExpression = getDEPairsSNVs(mutations, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print("done")

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.save(sys.argv[1] + '_codingSNVDEGs.npy', perPairDifferentialExpressionArrayFiltered)

#CNVs

perPairDifferentialExpression = getDEPairsSNVs(cnvs, geneSampleRef, expressionData, dict(), geneSampleExpr, negativeExpr)
print("done")

perPairDifferentialExpressionArray = np.empty([len(perPairDifferentialExpression), 2], dtype="object")
perPairDifferentialExpressionArray[:,0] = list(perPairDifferentialExpression.keys())
perPairDifferentialExpressionArray[:,1] = list(perPairDifferentialExpression.values())

from statsmodels.sandbox.stats.multicomp import multipletests
reject, pAdjusted, _, _ = multipletests(perPairDifferentialExpressionArray[:,1], method='bonferroni')

perPairDifferentialExpressionArrayFiltered = perPairDifferentialExpressionArray[reject]

np.save(sys.argv[1] + '_codingCNVDEGs.npy', perPairDifferentialExpressionArrayFiltered)

