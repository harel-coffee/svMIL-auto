"""
	
	First make a plot of just the expression in TADs for which we have a non-coding SV disrupting the boundary

"""

import settings
from inputParser import InputParser
import numpy as np
import sys
import glob
import re
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from genomicShuffler import GenomicShuffler

#permutationRound = sys.argv[4]

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
geneNameConversionFile = sys.argv[1]
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
expressionFile = sys.argv[2]

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
print(expressionData)

# 
# #Get all SVs
# svDir = settings.files['svDir']
# svData = InputParser().getSVsFromFile_hmf(svDir)
# 
# #Filter out the coding effect SVs, we want to focus on non-coding SVs. 
# excludedSVs = np.loadtxt(settings.files['excludedSVs'], dtype='object')
# 
# 
# filteredSVs = []
# types = []
# for sv in svData:
# 
# 	svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName
# 	#if svEntry not in excludedSVs:
# 	#	filteredSVs.append(sv)
# 		
# 	#	if sv[8].svType not in types:
# 	#		types.append(sv[8].svType)
# 	filteredSVs.append(sv)
# 
# print(types)	
# 
# filteredSVs = np.array(filteredSVs, dtype='object')
# 
# np.save('filteredSVs_' + svType + '.npy', filteredSVs)

filteredSVs = np.load('filteredSVs.npy', allow_pickle=True, encoding='latin1')
print(filteredSVs.shape)


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


#to shuffle across patients, first transpose, the shuffle, then transpose back.
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
# #shuffledExpressionData = np.concatenate((genes, shuffledExpression), axis=1)
# 
# print(samples)
# print(shuffledExpressionData)
# 
# expressionData = shuffledExpressionData
# expressionDataRavel = expressionData[:,1:].ravel()
# np.random.shuffle(expressionDataRavel)
# expression = expressionDataRavel.reshape(expressionData[:,1:].shape)
# shuffledExpressionData = np.empty(expressionData.shape, dtype='object')
# shuffledExpressionData[:,0] = expressionData[:,0]
# shuffledExpressionData[:,1:] = expression
# expressionData = shuffledExpressionData

#Get a list of all genes * patients that have mutations.
# 
# #Get the CNVs per gene
# def getPatientsWithCNVGeneBased(cnvDir):
# 	
# 	tsvs = glob.glob(cnvDir + '/**/*.gene.tsv', recursive=True)
# 	
# 	cnvPatientsDel = dict()
# 	cnvPatientsAmp = dict()
# 	
# 	for tsv in tsvs:
# 		
# 		#get the samplename from the vcf
# 		sampleName = re.search('.*\/([A-Z\d]+)\.', tsv).group(1)
# 		
# 		if sampleName not in cnvPatientsAmp:
# 			cnvPatientsAmp[sampleName] = []
# 		if sampleName not in cnvPatientsDel:
# 			cnvPatientsDel[sampleName] = []
# 		
# 		#open the .gz file
# 		with open(tsv, 'r') as inF:
# 			
# 			lineCount = 0
# 			for line in inF:
# 
# 				if lineCount < 1: #skip header
# 					lineCount += 1
# 					continue
# 	
# 				splitLine = line.split("\t")
# 				
# 				gene = splitLine[3]
# 				
# 				
# 				if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
# 					continue
# 				
# 				if float(splitLine[5]) > 2.3: 
# 					
# 					cnvPatientsAmp[sampleName].append(gene)
# 				elif float(splitLine[5]) < 1.7:
# 					
# 					cnvPatientsDel[sampleName].append(gene)
# 					
# 	return cnvPatientsAmp, cnvPatientsDel
# 
# cnvPatientsAmp, cnvPatientsDel = getPatientsWithCNVGeneBased(sys.argv[3])
# 
# #Get the SNVs per gene
# 
# def getPatientsWithSNVs(snvDir):
# 	import gzip
# 	#search through the SNVs and link these to genes.
# 	vcfs = glob.glob(snvDir + '/**/*.somatic.vcf.gz', recursive=True)
# 
# 	patientsWithSNVs = dict()
# 	for vcf in vcfs:
# 		
# 		#get the samplename from the vcf
# 		sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
# 		
# 		
# 		#open the .gz file
# 		with gzip.open(vcf, 'rb') as inF:
# 			
# 			for line in inF:
# 				line = line.strip().decode('utf-8')
# 
# 				if re.search('^#', line): #skip header
# 					continue
# 				
# 				#skip the SV if it did not pass.
# 				splitLine = line.split("\t")
# 				filterInfo = splitLine[6]
# 				if filterInfo != 'PASS':
# 					continue
# 		
# 				#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
# 				infoField = splitLine[7]
# 				
# 				geneSearch = re.search('(ENSG\d+)', infoField)
# 				if geneSearch:
# 					geneMatch = re.search('(ENSG\d+)', infoField).group(1)
# 					#skip genes for which we do not know the name
# 					if geneMatch not in geneNameConversionMap:
# 						continue
# 					geneName = geneNameConversionMap[geneMatch]
# 					
# 					if sampleName not in patientsWithSNVs:
# 						patientsWithSNVs[sampleName] = []
# 					patientsWithSNVs[sampleName].append(geneName)
# 
# 	return patientsWithSNVs
# 
# snvPatients = getPatientsWithSNVs(sys.argv[3])
# 
# # #Get the SV per gene
# def getPatientsWithSVs(svDir, allGenes):
# 	
# 	#Get all parsed and annotated SV type files from the main dir
# 	
# 	vcfs = glob.glob(svDir + '/**/*.svTypes.passed', recursive=True)
# 	
# 	svPatientsDel = dict()
# 	svPatientsDup = dict()
# 	svPatientsInv = dict()
# 	svPatientsItx = dict()
# 	
# 	for vcf in vcfs:
# 		
# 		#get the samplename from the vcf
# 		sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
# 		if sampleName not in svPatientsDel:
# 			svPatientsDel[sampleName] = []
# 		if sampleName not in svPatientsDup:
# 			svPatientsDup[sampleName] = []
# 		if sampleName not in svPatientsInv:
# 			svPatientsInv[sampleName] = []
# 		if sampleName not in svPatientsItx:
# 			svPatientsItx[sampleName] = []	
# 
# 
# 		#open the .gz file
# 		with open(vcf, 'r') as inF:
# 			
# 			for line in inF:
# 
# 				if re.search('^#', line): #skip header
# 					continue
# 				
# 				#skip the SV if it did not pass.
# 				splitLine = line.split("\t")
# 				filterInfo = splitLine[6]
# 				if filterInfo != 'PASS':
# 					continue
# 				
# 				#Check if the SV is a deletion
# 				infoField = splitLine[7]
# 				splitInfoField = infoField.split(";")
# 				svType = ''
# 				for field in splitInfoField:
# 					
# 					splitField = field.split("=")
# 					if splitField[0] == 'SIMPLE_TYPE':
# 						svType = splitField[1]
# 				
# 				#skip non-deletions
# 				if svType not in ['DEL', 'DUP', 'INV', 'ITX']:
# 				#if svType not in ['DUP']:
# 					continue
# 				
# 				chr1 = splitLine[0]
# 				pos1 = int(splitLine[1])
# 				pos2Info = splitLine[4]
# 				pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))
# 				chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)
# 				
# 				s1 = pos1
# 				e1 = pos1
# 				s2 = pos2
# 				e2 = pos2
# 				orderedChr1 = chr1
# 				orderedChr2 = chr2
# 				
# 				#switch chromosomes if necessary
# 				if chr1 != chr2:
# 					if chr1 == 'Y' and chr2 == 'X':
# 						orderedChr1 = chr2
# 						orderedChr2 = chr1
# 					if (chr1 == 'X' or chr1 == 'Y' or chr1 == 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
# 						orderedChr1 = chr2
# 						orderedChr2 = chr1
# 					if (chr1 != 'X' and chr1 != 'Y' and chr1 != 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
# 						if int(chr1) > int(chr2):
# 							orderedChr1 = chr2
# 							orderedChr2 = chr1
# 					if (chr1 in ['X', 'Y', 'MT']) and (chr2 in ['X', 'Y', 'MT']): #order these as well
# 						if chr1 == 'Y' and chr2 == 'X':
# 							orderedChr1 = chr2
# 							orderedChr2 = chr1
# 						if chr1 == 'MT' and chr2 in ['X', 'Y']:
# 							orderedChr1 = chr2
# 							orderedChr2 = chr1
# 
# 					
# 					#always switch the coordinates as well if chromosomes are switched.
# 					if orderedChr1 == chr2:
# 						s1 = pos2
# 						e1 = pos2
# 						s2  = pos1
# 						e2 = pos1	
# 				
# 				else: #if the chr are the same but the positions are reversed, change these as well. 
# 					if pos2 < pos1:
# 						s1 = pos2
# 						e1 = pos2
# 						s2  = pos1
# 						e2 = pos1	
# 
# 				chr1 = 'chr' + orderedChr1
# 				chr2 = 'chr' + orderedChr2
# 			
# 				#Check which genes are overlapped by this SV.
# 				#Keep track of the disrupted genes in the patient. 
# 				
# 				#intrachromosomal SV
# 				if chr1 == chr2:
# 					
# 					geneChrSubset = allGenes[allGenes[:,0] == chr1]
# 					
# 					geneMatches = geneChrSubset[(geneChrSubset[:,1] <= e2) * (geneChrSubset[:,2] >= s1)]
# 					
# 					if svType == 'DEL':
# 						for match in geneMatches:
# 							svPatientsDel[sampleName].append(match[3].name)
# 							
# 							
# 					elif svType == 'DUP':
# 						for match in geneMatches:
# 							svPatientsDup[sampleName].append(match[3].name)
# 					elif svType == 'INV':
# 						for match in geneMatches:
# 							svPatientsInv[sampleName].append(match[3].name)
# 						
# 				else:
# 					
# 					#find breakpoints in the gene for each side of the SV
# 					geneChr1Subset = allGenes[allGenes[:,0] == chr1]
# 					geneChr2Subset = allGenes[allGenes[:,0] == chr2]
# 					
# 					#check if the bp start is within the gene. 
# 					geneChr1Matches = geneChr1Subset[(s1 >= geneChr1Subset[:,1]) * (s1 <= geneChr1Subset[:,2])]
# 					geneChr2Matches = geneChr2Subset[(s2 >= geneChr2Subset[:,1]) * (s2 <= geneChr2Subset[:,2])]
# 					
# 					for match in geneChr1Matches:
# 						svPatientsItx[sampleName].append(match[3].name)
# 						
# 						
# 						
# 					for match in geneChr2Matches:
# 						svPatientsItx[sampleName].append(match[3].name)
# 						
# 						
# 			
# 	return svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx
# 
# svPatientsDel, svPatientsDup, svPatientsInv, svPatientsItx = getPatientsWithSVs(sys.argv[3], allGenes)

# np.save('svPatientsDel.npy', svPatientsDel)
# np.save('svPatientsDup.npy', svPatientsDup)
# np.save('svPatientsInv.npy', svPatientsInv)
# np.save('svPatientsItx.npy', svPatientsItx)
#np.save('cnvPatientsDel.npy', cnvPatientsDel)
#np.save('cnvPatientsAmp.npy', cnvPatientsAmp)

# np.save('svPatients_dupOnly.npy', svPatients)
# np.save('cnvPatients_ampOnly.npy', cnvPatients)
# np.save('snvPatients.npy', snvPatients)
# exit()
#svPatients = np.load('svPatients.npy', allow_pickle=True, encoding='latin1').item()
snvPatients = np.load('snvPatients.npy', allow_pickle=True, encoding='latin1').item()
#cnvPatients = np.load('cnvPatients.npy', allow_pickle=True, encoding='latin1').item()

svPatientsDel = np.load('svPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load('svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svPatientsInv = np.load('svPatientsInv.npy', allow_pickle=True, encoding='latin1').item()
svPatientsItx = np.load('svPatientsItx.npy', allow_pickle=True, encoding='latin1').item()

cnvPatientsDel = np.load('cnvPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
cnvPatientsAmp = np.load('cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()


disruptedTadExpression = []
nonDisruptedTadExpression = []

disruptedPairs = dict() #store for each patient a dict with genes, and in there, the expression of the gene. 
nonDisruptedPairs = dict()

for patient in samples:
	
	if patient == '':
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

	if tad != 'chr8_113800000_115100000' and tad != 'chr8_118325000_119725000':
		continue
	print(tad)	
	
	patientCount = dict()
		
	for sv in tadDisruptions[tad]:
		
		patient = sv[0]
		if patient not in patientCount:
			patientCount[patient] = []
		patientCount[patient].append(sv[1][8].svType)
	
	# for patient in patientCount:
	# 	
	# 	if len(patientCount[patient]) > 1:
	# 		print(tad)
	# 		print(patient)
	# 		print(patientCount[patient])
	# 		exit()
	# continue
	
	
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
		
		print(gene[3].name)
		print(gene)
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
			
			if patient != 'CPCT02050146T':
				continue
			
			print(sv)
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
			
			print('adding')
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
exit()
tadPositiveAndNegativeSet = np.array(tadPositiveAndNegativeSet, dtype='object')
np.savetxt('tadPositiveAndNegativeSet.txt', tadPositiveAndNegativeSet, fmt='%s', delimiter='\t')

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
		
		if z > 8:
			print(gene)
			print(patient)
			print(z)
			print(expr)
			print(negExpr)
			print(np.mean(negExpr))
			print(np.std(negExpr))
			
			#print overall distribution of this gene
		posExpr = []
		for patientName in disruptedPairs:
			
			
			if gene in disruptedPairs[patientName]:
				posExpr.append(disruptedPairs[patientName][gene])
			
		zScoresPerGene[gene][patient] = z
		scorePairs.append(patient + '_' + gene)
		
		#pValue = stats.norm.sf(abs(z))*2
		#pValues.append([patient + '_' + gene, pValue, expr, np.mean(negExpr), np.std(negExpr)])


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

#np.savetxt('tadDisr/zScores_random_degs_' + str(permutationRound) + '.txt', zScores, fmt='%s', delimiter='\t')
#np.savetxt('tadDisr/pValues_shuffled_' + str(permutationRound) + '.txt', signPatients, fmt='%s', delimiter='\t')
np.savetxt('pValues_allGenes_smallTads_dupFilter.txt', signPatients, fmt='%s', delimiter='\t')
#np.savetxt('pValues_shuffled.txt', signPatients, fmt='%s', delimiter='\t')
#np.savetxt('pValues.txt', signPatients, fmt='%s', delimiter='\t')
exit()
# import matplotlib.pyplot as plt
# 
# zScores = np.loadtxt('zScores.txt', dtype='object')
# 
# z = zScores[(zScores[:,1] != 'nan') * (zScores[:,1] != 'inf')]
# z = [float(i) for i in z[:,1]]
# 
# 
# plt.hist(np.log(z))
# plt.show()


			
#print(len(disruptedTadExpression))
#print(len(nonDisruptedTadExpression))

#disruptedTadExpression = np.array(disruptedTadExpression, dtype='float')
#nonDisruptedTadExpression = np.array(nonDisruptedTadExpression, dtype='float')

#np.savetxt('disruptedTadExpression_del.txt', disruptedTadExpression, fmt='%s', delimiter='\t')
#np.savetxt('nonDisruptedTadExpression_del.txt', nonDisruptedTadExpression, fmt='%s', delimiter='\t')

import matplotlib.pyplot as plt

zScores = np.loadtxt('zScores.txt', dtype='object')
zScoresRandom = np.loadtxt('zScores_random_degs.txt', dtype='object')

zScores = zScores[zScores[:,1] != 'inf']
zScoresRandom = zScoresRandom[zScoresRandom[:,1] != 'inf']

zScores = zScores[:,1].astype(float)
zScoresRandom = zScoresRandom[:,1].astype(float)
plt.hist(np.log(zScores))
plt.show()
plt.clf()
plt.hist(np.log(zScoresRandom))
plt.show()
plt.clf()

exit()

disruptedTadExpression = np.loadtxt('disruptedTadExpression_itx.txt', dtype='float')
nonDisruptedTadExpression = np.loadtxt('nonDisruptedTadExpression_itx.txt', dtype='float')


# 		
# plt.boxplot(disruptedTadExpression)
# plt.show()
# plt.clf()
# plt.boxplot(nonDisruptedTadExpression)
# plt.show()
# 		
# plt.boxplot(np.log(disruptedTadExpression))
# plt.show()
# plt.clf()
# plt.boxplot(np.log(nonDisruptedTadExpression))
# plt.show()

filteredDisruptedExpression = disruptedTadExpression[disruptedTadExpression != 0]
filteredNonDisruptedExpression = nonDisruptedTadExpression[nonDisruptedTadExpression != 0]

# plt.hist(np.log(filteredDisruptedExpression))
# plt.show()
# plt.clf()
# plt.hist(np.log(filteredNonDisruptedExpression))
# plt.show()
# plt.clf()

plt.boxplot(np.log(filteredDisruptedExpression))
plt.show()
plt.clf()
plt.boxplot(np.log(filteredNonDisruptedExpression))
plt.show()

def tTest(data1,data2):
	
	z = (np.mean(data1) - np.mean(data2)) / float(np.std(data2))
	pValue = stats.norm.sf(abs(z))
	
	return z, pValue


print('P-value of all genes compared to non-disrupted: ', tTest(disruptedTadExpression, nonDisruptedTadExpression))
print('P-value of all genes compared to non-disrupted without 0: ', tTest(filteredDisruptedExpression, filteredNonDisruptedExpression))




