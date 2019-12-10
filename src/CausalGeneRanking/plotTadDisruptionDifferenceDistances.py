"""
	Make an actual plot of the tad disruptions, comparing genes in similar locations to each other. 

"""

import sys
import numpy as np
import settings
from inputParser import InputParser

#Determine the number of bins that we will use for each TAD
bins = 10 #10 on each side.

#First collect the genes that have a z-score (filtered for mutation effects) and get their positions within the TAD
zScores = np.loadtxt('zScores.txt', dtype='object')

splitZScores = []
for zScore in zScores:
	splitScore = zScore[0].split("_")
	splitZScores.append([splitScore[0], splitScore[1], zScore[1]])
		


zScores = np.array(splitZScores, dtype='object')

geneNameConversionFile = sys.argv[1]
allGenes = []
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
		
		if geneName in zScores[:,1]:
			allGenes.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), ensgId, geneName])
		
allGenes = np.array(allGenes, dtype='object')


#then go through the TADs that are disrupted by a non-coding SV. 

#Get all SVs
svDir = settings.files['svDir']
svData = InputParser().getSVsFromFile_hmf(svDir)

#Filter out the coding effect SVs, we want to focus on non-coding SVs. 
excludedSVs = np.loadtxt(settings.files['excludedSVs'], dtype='object')

#svType = 'DEL'

filteredSVs = []
types = []
for sv in svData:

	#if sv[8].svType != svType:
	#	continue
	
	svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName
	if svEntry not in excludedSVs:
		filteredSVs.append(sv)
		
		if sv[8].svType not in types:
			types.append(sv[8].svType)

	

filteredSVs = np.array(filteredSVs, dtype='object')

patients = np.unique(filteredSVs[:,7])

print('types: ', types)

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
	startMatches = (intraSVSubset[:,1] <= tad[1]) * (tad[1] <= intraSVSubset[:,5]) * (tad[2] >= intraSVSubset[:,5])
	#Then also check if there are SVs that start before the TAD end, but end outside the TAD. 
	endMatches = (intraSVSubset[:,5] >= tad[2]) * (tad[1] <= intraSVSubset[:,1]) * (tad[2] >= intraSVSubset[:,1])
	
	#either of these must be true for the TAD to be disrupted by an intra-chromosomal SV.
	allMatches = intraSVSubset[startMatches + endMatches]
	
	for match in allMatches:
		
		if match[7] not in tadDisruptions[tadStr]:
			tadDisruptions[tadStr].append(match[7])
	
	#then repeat for interchromosomal SVs.
	#check if there is any bp in this TAD to count it as disrupted.
	svChr2Subset = filteredSVs[filteredSVs[:,3] == tad[0]]
	#interSVSubset = svChr2Subset[svChr2Subset[:,0] == svChr2Subset[:,3]]
	
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
		
		if match[7] not in tadDisruptions[tadStr]:
			tadDisruptions[tadStr].append(match[7])
	
	for match in allChr2Matches:
		
		if match[7] not in tadDisruptions[tadStr]:
			tadDisruptions[tadStr].append(match[7])

binZScores = dict()
for binInd in range(0, bins):
		
	if binInd not in binZScores:
		binZScores[binInd] = []

for tad in tadDisruptions:
	
	splitTad = tad.split('_')
	
	#Make a mapping for positions to the right bin.
	
	#determine the size and how large each bin should be
	binSize = (float(splitTad[2]) - float(splitTad[1])) / bins
	
	
	currentStart = float(splitTad[1]) #start at the TAD start
	binStarts = [currentStart] #list at which position each bin should start.
	for binInd in range(0, bins):
		
		currentStart += binSize
		binStarts.append(currentStart)

		
	#Go through the genes; find the genes that will be in this bin
	geneChrSubset = allGenes[allGenes[:,0] == splitTad[0]]
	
	for binInd in range(0, len(binStarts)-1):
		
		#get the genes in this bin
		genes = geneChrSubset[(geneChrSubset[:,2] >= binStarts[binInd]) * (geneChrSubset[:,1] <= binStarts[binInd+1])]
		
		#get the z-scores of these genes
		
		for gene in genes:
			geneName = gene[4]
			
			geneZScores = zScores[zScores[:,1] == geneName]
			
			#keep the z-scores separate for each patient
			for patient in range(0, len(geneZScores[:,0])):
				if geneZScores[patient,2] != 'nan':
					binZScores[binInd].append(np.log(float(geneZScores[patient,2])))

#Make a series of boxplots

import matplotlib.pyplot as plt

allData = []
for binInd in range(0, bins):
	allData.append(binZScores[binInd])

plt.boxplot(allData)
#plt.ylim([0, 1900])
plt.show()
	








