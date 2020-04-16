"""
	For each SV, determine if it starts or ends in a TAD.
	For these TADs, count how many genes are in that TAD.
	Make a plot where we show the number of SVs vs the number of genes that these disrupt
	E.g., 5 genes are disrupted 100 times. 

"""

import sys
import numpy as np
import os


path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVGenePairs/')

from inputParser import InputParser
import settings

outDir = sys.argv[2]
finalOutDir = outDir + '/figureS1/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)


#1. Read the SVs
svDir = settings.files['svDir']
svData = InputParser().getSVsFromFile_hmf(svDir)

#2. read the genes
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#3. Read the TADs

tadFile = settings.files['tadFile']
tadData = InputParser().getTADsFromFile(tadFile)

#4. Map the genes to the TADs
def mapGenesToTads(genes, tadData):
		"""
			For computing effects of disruptions on genes, it is convenient to know which genes are located in which TADs.
			Find out for every TAD which genes are located inside of these, and map them to the TADs. 
		
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
		"""
		
		previousChromosome = 0
		for tad in tadData:
		
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				geneChrSubset = genes[np.where(genes[:,0] == tad[0])] #TADs are intrachromosomal, so looking at 1 chromosome is sufficient
			
			
			#Find all genes that are within the TAD.
			#Because some genes are in two TADs, and the TAD definition is likely not entirely correct, we will count overlap of TADs with any part of the gene even if it is just 1 bp for now.
			#So the start of the gene must be before the end of the TAD, and the end of the gene after the start of the TAD. 
			startMatches = geneChrSubset[:,1] <= tad[2]
			endMatches = geneChrSubset[:,2] >= tad[1]
			
			allMatches = startMatches * endMatches
			matchingGenes = geneChrSubset[allMatches,:]
		
			#Add these genes to the TADs if any.
			if matchingGenes.shape[0] < 1:
				
				#If there are no genes, we can for now add the genes immediately to the left and right of the TAD.
				#With better data, this step can be removed, but since the TADs are so sparse, it is a temporary solution.
				
				#1. Compute the distance from either TAD boundary to all genes
				#We look at all genes to the left, so we can use the end of a gene to the start of a TAD, and the start of a gene to the end of a TAD. 
				#2. Get the gene with the smallest distance to each boundary
				
				startDistances = tad[1] - geneChrSubset[:,2]
				#Exclude all negative elements, these are on the right of the TAD while we want genes on the left of the TAD
				negativeElements = startDistances < 0

				startDistances[negativeElements] = float("inf") #skip the ones that are on the wrong side of the TAD boundary, so set them to inf and they will never be selected. 
				if startDistances[negativeElements].shape[0] == startDistances.shape[0]: #but if everything is inf, we skip this TAD.
					continue
								
				#Repeat but then for the other TAD end. 
				endDistances = geneChrSubset[:,1] - tad[2]
				#Exclude all negative elements
				negativeElements = endDistances < 0

				endDistances[negativeElements] = float("inf") #skip the ones that are on the wrong side of the TAD boundary, so set them to inf and they will never be selected. 
				if endDistances[negativeElements].shape[0] == endDistances.shape[0]: #but if everything is inf, we skip this TAD.
					continue
				
			else:
				tad[3].setGenes(matchingGenes[:,3]) #Add the eQTL objects to the TAD objects. 
	
		return tadData	
tadData = mapGenesToTads(genes, tadData)

#For each SV, determine which TADs it starts and ends in
#Count the number of genes in that TAD
disruptions = [] #store the number of genes that each SV disrupts.
disruptionDict = dict() #store by categories
disruptionDict['1'] = 0
disruptionDict['2-5'] = 0
disruptionDict['6-10'] = 0
disruptionDict['11-20'] = 0
disruptionDict['>20'] = 0
for sv in svData:
	
	if sv[0] == sv[3]: #intrachromosomal SV
		
		geneCount = 0 
		
		tadChrSubset = tadData[tadData[:,0] == sv[0]]
		
		#Get the TAD that the SV starts in
		startMatches = (sv[1] >= tadChrSubset[:,1]) * (sv[1] <= tadChrSubset[:,2])
		
		
		startTAD = tadChrSubset[startMatches]
		
		if startTAD.shape[0] > 0:
			geneCount += len(startTAD[0][3].genes)
		
		#Get the TAD that the SV ends in
		endMatches = (sv[5] >= tadChrSubset[:,1]) * (sv[5] <= tadChrSubset[:,2])
		
		endTAD = tadChrSubset[endMatches]
		if endTAD.shape[0] > 0:
			geneCount += len(endTAD[0][3].genes)
		
		#Count the number of genes in the start & end TAD
		
		if geneCount > 0:
			disruptions.append(geneCount)
			
		if geneCount == 1:
			disruptionDict['1'] += 1
		if geneCount > 1 and geneCount < 6:
			disruptionDict['2-5'] += 1
		if geneCount > 5 and geneCount < 11:
			disruptionDict['6-10'] += 1
		if geneCount > 10 and geneCount < 21:
			disruptionDict['11-20'] += 1
		if geneCount > 20:
			disruptionDict['>20'] += 1

	else: #interchromosomal
		geneCount = 0 
		
		tadChrSubset = tadData[tadData[:,0] == sv[0]]
		
		#Get the TAD that the SV starts in
		startMatches = (sv[1] >= tadChrSubset[:,1]) * (sv[1] <= tadChrSubset[:,2])
		
		
		startTAD = tadChrSubset[startMatches]
		
		if startTAD.shape[0] > 0:
			geneCount += len(startTAD[0][3].genes)
			
		tadChr2Subset = tadData[tadData[:,0] == sv[3]]
		
		#Get the TAD that the SV ends in
		endMatches = (sv[5] >= tadChr2Subset[:,1]) * (sv[5] <= tadChr2Subset[:,2])
		
		endTAD = tadChr2Subset[endMatches]
		if endTAD.shape[0] > 0:
			geneCount += len(endTAD[0][3].genes)
		
		#Count the number of genes in the start & end TAD
		
		if geneCount > 0:
			disruptions.append(geneCount)
		
		if geneCount == 1:
			disruptionDict['1'] += 1
		if geneCount > 1 and geneCount < 6:
			disruptionDict['2-5'] += 1
		if geneCount > 5 and geneCount < 11:
			disruptionDict['6-10'] += 1
		if geneCount > 10 and geneCount < 21:
			disruptionDict['11-20'] += 1
		if geneCount > 20:
			disruptionDict['>20'] += 1

print(disruptions)

#Make a histogram showing how frequently how many genes are disrupted
import matplotlib.pyplot as plt

labels = ['1', '2-5', '6-10', '11-20', '>20']
values = [disruptionDict['1'], disruptionDict['2-5'], disruptionDict['6-10'], disruptionDict['11-20'], disruptionDict['>20']]
		
fig, ax = plt.subplots()
y = np.arange(len(labels))
ax.bar(y, values)
ax.set_xticks(y)
ax.set_xticklabels(labels)
ax.set_xlabel('Number of genes in TADs with SV breakpoints')
ax.set_ylabel('Number of SVs')
	
plt.savefig(finalOutDir + '/figS1.svg')



