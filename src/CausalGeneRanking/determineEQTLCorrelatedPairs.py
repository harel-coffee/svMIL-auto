"""
	Script to determine if there are any pairs of eQTLs that are frequently affected together by SVs. 
	
"""

import sys
import numpy as np
import random
import pickle as pkl
import os

from inputParser import InputParser
import settings

from sv import SV
from gene import Gene
from eQTL import EQTL
from variantShuffler import VariantShuffler


shuffle = sys.argv[1]
permutationInd = sys.argv[2]
uuid = sys.argv[3]

#Read the SV input file
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes. 

#Combine the genes for now
allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#Then read all the eQTLs
def mapEQTLsToGenes(eQTL, geneDict, geneSymbol):
		"""
			Map the right gene object to the eQTL object. 
		"""
		
		geneDict[geneSymbol].addEQTL(eQTL)
		eQTL.addGene(geneDict[geneSymbol])
		
def getEQTLsFromFile(eQTLFile, genes):
		#Filter the eQTLs that do not have a match
		geneDict = dict()
		
		for gene in genes:
			if gene not in geneDict:
				geneDict[gene.name] = gene		
		
		eQTLs = []
		eQTLGenes = dict()
		with open(eQTLFile, 'rb') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				
				if splitLine[3] not in geneDict:
					continue
				
				#The mapping information is in the file, so we can already do it here
				eQTLGenes[splitLine[3]] = 0
						
				#Add the chr notation for uniformity. 		
				eQTLs.append(["chr" + splitLine[0], int(splitLine[1]), int(splitLine[2]), splitLine[3]]) #Keep the eQTL information raw as well for quick overlapping. 
		
		return np.array(eQTLs, dtype='object'), eQTLGenes

eQTLFile = "../../data/eQTLs/breast_eQTLs.txt" #These files will need to go to the settings!
[eQTLData, eQTLGenes] = getEQTLsFromFile(eQTLFile, allGenes[:,3])

#Determine all pairs of genes (that have eQTLs)
genePairs = dict()
for gene in eQTLGenes:
	for gene2 in eQTLGenes:
		genePairs[gene + "_" + gene2] = 0


#Count how often these pairs are affected together by an SV
#For each SV, check which eQTLs it overlaps
#Make pairs of each of the eQTLs
#Increase their counter

#Use parameter to determine if SVs should be shuffled or not
svFile = settings.files['svFile']
svData = InputParser().getSVsFromFile(svFile)


if shuffle == "True":
	print "shuffling SVs"
	variantShuffler = VariantShuffler()
	svData = variantShuffler.shuffleSVs(svData)

for sv in svData:
	
	#Match all eQTLs on either chromosome 1 or chromosome 2
	#If chr1, use s1 and e1
	#If chr2, use s2 and e2.
	
	chr1EQTLs = eQTLData[eQTLData[:,0] == sv[0],:]
	chr2EQTLs = eQTLData[eQTLData[:,0] == sv[3],:]
	
	
	#Find the eQTLs where the position is larger than the sv start, and smaller than the SV end.
	
	#Start and end is the same for eQTLs
	chr1StartMatches = chr1EQTLs[:,1] >= sv[1]
	chr1EndMatches = chr1EQTLs[:,1] <= sv[2]
	
	chr1MatchesInd = chr1StartMatches * chr1EndMatches
	chr1Matches = chr1EQTLs[chr1MatchesInd]
	
	chr2StartMatches = chr2EQTLs[:,1] >= sv[1]
	chr2EndMatches = chr2EQTLs[:,1] <= sv[2]
	
	chr2MatchesInd = chr2StartMatches * chr2EndMatches
	chr2Matches = chr2EQTLs[chr2MatchesInd]
	
	allMatchingEQTLs = np.concatenate((chr1Matches, chr2Matches), axis=0)
	
	uniqueGenes = dict()
	
	for matchingEQTL in allMatchingEQTLs:
		uniqueGenes[matchingEQTL[3]] = 0
	
	#Determine the pairs.
	for uniqueGeneInd in range(0, len(uniqueGenes)):
		uniqueGene = uniqueGenes.keys()[uniqueGeneInd]
		for uniqueGene2Ind in range(uniqueGeneInd+1, len(uniqueGenes)):
			uniqueGene2 = uniqueGenes.keys()[uniqueGene2Ind]
			if uniqueGene is not uniqueGene2:
				genePairs[uniqueGene + "_" + uniqueGene2] += 1
	
	


#Rank the gene pairs
genePairValues = np.array(genePairs.values())
genePairKeys = genePairs.keys()

#sort the gene pair scores
sortedGenesInd = np.argsort(genePairValues)[::-1]

#Output the gene pair scores to a file.


fileType = "realSVs_counts.txt"
if shuffle == "True":
	print "writing to the shuffled out file"
	fileType = "shuffledSVs_counts_" + permutationInd + ".txt"
outFile = "./RankedGenes/" + uuid + "/" + fileType


with open(outFile, 'w') as outF:
	
	for geneInd in range(0, len(sortedGenesInd)):
		sortedGeneInd = sortedGenesInd[geneInd]
		line = genePairKeys[sortedGeneInd] + "\t" + str(genePairValues[sortedGeneInd]) + "\n"
		outF.write(line)
	


# 
# for geneInd in range(0, 100): #show only the top 100
# 	sortedGeneInd = sortedGenesInd[geneInd]
# 
# 	print "gene: ", genePairKeys[sortedGeneInd], " score: ", genePairValues[sortedGeneInd]
# 	
# 
# 
# #Read the actual gene scores.
# #Are there any correlations between the scores of the genes that are frequently found together as a pair?
# 
# geneScoreFile = "RankedGenes/0/breast/realSVs_geneScores.txt"
# 
# geneScores = dict()
# with open(geneScoreFile, 'r') as inF:
# 	
# 	for line in inF:
# 		splitLine = line.split("\t")
# 		geneName = splitLine[0]
# 		geneScore = float(splitLine[1])
# 		geneScores[geneName] = geneScore
# 		
# 		
# #Link the scores to the gene pairs
# pair1Scores = [] #first partner of the pair
# pair2Scores = [] #2nd partner of the pair
# 
# for geneInd in range(0,10):
# 	sortedGeneInd = sortedGenesInd[geneInd]
# 	
# 	splitPair = genePairKeys[sortedGeneInd].split("_")
# 	pair1 = splitPair[0]
# 	pair2 = splitPair[1]
# 	
# 	pair1Scores.append(geneScores[pair1])
# 	pair2Scores.append(geneScores[pair2])
# 	
# 	print "gene 1: ", pair1, " score: ", geneScores[pair1], " gene 2: ", pair2, " score: ", geneScores[pair2]
# 
# pair1Scores = np.array(pair1Scores)
# pair2Scores = np.array(pair2Scores)
# 
# print "correlation in highest ranking group: ", np.corrcoef(pair1Scores, pair2Scores)




	




