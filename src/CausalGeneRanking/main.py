"""
	Set of scripts intended to do ranking of (causal) genes based on their neighborhood and the influence that SVs have in this neighborhood. 

	The idea is to look at all known causal genes. The neighborhood may consist of eQTLs that have an effect on this gene, or TADs directly neighboring the gene.
	If we annotate causal genes with these effects, then we can check for potential effects if anything in the neighborhood is affected by SVs.
	Using some magical statistics, we can then make a prediction on how much more likely the effects on the neighborhood are to be disruptive to the causal genes than expected by random chance.
	Later, this format can be extended to non-causal genes as well. 


	The setup of these scripts will initally be (it will not be too pretty and perfect at first):
	
	- The main script where the input (genes) are parsed and the relevant scripts are called
	- The getNeighborhood class, where things such as neighboring TADs, related eQTLs and overlapping SVs are obtained
	- The RankGenes class, where the genes are ranked by how likely they are causal for the disease depending on the SVs affecting their neighborhood
	- The output will be a matrix with all causal genes in the columns, and the relevant SVs in the rows.
	
	
	Using a gene-based approach will likely be quicker than an SV-based approach, and we can get the relevant SVs through the genes. If an SV is never affecting any of our features defined as interesting, there is no
	need to look at that SV at all. This idea may change along the way.
	

"""
#############
###Imports###

import sys
import numpy as np
import random
import pickle as pkl
import os


from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking
from inputParser import InputParser

###############
###Main code###


#0. Collect all the relevant parameters here for clarity
causalGeneFile = sys.argv[1]
uuid = sys.argv[2] #This uuid will normally be provided by the sh script when running in parallel
permutationYN = sys.argv[3] #True or False depending on if we want to permute or not
mode = sys.argv[4] #Either SV or SNV, then the relevant functions for this data type will be called. For now, a combination of data types is not yet implemented. 


#1. Read and parse the causal genes

causalGenes = InputParser().readCausalGeneFile(causalGeneFile)
uniqueCancerTypes = []

#The combination of SVs and SNVs will come afterwards, because then we will need to map the names of the cancer types correctly. 

#2. Read the SVs or SNVs depending on the mode.

variantData = []
if mode == "SV":		
	svFile = "../../data/TPTNTestSet/TP.txt" #should be a setting, or provide the data as parameter
	variantData = InputParser().getSVsFromFile(svFile)

if mode == "SNV":
	snvFile = "../../data/SNVs/cosmicNCV.txt" #use a simple subset for now because the original file with NC SNVs is very large
	variantData = InputParser().getSNVsFromFile(snvFile)

#3. If this is a permutation run, we wish to shuffle these SVs.
 #Check if this run is a permutation or not. The output file name depends on this
if permutationYN == "True":
	print "Shuffling variants"
	
	#Shuffle the variants, provide the mode such that the function knows how to permute
	
	#svData = shuffleSVs(svData)

#2. Get the neighborhood for these genes
print "Defining the neighborhood for the causal genes and the SVs"
NeighborhoodDefiner(causalGenes, variantData, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)


#3. Do simple ranking of the genes and report the causal SVs
print "Ranking the genes for the SVs"
geneRanking = GeneRanking(causalGenes[:,3], mode)

#Output the ranking scores to a file (should probably also be its own class or at least a function)

rankedGeneScoreDir = "./RankedGenes/" #This should be in the settings later
if not os.path.exists(rankedGeneScoreDir):
    os.makedirs(rankedGeneScoreDir)
if not os.path.exists(rankedGeneScoreDir + "/" + uuid):
	os.makedirs(rankedGeneScoreDir + "/" + uuid) #this should be unique, so I now avoid checking if the directory exists. Could later be a thing from the sh file 


#Obtain a numpy matrix with the scores per gene
#Format: a file per cancer type.
#Each row corresponds to a gene. Each gene will have a score for the eQTLs, TADs and Gene itself.

for cancerType in geneRanking.scores:
	
	
	cancerTypeScores = geneRanking.scores[cancerType]
	
	
	perGeneScores = np.empty([len(causalGenes), 4], dtype="object") #store by gene name because it is easiest to match back later
	
	for row in range(0, cancerTypeScores.shape[0]):
		gene = cancerTypeScores[row][0]
		geneName = gene.name
		
		geneScore = cancerTypeScores[row,1]
		eQTLScore = cancerTypeScores[row,2]
		tadScore = cancerTypeScores[row,3]
		
		perGeneScores[row][0] = geneName
		perGeneScores[row][1] = geneScore
		perGeneScores[row][2] = eQTLScore
		perGeneScores[row][3] = tadScore

	
	cancerTypeFolder = rankedGeneScoreDir + "/" + uuid + "/" + cancerType
	if not os.path.exists(cancerTypeFolder):
		os.makedirs(cancerTypeFolder)

	if permutationYN == "True":
		permutationRound = sys.argv[4]	
		outfileName = cancerTypeFolder + "/permutedSVs_" + permutationRound + "_geneScores.txt"
	else:
		outfileName = cancerTypeFolder + "/realSVs_geneScores.txt"
		
		
	#Write to numpy output file	
	np.savetxt(outfileName, perGeneScores, delimiter='\t', fmt='%s')
