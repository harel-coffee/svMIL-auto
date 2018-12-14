"""
	Set of scripts intended to do ranking of (causal) genes based on their neighborhood and the influence that SVs have in this neighborhood. 

	The idea is to look at all known causal genes. The neighborhood may consist of eQTLs that have an effect on this gene, or TADs directly neighboring the gene.
	If we annotate causal genes with these effects, then we can check for potential effects if anything in the neighborhood is affected by SVs.
	Using some magical statistics, we can then make a prediction on how much more likely the effects on the neighborhood are to be disruptive to the causal genes than expected by random chance.
	Later, this format can be extended to non-causal genes as well. 


	The setup of these scripts will initally be (it will not be too pretty and perfect at first):
	
	- The main script where the input (genes) are parsed and the relevant scripts are called
	- The neighborhoodDefiner, which takes SNVs or SVs as input (or both) and links these to the genes and elements in the neighborhood of that gene that these variants disrupt
	- The geneRanking script, which takes annotated neighborhoods of genes as input, and then computes a score for each layer (i.e. genes, TADs, eQTLs)
	- The above 3 scripts need to be repeated 1000 times (or more) with permutations to compute the scores for genes when the variants are randomly distributed
	
	To run these scripts with permutations, the starting point is runRankingWithPermutations.sh. It does not require any parameters (set the file locations in settings.py), and will run 1 normal scoring run and 1000 permutations on the HPC.
	Then when all permutations are completed, you will need to run computePValuesPerGene.py. This script reads a given output directory containing all gene scores for the normal run and 1000 permutation runs. It will compute
	a p-value for each layer and rank the causal genes by which have significant p-values in the most layers.
	
	To run without permutations, this script main.py can be run as: main.py "runName" N, so for example "main.py ABC N" will run the code once without permutations, and write output to a folder named ABC in the
	RankedGenes subfolder. 
	
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
from geneRankingAlphaBeta import GeneRanking
from inputParser import InputParser
from variantShuffler import VariantShuffler
from channelVisualizer import ChannelVisualizer
import settings


###############
###Main code###


#0. Collect all the relevant parameters here for clarity
uuid = sys.argv[1] #This uuid will normally be provided by the sh script when running in parallel
permutationYN = sys.argv[2] #True or False depending on if we want to permute or not
mode = settings.general['mode'] #Either SV or SNV, then the relevant functions for this data type will be called. For now, a combination of data types is not yet implemented. 
#permutationRound is parameter 5, only used when running on the HPC

import pickle

filehandler = open("GenesAndNeighborhoods.pkl", 'rb')
causalGenes = pickle.load(filehandler)
filehandler.close()


#1. Read and parse the causal genes
# 
# causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
# nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
# 
# uniqueCancerTypes = []
# 
# #Combine the genes for now
# causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
# 
# #The combination of SVs and SNVs will come afterwards, because then we will need to map the names of the cancer types correctly. 
# 
# #2. Read the SVs or SNVs depending on the mode.
# 
# variantData = []
# if mode == "SV":
# 	print "Reading SV data"
# 	svFile = settings.files['svFile']
# 	svData = InputParser().getSVsFromFile(svFile)
# 	
# 
# if mode == "SNV":
# 	print "Reading SNV data"
# 	snvFile = settings.files['snvFile']
# 	snvData = InputParser().getSNVsFromFile(snvFile)
# 	
# #This can be done better with an array of parameters,but this is quick and dirty for now	
# if mode == "SV+SNV":
# 	print "Reading SV data"
# 	svFile = settings.files['svFile']
# 	svData = InputParser().getSVsFromFile(svFile)
# 	print "Reading SNV data"
# 	snvFile = settings.files['snvFile']
# 	snvData = InputParser().getSNVsFromFile(snvFile)
# 	
# 	
# #3. If this is a permutation run, we wish to shuffle these SVs.
#  #Check if this run is a permutation or not. The output file name depends on this
# if permutationYN == "True":
# 	print "Shuffling variants"
# 	variantShuffler = VariantShuffler()
# 	#Shuffle the variants, provide the mode such that the function knows how to permute
# 	if mode == "SV":
# 		svData = variantShuffler.shuffleSVs(svData)
# 	if mode == "SNV":
# 		snvData = variantShuffler.shuffleSNVs(snvData)
# 	if mode == "SV+SNV":
# 		svData = variantShuffler.shuffleSVs(svData)
# 		snvData = variantShuffler.shuffleSNVs(snvData)
# 
# #Number of patients
# #print len(np.unique(svData[:,7]))
# 
# 		
# #2. Get the neighborhood for these genes
# if mode == "SV":
# 	print "Defining the neighborhood for the causal genes and the SVs"
# 	NeighborhoodDefiner(causalGenes, svData, None, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
# if mode == "SNV":
# 	print "Defining the neighborhood for the causal genes and the SNVs"
# 	NeighborhoodDefiner(causalGenes, None, snvData, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
# if mode == "SV+SNV":
# 	print "Defining the neighborhood for the causal genes and the SVs and SNVs"
# 	NeighborhoodDefiner(causalGenes, svData, snvData, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
# 	
# 	
#3. Do simple ranking of the genes and report the causal SVs
print "Ranking the genes for the variants"
#geneRanking = GeneRanking(causalGenes[:,3], mode)
#Skip the ranking for now and instead do exploration

#Save the causal genes up until here and load them for faster development
# import pickle
# 
# filehandler = open("GenesAndNeighborhoods.pkl", 'wb')
# pickle.dump(causalGenes, filehandler)
# filehandler.close()
# exit()

ChannelVisualizer(causalGenes[:,3], mode)

exit()


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
	
	
	perGeneScores = np.empty([len(causalGenes), 5], dtype="object") #store by gene name because it is easiest to match back later
	
	for row in range(0, cancerTypeScores.shape[0]):
		gene = cancerTypeScores[row][0]
		geneName = gene.name
		
		# geneScore = cancerTypeScores[row,1]
		eQTLScore = cancerTypeScores[row,1]
		alpha = cancerTypeScores[row,2]
		beta = cancerTypeScores[row,3]
		eQTLCount = cancerTypeScores[row,4]
		# tadScore = cancerTypeScores[row,3]
		# interactionScore = cancerTypeScores[row,4]
		
		perGeneScores[row][0] = geneName
		perGeneScores[row][1] = eQTLScore
		perGeneScores[row][2] = alpha
		perGeneScores[row][3] = beta
		perGeneScores[row][4] = eQTLCount
		# perGeneScores[row][2] = eQTLScore
		# perGeneScores[row][3] = tadScore
		# perGeneScores[row][4] = interactionScore

	
	cancerTypeFolder = rankedGeneScoreDir + "/" + uuid + "/" + cancerType
	if not os.path.exists(cancerTypeFolder):
		os.makedirs(cancerTypeFolder)

	if permutationYN == "True":
		permutationRound = sys.argv[3]	
		outfileName = cancerTypeFolder + "/permutedSVs_" + permutationRound + "_geneScores.txt"
	else:
		outfileName = cancerTypeFolder + "/realSVs_geneScores.txt"
		
		
	#Write to numpy output file	
	np.savetxt(outfileName, perGeneScores, delimiter='\t', fmt='%s')
