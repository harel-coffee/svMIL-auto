"""
	Set of scripts intended to do ranking of (causal) genes based on their neighborhood and the influence that SVs have in this neighborhood. 

	The idea is to look at all genes. The neighborhood may consist of eQTLs that have an effect on this gene, or TADs directly neighboring the gene.
	If we annotate causal genes with these effects, then we can check for potential effects on the gene (expression) if anything in the neighborhood is affected by SVs.
	Using some magical statistics, we can then make a prediction on how much more likely the effects on the neighborhood are to be disruptive to the genes than expected by random chance.

	The setup of these scripts will initally be (it will not be too pretty and perfect at first):
	
	- The main script where the input (genes) are parsed and the relevant scripts are called
	- The neighborhoodDefiner, which takes SNVs or SVs as input (or both) and links these to the genes and elements in the neighborhood of that gene that these variants disrupt
	- The geneRanking script, which takes annotated neighborhoods of genes as input, and then computes a score for each layer (i.e. genes, TADs, eQTLs)
	- The above 3 scripts need to be repeated 100 times (or more) with permutations to compute the scores for genes when the variants are randomly distributed
	
	To run these scripts with permutations, the starting point is runRankingWithPermutations.sh. It does not require any parameters (set the file locations in settings.py), and will run 1 normal scoring run and 100 permutations on the HPC.
	Then when all permutations are completed, you will need to run computePValuesPerGene.py. This script reads a given output directory containing all gene scores for the normal run and 100 permutation runs. It will compute
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
import re
import time

from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking
from inputParser import InputParser
from genomicShuffler import GenomicShuffler
from channelVisualizer import ChannelVisualizer
from outputWriter import OutputWriter
# from genome import Genome
import settings


###############
###Main code###

startTime = time.time() #Keep run time of the program

#0. Collect all the relevant parameters here for clarity
uuid = sys.argv[1] #This uuid will normally be provided by the sh script when running in parallel. It is the folder name that will be created in RankedGenes, and where the output will be written to. 
permutationYN = sys.argv[2] #True or False depending on if we want to permute or not
mode = settings.general['mode'] #Either SV or SNV, then the relevant functions for this data type will be called.
#permutationRound is parameter 5, only used when running on the HPC

#1. Read and parse the causal genes and the nonCausal genes. For now, keep them separate to test on causal/non-causal genes separately
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#provide list of SVs that should be excluded from the model. 
excludedSVs = np.loadtxt(settings.files['excludedSVs'], dtype='object')

#2. Read the SVs or SNVs depending on the mode.
variantData = []
if mode == "SV":
	print "Reading SV data"
	svFile = settings.files['svFile']
	svData = InputParser().getSVsFromFile(svFile, "all", excludedSVs)
	
	
#3. If this is a permutation run, we wish to shuffle these SVs or SNVs.
if permutationYN == "True":
	print "Shuffling variants"
	genomicShuffler = GenomicShuffler()
	#Shuffle the variants, provide the mode such that the function knows how to permute
	if mode == "SV":
		svData = genomicShuffler.shuffleSVs(svData)
	if mode == "SNV":
		snvData = genomicShuffler.shuffleSNVs(snvData)
	if mode == "SV+SNV":
		svData = genomicShuffler.shuffleSVs(svData)
		snvData = genomicShuffler.shuffleSNVs(snvData)

permutationRound = ""
if permutationYN == "True":
	permutationRound = sys.argv[3]

#2. Get the neighborhood for these genes based on the SVs or SNVs
if mode == "SV":
	print "Defining the neighborhood for the causal genes and the SVs"
	NeighborhoodDefiner(causalGenes, svData, None, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
if mode == "SNV":
	print "Defining the neighborhood for the causal genes and the SNVs"
	NeighborhoodDefiner(causalGenes, None, snvData, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
if mode == "SV+SNV":
	print "Defining the neighborhood for the causal genes and the SVs and SNVs"
	NeighborhoodDefiner(causalGenes, svData, snvData, mode) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
	

#3. Do ranking of the genes and report the causal SVs
print "Ranking the genes for the variants"
geneRanking = GeneRanking(causalGenes[:,3], svData, mode, sys.argv[1], permutationRound)

#Save the causal genes up until here and load them for faster development of the deep learning part
# import pickle
# 
# filehandler = open("GenesAndNeighborhoods.pkl", 'wb')
# pickle.dump(causalGenes, filehandler)
# filehandler.close()
#
# #For using deep learning, call the channel visualizer that makes channels for the disruptions and visualizes that. 
# ChannelVisualizer(causalGenes[:,3], mode, genome)
# 
# exit()


#Output the ranking scores to a file. 


OutputWriter().writeOutput(geneRanking, causalGenes, uuid, permutationYN, permutationRound)

endTime = time.time()
print "The program took ", endTime-startTime, " seconds to complete"

	