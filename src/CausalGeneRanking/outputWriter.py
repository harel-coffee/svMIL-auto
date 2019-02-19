import numpy as np
import os
import settings

class OutputWriter:
	"""
		The purpose of this class is to write the final gene ranking to an output file. 
	"""
	
	def setPaths(self, rankedGeneScoreDir, uuid):
		"""
			Create the paths and directories that we will write output to.
			
			rankedGeneScoreDir: string with the folder name that should be created directly inside the current folder.
			uuid: string with the folder name that will be created inside the rankedGeneScoreDir, where the gene score files will be written to. 
		"""
		
		if not os.path.exists(rankedGeneScoreDir):
			os.makedirs(rankedGeneScoreDir)
		if not os.path.exists(rankedGeneScoreDir + "/" + uuid):
			os.makedirs(rankedGeneScoreDir + "/" + uuid) #this should be unique, so I now avoid checking if the directory exists. Could later be a thing from the sh file 
		
	
	def writeOutput(self, geneRanking, genes, uuid, permutationYN, permutationRound):
		"""
			Writes the gene ranking to a file. Genes can be ranked by every score, depending on which column is specified here.
			Every row in the file will correspond to a gene. Every column is a different data type.
			
			Current columns format:
			0 - Gene name
			1 - Gene score (number of samples with an SV overlapping that gene)
			2 - eQTL gains (number of samples with a gain of eQTL for that gene)
			3 - eQTL loss (number of samples with a loss of eQTL for that gene)
			4 - enhancer gains (number of samples with a gain of enhancers for that gene)
			5 - enhancer losses (number of asmples with a loss of enahncers for that gene)
			
			geneRanking: (object) GeneRanking object with the scores associated to each gene.
			genes: (numpy array) NP array with all the genes and their information
			uuid: (string) ID of the run, where the results will be written to (see setPaths())
			permutationYN: (string) is this a permutation run, True or False? Determines if the scores will be written to a file with the permutation prefix or real prefix. 
			permutationRound: (int) permutation round number
			
		"""
		
		#First set the paths that 		
		rankedGeneScoreDir = settings.files['rankedGeneScoreDir'] 
		self.setPaths(rankedGeneScoreDir, uuid)
		
		#Obtain a numpy matrix with the scores per gene
		#Format: a file per cancer type.
		#Each row corresponds to a gene. Each gene will have a score for the eQTLs, TADs and Gene itself.
		for cancerType in geneRanking.scores:

			cancerTypeScores = geneRanking.scores[cancerType]

			#Store all the gene scores in here. 
			perGeneScores = np.empty([len(genes), 7], dtype="object") #store by gene name because it is easiest to match back later
			
			for row in range(0, cancerTypeScores.shape[0]):
				gene = cancerTypeScores[row][0]
				geneName = gene.name
				
				geneScore = cancerTypeScores[row,1]
				eQTLGainScore = cancerTypeScores[row,2]
				eQTLLossScore = cancerTypeScores[row,3]
				enhancerGainScore = cancerTypeScores[row,4]
				enhancerLossScore = cancerTypeScores[row,5]
				
				perGeneScores[row][0] = geneName
				perGeneScores[row][1] = geneScore
				perGeneScores[row][2] = eQTLGainScore
				perGeneScores[row][3] = eQTLLossScore
				perGeneScores[row][4] = enhancerGainScore
				perGeneScores[row][5] = enhancerLossScore
				perGeneScores[row][6] = enhancerGainScore + enhancerLossScore #data type to rank by. 

		
			#Also rank the output by highest total score (recurrence)
			perGeneScores = perGeneScores[perGeneScores[:,6].argsort()[::-1]] #Select the column  to rank by
			
			#Create the folder to write the output to specific for the current cancer type
			cancerTypeFolder = rankedGeneScoreDir + "/" + uuid + "/" + cancerType
			if not os.path.exists(cancerTypeFolder):
				os.makedirs(cancerTypeFolder)
		
			#change the output file name depending on if this is a permutation round or not
			if permutationYN == "True" or settings.general['shuffleTads'] == True:	
				outfileName = cancerTypeFolder + "/permutedSVs_" + permutationRound + "_geneScores.txt"
			else:
				outfileName = cancerTypeFolder + "/realSVs_geneScores.txt"
				
			#Write to numpy output file	
			np.savetxt(outfileName, perGeneScores, delimiter='\t', fmt='%s')
