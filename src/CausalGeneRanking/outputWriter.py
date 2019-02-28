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
			perGeneScores = np.empty([len(genes), 31], dtype="object") #store by gene name because it is easiest to match back later
			
			for row in range(0, cancerTypeScores.shape[0]):
				gene = cancerTypeScores[row][0]
				geneName = gene.name
				geneChr = gene.chromosome
				geneStart = gene.start
				
				geneScore = cancerTypeScores[row,1]
				eQTLGainScore = cancerTypeScores[row,2]
				eQTLLossScore = cancerTypeScores[row,3]
				enhancerGainScore = cancerTypeScores[row,4]
				enhancerLossScore = cancerTypeScores[row,5]
				promoterGainScore = cancerTypeScores[row,6]
				promoterLossScore = cancerTypeScores[row,7]
				cpgGainScore = cancerTypeScores[row,8]
				cpgLossScore = cancerTypeScores[row,9]
				tfGainScore = cancerTypeScores[row,10]
				tfLossScore = cancerTypeScores[row,11]
				hicGainScore = cancerTypeScores[row,12]
				hicLossScore = cancerTypeScores[row,13]
				h3k9me3GainScore = cancerTypeScores[row,14]
				h3k9me3LossScore = cancerTypeScores[row,15]
				h3k4me3GainScore = cancerTypeScores[row,16]
				h3k4me3LossScore = cancerTypeScores[row,17]
				h3k27acGainScore = cancerTypeScores[row,18]
				h3k27acLossScore = cancerTypeScores[row,19]
				h3k27me3GainScore = cancerTypeScores[row,20]
				h3k27me3LossScore = cancerTypeScores[row,21]
				h3k4me1GainScore = cancerTypeScores[row,22]
				h3k4me1LossScore = cancerTypeScores[row,23]
				h3k36me3GainScore = cancerTypeScores[row,24]
				h3k36me3LossScore = cancerTypeScores[row,25]
				dnaseIGainScore = cancerTypeScores[row,26]
				dnaseILossScore = cancerTypeScores[row,27]
				
				
				perGeneScores[row][0] = geneName
				perGeneScores[row][1] = geneChr
				perGeneScores[row][2] = geneStart
				
				perGeneScores[row][3] = geneScore
				perGeneScores[row][4] = eQTLGainScore
				perGeneScores[row][5] = eQTLLossScore
				perGeneScores[row][6] = enhancerGainScore
				perGeneScores[row][7] = enhancerLossScore
				perGeneScores[row][8] = promoterGainScore
				perGeneScores[row][9] = promoterLossScore
				perGeneScores[row][10] = cpgGainScore
				perGeneScores[row][11] = cpgLossScore
				perGeneScores[row][12] = tfGainScore
				perGeneScores[row][13] = tfLossScore
				perGeneScores[row][14] = hicGainScore
				perGeneScores[row][15] = hicLossScore
				perGeneScores[row][16] = h3k9me3GainScore
				perGeneScores[row][17] = h3k9me3LossScore
				perGeneScores[row][18] = h3k4me3GainScore
				perGeneScores[row][19] = h3k4me3LossScore
				perGeneScores[row][20] = h3k27acGainScore
				perGeneScores[row][21] = h3k27acLossScore
				perGeneScores[row][22] = h3k27me3GainScore
				perGeneScores[row][23] = h3k27me3LossScore
				perGeneScores[row][24] = h3k4me1GainScore
				perGeneScores[row][25] = h3k4me1LossScore
				perGeneScores[row][26] = h3k36me3GainScore
				perGeneScores[row][27] = h3k36me3LossScore
				perGeneScores[row][28] = dnaseIGainScore
				perGeneScores[row][29] = dnaseILossScore
				
				
				perGeneScores[row][30] = np.sum(perGeneScores[row][4:29])

		
			#Also rank the output by highest total score (recurrence)
			perGeneScores = perGeneScores[perGeneScores[:,30].argsort()[::-1]] #Select the column  to rank by
			
			#Create the folder to write the output to specific for the current cancer type
			cancerTypeFolder = rankedGeneScoreDir + "/" + uuid + "/" + cancerType
			if not os.path.exists(cancerTypeFolder):
				os.makedirs(cancerTypeFolder)
		
			#change the output file name depending on if this is a permutation round or not
			if permutationYN == "True" or settings.general['shuffleTads'] == True:	
				outfileName = cancerTypeFolder + "/permutedSVs_" + permutationRound + "_geneScores.txt"
			else:
				outfileName = cancerTypeFolder + "/realSVs_geneScores_chr.txt"
			
			header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\ttotal"
				
			#Write to numpy output file	
			np.savetxt(outfileName, perGeneScores, delimiter='\t', fmt='%s', header=header)
