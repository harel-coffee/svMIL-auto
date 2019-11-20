from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import settings
import sys
import os
from six.moves import range
import re
import pickle as pkl

class GeneRanking:
	"""
		Class responsible for ranking genes by their causality given the SVs in their neighborhood.
		
		Applying the following rules:
		
		For each cancer type, gather all samples in which SVs have been detected.
		
		- All scores are stored in a SV x gene matrix, to allow SVs to affect multiple genes. 
		- If a gene is directly affected by an SV, the score is 1.
		- Currently, we look sample-based, meaning that we only count the number of SVs affecting that gene once per sample. 
		- For each sample of a cancer type where we compute these scores, we take the sum across the samples.
		- The genes are then ranked by their final score. The gene with the highest score is most likely causal.
	
		TO DO:
		- SNVs can currently not be used in the ranking. 
	
	"""
	
	def __init__(self, genes, svData, mode, runId, permutationRound):
		"""
			genes: (numpy array) array with the genes and their information. chr	start	end	Gene (object)
			
			TO DO:
			- svData and mode are currently not used, need to be added later if we include SNVs. 
			
		"""
		self.scoreGenes(genes, svData, runId, permutationRound)
	
			
	def scoreGenes(self, genes, svData, runId, permutationRound):
		"""
			Score the provided genes. Currently we score by:
				- SVs in genes
				- Gains/losses of eQTLs
				- Gains/losses of enhancers
				
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject
		"""
		
		cosmicGenes = [] 
		with open(settings.files['causalGenesFile'], 'r') as geneFile:
			
			lineCount = 0
			header = []
			for line in geneFile:
				splitLine = line.split("\t")
				#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
				if lineCount < 1:
		
					header = splitLine
					lineCount += 1
					continue
					
				#Obtain the gene name and gene position
				
				geneSymbolInd = header.index('Gene Symbol')
				cosmicGenes.append(splitLine[geneSymbolInd])
		
		sampleMap = dict() #samples and their index in the final scoring matrix. 
		geneMap = dict() #genes and their index in the final scoring matrix
		reverseGeneMap = dict() #also keep a map where we can search by index to later obtain back the gene from the matrix.
		reverseSampleMap = dict()
		scores = dict() #final scores per gene
		#1. Get all unique cancer types and map the gene objects to the right cancer type
		cancerTypes = np.unique(svData[:,6])
		#2. Make the sample map
		samples = np.unique(svData[:,7])
		for sampleInd in range(0, len(samples)):
			sampleMap[samples[sampleInd]] = sampleInd
		for sample in sampleMap:
			index = sampleMap[sample]
			reverseSampleMap[index] = sample
		
		#Make the gene maps
		geneIndex = 0
		for gene in genes:
			if gene not in geneMap:
				geneMap[gene] = geneIndex
				geneIndex += 1
		for gene in geneMap:
			index = geneMap[gene]
			reverseGeneMap[index] = gene
		
		
		print("doing the scoring")
		#Score per cancer type individually
		for cancerType in cancerTypes:
			print("current cancer type: ", cancerType)

			###1. Scoring per SV-gene pair

			#Do the scoring of the genes
			#We make a scoring matrix of patients x genes. Each gene has a score in each patient of if an SV overlaps with that element in the neighborhood of the gene yes/no.
			#To get the total score for a gene, we can sum across all patients.
			eQTLLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, dict(), [], "eQTL")
			enhancerLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "enhancer")
			promoterLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "promoter")
			cpgLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "cpg")
			tfLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "tf")
			hicLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "hic")
			h3k9me3LossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k9me3")
			h3k4me3LossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k4me3")
			h3k27acLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k27ac")
			h3k27me3LossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k27me3")
			h3k4me1LossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k4me1")
			h3k36me3LossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "h3k36me3")
			dnaseILossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "dnaseI")
			ctcfLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "CTCF")
			ctcfEnhancerLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "CTCF+Enhancer")
			ctcfPromoterLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "CTCF+Promoter")
			chromhmmEnhancerLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Enhancer")
			heterochromatinLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Heterochromatin")
			poisedPromoterLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Poised_Promoter")
			chromhmmPromoterLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Promoter")
			repeatLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Repeat")
			repressedLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Repressed")
			transcribedLossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, "Transcribed")
			
			
			eQTLGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "eQTL")
			enhancerGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "enhancer")
			promoterGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "promoter")
			cpgGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "cpg")
			tfGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "tf")
			hicGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "hic")
			h3k9me3GainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k9me3")
			h3k4me3GainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k4me3")
			h3k27acGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k27ac")
			h3k27me3GainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k27me3")
			h3k4me1GainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k4me1")
			h3k36me3GainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "h3k36me3")
			dnaseIGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "dnaseI")
			ctcfGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "CTCF")
			ctcfEnhancerGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "CTCF+Enhancer")
			ctcfPromoterGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "CTCF+Promoter")
			chromhmmEnhancerGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Enhancer")
			heterochromatinGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Heterochromatin")
			poisedPromoterGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Poised_Promoter")
			chromhmmPromoterGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Promoter")
			repeatGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Repeat")
			repressedGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Repressed")
			transcribedGainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, "Transcribed")
			
			#Iterate across the svGeneMap
			#Write the scoring matrix to a file, storing the entire matrix in memory may be too heavy
			pairScores = np.zeros([len(svGeneIndices), 51])
			pairIds = []
			for ind in range(0, len(svGeneIndices)):
				sv = svGeneIndices[ind]
				
				pairIds.append(sv)
				
				
				if sv in eQTLLossScores:
					pairScores[ind,0] = eQTLLossScores[sv]
				if sv in enhancerLossScores:
					pairScores[ind,1] = enhancerLossScores[sv]
				if sv in promoterLossScores:
					pairScores[ind,2] = promoterLossScores[sv]
				if sv in cpgLossScores:
					pairScores[ind,3] = cpgLossScores[sv]
				if sv in tfLossScores:
					pairScores[ind,4] = tfLossScores[sv]
				if sv in hicLossScores:
					pairScores[ind,5] = hicLossScores[sv]
				
				if sv in h3k9me3LossScores:
					pairScores[ind,6] = h3k9me3LossScores[sv]
				if sv in h3k4me3LossScores:
					pairScores[ind,7] = h3k4me3LossScores[sv]
				if sv in h3k27acLossScores:
					pairScores[ind,8] = h3k27acLossScores[sv]
				if sv in h3k27me3LossScores:
					pairScores[ind,9] = h3k27me3LossScores[sv]
				if sv in h3k4me1LossScores:
					pairScores[ind,10] = h3k4me1LossScores[sv]
				if sv in h3k36me3LossScores:
					pairScores[ind,11] = h3k36me3LossScores[sv]
					
				if sv in dnaseILossScores:
					pairScores[ind,12] = dnaseILossScores[sv]
				
				if sv in ctcfLossScores:
					pairScores[ind,13] = ctcfLossScores[sv]
				if sv in ctcfEnhancerLossScores:
					pairScores[ind,14] = ctcfEnhancerLossScores[sv]
				if sv in ctcfPromoterLossScores:
					pairScores[ind,15] = ctcfPromoterLossScores[sv]
				if sv in chromhmmEnhancerLossScores:
					pairScores[ind,16] = chromhmmEnhancerLossScores[sv]
				if sv in heterochromatinLossScores:
					pairScores[ind,17] = heterochromatinLossScores[sv]
				if sv in poisedPromoterLossScores:
					pairScores[ind,18] = poisedPromoterLossScores[sv]
				if sv in chromhmmPromoterLossScores:
					pairScores[ind,19] = chromhmmPromoterLossScores[sv]
				if sv in repeatLossScores:
					pairScores[ind,20] = repeatLossScores[sv]
				if sv in repressedLossScores:
					pairScores[ind,21] = repressedLossScores[sv]
				if sv in transcribedLossScores:
					pairScores[ind,22] = transcribedLossScores[sv]
					
				#Repeat for gains
					
				if sv in eQTLGainScores:
					pairScores[ind,23] = eQTLGainScores[sv]
				if sv in enhancerGainScores:
					pairScores[ind,24] = enhancerGainScores[sv]
				if sv in promoterGainScores:
					pairScores[ind,25] = promoterGainScores[sv]
				if sv in cpgGainScores:
					pairScores[ind,26] = cpgGainScores[sv]
				if sv in tfGainScores:
					pairScores[ind,27] = tfGainScores[sv]
				if sv in hicGainScores:
					pairScores[ind,28] = hicGainScores[sv]
				
				if sv in h3k9me3GainScores:
					pairScores[ind,29] = h3k9me3GainScores[sv]
				if sv in h3k4me3GainScores:
					pairScores[ind,30] = h3k4me3GainScores[sv]
				if sv in h3k27acGainScores:
					pairScores[ind,31] = h3k27acGainScores[sv]
				if sv in h3k27me3GainScores:
					pairScores[ind,32] = h3k27me3GainScores[sv]
				if sv in h3k4me1GainScores:
					pairScores[ind,33] = h3k4me1GainScores[sv]
				if sv in h3k36me3GainScores:
					pairScores[ind,34] = h3k36me3GainScores[sv]
					
				if sv in dnaseIGainScores:
					pairScores[ind,35] = dnaseIGainScores[sv]	
				
				if sv in ctcfGainScores:
					pairScores[ind,36] = ctcfGainScores[sv]
				if sv in ctcfEnhancerGainScores:
					pairScores[ind,37] = ctcfEnhancerGainScores[sv]
				if sv in ctcfPromoterGainScores:
					pairScores[ind,38] = ctcfPromoterGainScores[sv]
				if sv in chromhmmEnhancerGainScores:
					pairScores[ind,39] = chromhmmEnhancerGainScores[sv]
				if sv in heterochromatinGainScores:
					pairScores[ind,40] = heterochromatinGainScores[sv]
				if sv in poisedPromoterGainScores:
					pairScores[ind,41] = poisedPromoterGainScores[sv]
				if sv in chromhmmPromoterGainScores:
					pairScores[ind,42] = chromhmmPromoterGainScores[sv]
				if sv in repeatGainScores:
					pairScores[ind,43] = repeatGainScores[sv]
				if sv in repressedGainScores:
					pairScores[ind,44] = repressedGainScores[sv]
				if sv in transcribedGainScores:
					pairScores[ind,45] = transcribedGainScores[sv]
				
				
				#pairScores[ind,26] = np.sum(pairScores[ind,0:25])
				translocation = 0
				deletion = 0
				duplication = 0
				inversion = 0
				#Add some features for SV type
				splitSV = sv.split("_")
				svType = "_".join(splitSV[8:])
				if re.search("trans", svType, re.IGNORECASE):
					translocation = 1
				if re.search("del", svType, re.IGNORECASE):
					deletion = 1
				if re.search("tandem.dup", svType, re.IGNORECASE):
					duplication = 1
				if re.search("invers", svType, re.IGNORECASE):
					inversion = 1
				
				pairScores[ind,46] = deletion
				pairScores[ind,47] = duplication
				pairScores[ind,48] = inversion
				pairScores[ind,49] = translocation
				
				#check cosmic
				if splitSV[0] in cosmicGenes:
					pairScores[ind,50] = 1
				else:
					pairScores[ind,50] = 0
				
				
			print(pairScores)	
			pairIds = np.array(pairIds)
			
			if not os.path.exists(settings.files['rankedGeneScoreDir']):
				os.makedirs(settings.files['rankedGeneScoreDir'])
			if not os.path.exists(settings.files['rankedGeneScoreDir'] + "/" + runId):
				os.makedirs(settings.files['rankedGeneScoreDir'] + "/" + runId)
			if not os.path.exists(settings.files['rankedGeneScoreDir'] + "/" + runId + '/' + cancerType):
				os.makedirs(settings.files['rankedGeneScoreDir'] + "/" + runId + '/' + cancerType)
			
			pairScoresWithPairIds = np.empty([len(svGeneIndices), 52], dtype="object")
			pairScoresWithPairIds[:,0] = pairIds
			pairScoresWithPairIds[:,1:52] = pairScores
			
			#pairScoresWithPairIds = pairScoresWithPairIds[pairScoresWithPairIds[:,27].argsort()[::-1]] #Select the column  to rank by
			print(pairScoresWithPairIds)
			np.savetxt(settings.files['rankedGeneScoreDir'] + '/' + runId+ '/' + cancerType + "/nonCoding_geneSVPairs.txt_" + str(permutationRound), pairScoresWithPairIds, delimiter='\t', fmt='%s')
			
			#Also output the coding pairs
			codingPairs = []
			for gene in genes:
				for sv in gene.SVs:
					
					codingPairs.append(gene.name + "_" + sv)
			codingPairs = np.array(codingPairs, dtype="object")
			np.savetxt(settings.files['rankedGeneScoreDir'] + '/' + runId + '/' + cancerType + "/coding_geneSVPairs.txt_" + str(permutationRound), codingPairs, delimiter='\t', fmt='%s')		
				
			
			features = ['eQTL', 'enhancer', 'promoter', 'cpg', 'tf', 'hic', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'dnaseI',
						'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed']
			matrices = []
			geneScoringMatrix, geneGeneMatrix = self.scoreBySVsInGenes(genes, sampleMap, geneMap)
			for feature in features:
				scoringMatrix, geneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, feature)
				matrices.append(scoringMatrix)
				scoringMatrix, geneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, feature)
				matrices.append(scoringMatrix)
			
			
			# print("scoring genes:")
			# geneScoringMatrix, geneGeneMatrix = self.scoreBySVsInGenes(genes, sampleMap, geneMap)
			# print("scoring eQTL: ")
			# eQTLGainsScoringMatrix, eQTLGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "eQTL")
			# eQTLLossesScoringMatrix, eQTLLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "eQTL")
			# print("scoring enhancers: ")
			# enhancerGainsScoringMatrix, enhancerGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "enhancer")
			# enhancerLossesScoringMatrix, enhancerLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "enhancer")
			# print("scoring promoters: ")
			# promoterGainsScoringMatrix, promoterGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "promoter")
			# promoterLossesScoringMatrix, promoterLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "promoter")
			# print("scoring cpg: ")
			# cpgGainsScoringMatrix, cpgGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "cpg")
			# cpgLossesScoringMatrix, cpgLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "cpg")
			# print("scoring tfs: ")
			# tfGainsScoringMatrix, tfGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "tf")
			# tfLossesScoringMatrix, tfLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "tf")
			# print("scoring hic: ")
			# hicGainsScoringMatrix, hicGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "hic")
			# hicLossesScoringMatrix, hicLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "hic")
			# 
			# print("scoring histone marks: ")
			# h3k9me3GainsScoringMatrix, h3k9me3GainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k9me3")
			# h3k9me3LossesScoringMatrix, h3k9me3LossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k9me3")
			# h3k4me3GainsScoringMatrix, h3k4me3GainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k4me3")
			# h3k4me3LossesScoringMatrix, h3k4me3LossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k4me3")
			# h3k27acGainsScoringMatrix, h3k27acGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k27ac")
			# h3k27acLossesScoringMatrix, h3k27acLossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k27ac")
			# h3k27me3GainsScoringMatrix, h3k27me3GainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k27me3")
			# h3k27me3LossesScoringMatrix, h3k27me3LossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k27me3")
			# h3k4me1GainsScoringMatrix, h3k4me1GainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k4me1")
			# h3k4me1LossesScoringMatrix, h3k4me1LossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k4me1")
			# h3k36me3GainsScoringMatrix, h3k36me3GainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k36me3")
			# h3k36me3LossesScoringMatrix, h3k36me3LossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k36me3")
			# 
			# print("scoring dnaseI: ")
			# dnaseIGainsScoringMatrix, dnaseIGainsGeneMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "dnaseI")
			# dnaseILossesScoringMatrix, dnaseILossesGeneMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "dnaseI")
			# 
			#Add the chromhmm states as well
			
			#The final scoring matrix selected here determines how the genes will be ranked. For testing, this varies. Eventually, this will be a combined score that we rank by. 
			scoringMatrix = geneScoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.
			
			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort by highest final score and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]
			# 
			# #Now map the indices of the scoring matrix back to the actual genes, and report the scores in the different layers per gene. 
			geneScores = []
			# patientsTotalScoreMatrix = np.zeros([len(geneMap)+1, len(sampleMap)+1], dtype="object") #have a matrix where we store the total scores per gene per patient
			# genePatientPairScores = []
			# for sample in sampleMap:
			# 	sampleInd = sampleMap[sample]
			# 	patientsTotalScoreMatrix[0,sampleInd+1] = sample
			# print("	collecting scores per gene: ")
			
			for geneInd in sortedGenesInd:
				
				gene = reverseGeneMap[geneInd] #Get the gene back from the scoring matrix by index
				
				sampleIndices = []
				
				for matrix in matrices:
					sampleIndices += list(np.where(matrix[:,geneInd] == 1)[0])
					# 
					# #Get a unique list of samples in which the gene is affected (combined across features)
					# sampleIndices += list(np.where(eQTLGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(eQTLLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(enhancerGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(enhancerLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(promoterGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(promoterLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(cpgGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(cpgLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(tfGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(tfLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(hicGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(hicLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k9me3GainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k9me3LossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k4me3GainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k4me3LossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k27acGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k27acLossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k27me3GainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k27me3LossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k4me1GainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k4me1LossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k36me3GainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(h3k36me3LossesScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(dnaseIGainsScoringMatrix[:,geneInd] == 1)[0])
					# sampleIndices += list(np.where(dnaseILossesScoringMatrix[:,geneInd] == 1)[0])
				
				sampleIndices = np.unique(sampleIndices)
				
				samples = []
				for sampleInd in sampleIndices:
					samples.append(reverseSampleMap[sampleInd])
				if len(samples) < 1:
					samples = "None"
				else:
					samples = ",".join(samples)
				
				codingSampleIndices = list(np.where(geneScoringMatrix[:,geneInd] == 1)[0])
				
				codingSamples = []
				for sampleInd in codingSampleIndices:
					codingSamples.append(reverseSampleMap[sampleInd])
				if len(codingSamples) < 1:
					codingSamples = "None"
				else:
					codingSamples = ",".join(codingSamples)
				
				
				# geneScores.append([gene, geneGeneMatrix[geneInd], eQTLGainsGeneMatrix[geneInd], eQTLLossesGeneMatrix[geneInd],
				# 				   enhancerGainsGeneMatrix[geneInd], enhancerLossesGeneMatrix[geneInd], promoterGainsGeneMatrix[geneInd], promoterLossesGeneMatrix[geneInd],
				# 				   cpgGainsGeneMatrix[geneInd], cpgLossesGeneMatrix[geneInd], tfGainsGeneMatrix[geneInd], tfLossesGeneMatrix[geneInd],
				# 				   hicGainsGeneMatrix[geneInd], hicLossesGeneMatrix[geneInd],
				# 				   h3k9me3GainsGeneMatrix[geneInd], h3k9me3LossesGeneMatrix[geneInd], h3k4me3GainsGeneMatrix[geneInd], h3k4me3LossesGeneMatrix[geneInd],
				# 				   h3k27acGainsGeneMatrix[geneInd], h3k27acLossesGeneMatrix[geneInd], h3k27me3GainsGeneMatrix[geneInd], h3k27me3LossesGeneMatrix[geneInd],
				# 				   h3k4me1GainsGeneMatrix[geneInd], h3k4me1LossesGeneMatrix[geneInd], h3k36me3GainsGeneMatrix[geneInd], h3k36me3LossesGeneMatrix[geneInd],
				# 				   dnaseIGainsGeneMatrix[geneInd], dnaseILossesGeneMatrix[geneInd], samples])
				# 
				
				indivScores = [gene, np.sum(geneScoringMatrix[:,geneInd])]
				
				for matrix in matrices:
					indivScores.append(np.sum(matrix[:,geneInd]))
				
				geneScores.append(indivScores)
				
				indivScores.append(samples)
				indivScores.append(codingSamples)
				
				
				
				# geneScores.append([gene, np.sum(geneScoringMatrix[:,geneInd]), np.sum(eQTLGainsScoringMatrix[:,geneInd]), np.sum(eQTLLossesScoringMatrix[:,geneInd]),
				# 				   np.sum(enhancerGainsScoringMatrix[:,geneInd]), np.sum(enhancerLossesScoringMatrix[:,geneInd]), np.sum(promoterGainsScoringMatrix[:,geneInd]), np.sum(promoterLossesScoringMatrix[:,geneInd]),
				# 				   np.sum(cpgGainsScoringMatrix[:,geneInd]), np.sum(cpgLossesScoringMatrix[:,geneInd]), np.sum(tfGainsScoringMatrix[:,geneInd]), np.sum(tfLossesScoringMatrix[:,geneInd]),
				# 				   np.sum(hicGainsScoringMatrix[:,geneInd]), np.sum(hicLossesScoringMatrix[:,geneInd]),
				# 				   np.sum(h3k9me3GainsScoringMatrix[:,geneInd]), np.sum(h3k9me3LossesScoringMatrix[:,geneInd]), np.sum(h3k4me3GainsScoringMatrix[:,geneInd]), np.sum(h3k4me3LossesScoringMatrix[:,geneInd]),
				# 				   np.sum(h3k27acGainsScoringMatrix[:,geneInd]), np.sum(h3k27acLossesScoringMatrix[:,geneInd]), np.sum(h3k27me3GainsScoringMatrix[:,geneInd]), np.sum(h3k27me3LossesScoringMatrix[:,geneInd]),
				# 				   np.sum(h3k4me1GainsScoringMatrix[:,geneInd]), np.sum(h3k4me1LossesScoringMatrix[:,geneInd]), np.sum(h3k36me3GainsScoringMatrix[:,geneInd]), np.sum(h3k36me3LossesScoringMatrix[:,geneInd]),
				# 				   np.sum(dnaseIGainsScoringMatrix[:,geneInd]), np.sum(dnaseILossesScoringMatrix[:,geneInd]), samples, codingSamples])
				# 
				# 
				
				
			print("done")	
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
	
			###3. Make the feature file for MIL for each sv-gene pair
			#Each SV-gene pair is a bag. A bag can contain a variable set of isntances, which represent the gained/lost elements
			#The feature vector was pre-defined in the gene class for each instance.
			
			bags = dict()
			for geneInd in sortedGenesInd:
				gene = reverseGeneMap[geneInd] #Get the gene back from the scoring matrix by index
				
				
				for sv in gene.alteredElements:
					
					instances = []
					for element in gene.alteredElements[sv]:
						instances.append(gene.alteredElements[sv][element])
					
					if len(instances) > 0:
						bags[gene.name + "_" + sv] = instances
				
			#print(bags)
			print(len(bags))
		
			#output the bags to a file
			with open(settings.files['rankedGeneScoreDir'] + '/' + runId+ '/' + cancerType + '/bags.pkl', 'wb') as handle:
				pkl.dump(bags, handle, protocol=pkl.HIGHEST_PROTOCOL)

	
		self.scores = scores #Currently, we make it part of the object, but in principle this could be returned and that would be a bit nicer.

				
	def scoreBySVsInGenes(self, genes, sampleMap, geneMap):
		"""
			Determine for each gene how many samples have at least 1 SV overlapping with the gene. These were already mapped to the gene in the neighborhoodDefiner.
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			sampleMap: (dictionary) each sample is a key, and the value is the index of where the gene score is in the final scoring matrix.
			geneMap: (dictionary) each gene is a key, and the value is the index of where the gene score is in the final scoring matrix. 
			
			return
			scoringMatrix: (numpy array) matrix of samples x genes (samples in rows, genes in columns) with a value of 1 indicating that the sample has an SV overlapping that gene, or 0 if there are none. 
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		geneMatrix = np.zeros(len(geneMap))
		for geneInd in range(0, len(genes)): #first loop over the genes, then get their SVs
			gene = genes[geneInd]
			

			matrixGeneInd = geneMap[gene] #index of the gene in the scoring matrix
			
			#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
			
			#Perform additional check to see if the gene has SVs at all, otherwise we can skip it
			geneSVs = gene.SVs
			
			if geneSVs is None or len(geneSVs) < 1:
				continue
			
			for sv in geneSVs:

				#check the sample of this sv and get the right position of the samples in the scoring matrix
				splitSV = sv.split("_")
				sampleName = splitSV[len(splitSV)-1]
				sampleInd = sampleMap[sampleName]

				scoringMatrix[sampleInd][matrixGeneInd] = 1 #set the score to 1 to ensure that per sample we count only 1 SV.
				geneMatrix[geneInd] += 1
	
		return scoringMatrix, geneMatrix
	
	def scoreByElementGains(self, genes, sampleMap, geneMap, elementType):
		"""
			Determine for each gene how many elements are gained. 
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			sampleMap: (dictionary) each sample is a key, and the value is the index of where the gene score is in the final scoring matrix.
			geneMap: (dictionary) each gene is a key, and the value is the index of where the gene score is in the final scoring matrix.
			elementType: (string) type of the element that we should score the gains of. 
			
			return
			scoringMatrix: (numpy array) matrix of samples x genes (samples in rows, genes in columns) with a value of 1 indicating that the sample has an SV causing a gain of at least 1 element, or 0 if there are none. 
		"""

		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)]) #store the final scores
		geneMatrix = np.zeros(len(geneMap))
		
		if settings.general['gains'] == False:
			return scoringMatrix, geneMatrix
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]
			geneSVSamples = []
			for sv in gene.SVs:
				splitSV = sv.split("_")
				geneSVSamples.append(splitSV[len(splitSV)-1])

			matrixGeneInd = geneMap[gene]
			
			
			if len(gene.gainedElements) > 0: #if the gene has gained elements
				for sample in gene.gainedElements:
					
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					if settings.general['nonCoding'] == True and settings.general['coding'] == False:
						
						
							
						if sample in geneSVSamples:
							continue
					
					#check if we need to ignore a pair if it has an SNV in that patient.  
					if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 

						if sample in gene.SNVs:
							continue
					
					
					gain = False #Make sure that we only focus on gained elements of the provided type. 
					for element in gene.gainedElements[sample]:
						if element == elementType:
							gain = True
					if gain == True: 
						sampleInd = sampleMap[sample]

						scoringMatrix[sampleInd][matrixGeneInd] = 1
						geneMatrix[geneInd] += 1
				
		
		return scoringMatrix, geneMatrix
	
	def scoreByElementLosses(self, genes, sampleMap, geneMap, elementType):
		"""
			Determine for each gene how many elements are lost. 
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			sampleMap: (dictionary) each sample is a key, and the value is the index of where the gene score is in the final scoring matrix.
			geneMap: (dictionary) each gene is a key, and the value is the index of where the gene score is in the final scoring matrix.
			elementType: (string) type of the element that we should score the losses of. 
			
			return
			scoringMatrix: (numpy array) matrix of samples x genes (samples in rows, genes in columns) with a value of 1 indicating that the sample has an SV causing a loss of at least 1 element, or 0 if there are none. 
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		geneMatrix = np.zeros(len(geneMap))
		
		if settings.general['losses'] == False:
			return scoringMatrix, geneMatrix
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]
			geneSVSamples = []
			for sv in gene.SVs:
				splitSV = sv.split("_")
				geneSVSamples.append(splitSV[len(splitSV)-1])

			matrixGeneInd = geneMap[gene]
			
			
			if len(gene.lostElements) > 0:
				for sample in gene.lostElements:
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					if settings.general['nonCoding'] == True and settings.general['coding'] == False:
						
					
						if sample in geneSVSamples:
							continue
					
					#check if we need to ignore a pair if it has an SNV in that patient.  
					if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 

						if sample in gene.SNVs:
							continue
					
					loss = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.lostElements[sample]:
						if element == elementType:
							loss = True
					if loss == True:
						sampleInd = sampleMap[sample]
						
						scoringMatrix[sampleInd][matrixGeneInd] = 1
						geneMatrix[geneInd] += 1
		
		return scoringMatrix, geneMatrix

	def scoreByElementLossesSVs(self, genes, svGeneMap, svGeneIndices, elementType):
		"""
			Determine for each gene how many elements are lost. 
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			sampleMap: (dictionary) each sample is a key, and the value is the index of where the gene score is in the final scoring matrix.
			geneMap: (dictionary) each gene is a key, and the value is the index of where the gene score is in the final scoring matrix.
			elementType: (string) type of the element that we should score the losses of. 
			
			return
			scoringMatrix: (numpy array) matrix of samples x genes (samples in rows, genes in columns) with a value of 1 indicating that the sample has an SV causing a loss of at least 1 element, or 0 if there are none. 
		"""
		
		if settings.general['losses'] == False:
			return dict(), svGeneMap, svGeneIndices	
		
		#Idea: have sv-gene pairs on one axis, and a 1 or 0 for this particular feature in the first column
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			geneSVSamples = []
			for sv in gene.SVs:
				splitSV = sv.split("_")
				geneSVSamples.append(splitSV[len(splitSV)-1])
			
			
			
			if len(gene.lostElementsSVs) > 0:
				
				for sv in gene.lostElementsSVs:
					pairId = gene.name + "_" + sv
					
					geneSVSamples = []
					for geneSV in gene.SVs:
						splitSV = geneSV.split("_")
						geneSVSamples.append(splitSV[len(splitSV)-1])
					
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					if settings.general['nonCoding'] == True and settings.general['coding'] == False:
						
						splitSV = sv.split("_")

						if splitSV[len(splitSV)-1] in geneSVSamples:
							continue
					
					#check if we need to ignore a pair if it has an SNV in that patient.  
					if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 
						
						splitSV = sv.split("_")
						patientID = splitSV[6]
						if patientID in gene.SNVs:
							continue
					
					if pairId not in svGeneMap: #Check if we already used this index for a different feature
						svInd = len(svGeneIndices)
						svGeneMap[pairId] = svInd
						svGeneIndices.append(pairId)
					else:
						svInd = svGeneMap[pairId]
				
					loss = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.lostElementsSVs[sv]:
						if element == elementType:
							loss = True
					if loss == True:

						pairScores[pairId] = 1 #assume that each SV can disrupt a gene only once
						#pairScores[pairId] = gene.lostElementsSVs[sv][element]
						
		return pairScores, svGeneMap, svGeneIndices		
		
	def scoreByElementGainsSVs(self, genes, svGeneMap, svGeneIndices, elementType):
		"""
			Determine for each gene how many elements are lost. 
			
			genes:  (numpy array) array with the genes and their information. chr, start, end, geneObject
			sampleMap: (dictionary) each sample is a key, and the value is the index of where the gene score is in the final scoring matrix.
			geneMap: (dictionary) each gene is a key, and the value is the index of where the gene score is in the final scoring matrix.
			elementType: (string) type of the element that we should score the losses of. 
			
			return
			scoringMatrix: (numpy array) matrix of samples x genes (samples in rows, genes in columns) with a value of 1 indicating that the sample has an SV causing a loss of at least 1 element, or 0 if there are none. 
		"""
		
		if settings.general['gains'] == False:
			return dict(), svGeneMap, svGeneIndices	
		
		#Idea: have sv-gene pairs on one axis, and a 1 or 0 for this particular feature in the first column
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			
			
			if len(gene.gainedElementsSVs) > 0:
				
				for sv in gene.gainedElementsSVs:
					pairId = gene.name + "_" + sv
					
					geneSVSamples = []
					for geneSV in gene.SVs:
						splitSV = geneSV.split("_")
						geneSVSamples.append(splitSV[len(splitSV)-1])
					
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					if settings.general['nonCoding'] == True and settings.general['coding'] == False:
						
						splitSV = sv.split("_")

						if splitSV[len(splitSV)-1] in geneSVSamples:
							continue
					
					#check if we need to ignore a pair if it has an SNV in that patient.  
					if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 
						
						splitSV = sv.split("_")
						patientID = splitSV[6]
						if patientID in gene.SNVs:
							continue
						
					if pairId not in svGeneMap: #Check if we already used this index for a different feature
						svInd = len(svGeneIndices)
						svGeneMap[pairId] = svInd
						svGeneIndices.append(pairId)
					else:
						svInd = svGeneMap[pairId]
				
					loss = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.gainedElementsSVs[sv]:
						if element == elementType:
							loss = True
					if loss == True:
						pairScores[pairId] = 1 #assume that each SV can disrupt a gene only once
						#pairScores[pairId] = gene.gainedElementsSVs[sv][element]
						
		
		
		return pairScores, svGeneMap, svGeneIndices	

