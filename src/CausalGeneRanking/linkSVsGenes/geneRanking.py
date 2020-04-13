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
		for cancerType in cancerTypes:
			features = ['eQTL', 'enhancer', 'promoter', 'cpg', 'tf', 'hic', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'dnaseI', 'rnaPol',
							'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'superEnhancer', 'ctcf']
			svGeneMap = dict()
			svGeneIndices = []
			allLossScores = []
			allGainScores = []
			for feature in features:
				
				lossScores, svGeneMap, svGeneIndices = self.scoreByElementLossesSVs(genes, svGeneMap, svGeneIndices, feature)
				allLossScores.append(lossScores)
				
				gainScores, svGeneMap, svGeneIndices = self.scoreByElementGainsSVs(genes, svGeneMap, svGeneIndices, feature)
				allGainScores.append(gainScores)
			
			#get the strength features
			strengthFeatures = ['enhancer', 'ctcf', 'rnaPol', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
			
			allStrengthLossScores = []
			allStrengthGainScores = []
			for feature in strengthFeatures:
				lossScores = self.scoreByElementLossesStrengthsSVs(genes, feature)
				allStrengthLossScores.append(lossScores)
				
				gainScores = self.scoreByElementGainsStrengthsSVs(genes, feature)
				allStrengthGainScores.append(gainScores)
				
				
			delCount = 0
			dupCount = 0
			invCount = 0
			itxCount = 0
			pairScores = np.zeros([len(svGeneIndices), 79])
			pairIds = []
			for ind in range(0, len(svGeneIndices)):
				sv = svGeneIndices[ind]
				
				pairIds.append(sv)
				#print('gain/losses')
				for featureInd in range(0, len(features)):
					if sv in allLossScores[featureInd]:
						pairScores[ind,featureInd] = allLossScores[featureInd][sv]
					
					if sv in allGainScores[featureInd]:
						pairScores[ind,featureInd+len(features)] = allGainScores[featureInd][sv]
					
						
					#print(featureInd, features[featureInd], 'loss')
					#print(featureInd+len(features), features[featureInd], 'gain')
	
				#print('strengths')
				for featureInd in range(0, len(strengthFeatures)):
					if sv in allStrengthLossScores[featureInd]:
						pairScores[ind,(featureInd+len(features)*2)] = allStrengthLossScores[featureInd][sv]
					
					if sv in allStrengthGainScores[featureInd]:
						pairScores[ind,(featureInd+len(features)*2)+len(strengthFeatures)] = allStrengthGainScores[featureInd][sv]
					
					#print(featureInd+len(features)*2, strengthFeatures[featureInd], 'loss')
					#print((featureInd+len(features)*2)+len(strengthFeatures), strengthFeatures[featureInd], 'gain')
					
				translocation = 0
				deletion = 0
				duplication = 0
				inversion = 0
				#Add some features for SV type
				splitSV = sv.split("_")
				svType = splitSV[12]
				if re.search("ITX", svType, re.IGNORECASE):
					translocation = 1
					itxCount += 1
				elif re.search("DEL", svType, re.IGNORECASE):
					deletion = 1
					delCount += 1
				elif re.search("DUP", svType, re.IGNORECASE):
					duplication = 1
					dupCount += 1
				elif re.search("INV", svType, re.IGNORECASE):
					inversion = 1
					invCount += 1
				
				pairScores[ind,70] = deletion
				pairScores[ind,71] = duplication
				pairScores[ind,72] = inversion
				pairScores[ind,73] = translocation

			
			print('del pairs: ', delCount)
			print('inv pairs: ', invCount)
			print('dup pairs: ', dupCount)
			print('itx pairs: ', itxCount)
			
			pairIds = np.array(pairIds)
			
			outDir = sys.argv[5] + '/' + settings.files['rankedGeneScoreDir']

			if not os.path.exists(outDir):
				os.makedirs(outDir)
			if not os.path.exists(outDir + "/" + runId):
				os.makedirs(outDir + "/" + runId)

			
			pairScoresWithPairIds = np.empty([len(svGeneIndices), 80], dtype="object")
			pairScoresWithPairIds[:,0] = pairIds
			pairScoresWithPairIds[:,1:80] = pairScores
			
			np.savetxt(outDir + '/' + runId + "/nonCoding_geneSVPairs.txt_" + str(permutationRound), pairScoresWithPairIds, delimiter='\t', fmt='%s')
			
			#Also output the coding pairs
			codingPairs = []
			for gene in genes:
				for sv in gene.SVs:

					codingPairs.append(gene.name + "_" + sv)
			codingPairs = np.array(codingPairs, dtype="object")
			np.savetxt(outDir + '/' + runId + "/coding_geneSVPairs.txt_" + str(permutationRound), codingPairs, delimiter='\t', fmt='%s')		
			
			###3. Make the feature file for MIL for each sv-gene pair
			#Each SV-gene pair is a bag. A bag can contain a variable set of isntances, which represent the gained/lost elements
			#The feature vector was pre-defined in the gene class for each instance.

			bags = dict()
			for geneInd in range(0, genes.shape[0]):
				gene = reverseGeneMap[geneInd] #Get the gene back from the scoring matrix by index
				
				geneSVSamples = []
				for sv in gene.SVs:
					splitSV = sv.split("_")
					geneSVSamples.append(splitSV[len(splitSV)-1])
					
				for sv in gene.alteredElements:
					
					#filter out SNV/CNV effects
					splitSV = sv.split("_")
					sample = splitSV[6]
					
					instances = []
					for element in gene.alteredElements[sv]:
						instances.append(gene.alteredElements[sv][element])
					
					if len(instances) > 0:
						bags[gene.name + "_" + sv] = instances
				
			#print(bags)
			print(len(bags))
		
			#output the bags to a file
			with open(outDir + '/' + runId + '/bags.pkl', 'wb') as handle:
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
				sampleName = splitSV[len(splitSV)-2]
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
				geneSVSamples.append(splitSV[len(splitSV)-2])

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
					
					#check if there is a CNV in the patient. 
					if settings.general['nonCoding'] == True and settings.general['cnvs'] == False: 

						if sample in gene.CNVs:
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
				geneSVSamples.append(splitSV[len(splitSV)-2])

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
					
					if settings.general['nonCoding'] == True and settings.general['cnvs'] == False: 

						if sample in gene.CNVs:
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
				#geneSVSamples.append(splitSV[len(splitSV)-1])
				geneSVSamples.append(sv)
			
			
			if len(gene.lostElementsSVs) > 0:
				
				for sv in gene.lostElementsSVs:
					pairId = gene.name + "_" + sv
					
					geneSVSamples = []
					for geneSV in gene.SVs:
						splitSV = geneSV.split("_")
						geneSVSamples.append(splitSV[len(splitSV)-1])
					
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					# if settings.general['nonCoding'] == True and settings.general['coding'] == False:
					# 	
					# 	splitSV = sv.split("_")
					# 
					# 	#if splitSV[len(splitSV)-2] in geneSVSamples:
					# 	#	continue
					# 
					# 	if splitSV[7] != 'INV' and splitSV[7] != 'DUP':
					# 		
					# 		if sv in geneSVSamples:
					# 			continue
					# 
					#check if we need to ignore a pair if it has an SNV in that patient.  
					# if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 
					# 	
					# 	splitSV = sv.split("_")
					# 	patientID = splitSV[6]
					# 	if patientID in gene.SNVs:
					# 		continue
					# 
					# if settings.general['nonCoding'] == True and settings.general['cnvs'] == False: 
					# 	
					# 	splitSV = sv.split("_")
					# 	patientID = splitSV[6]
					# 	if patientID in gene.CNVs:
					# 		continue
					# 
					
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
						#pairScores[pairId] = len(gene.lostElementsSVs[sv][element])
						
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
						#geneSVSamples.append(splitSV[len(splitSV)-1])
						geneSVSamples.append(geneSV)
						
						#print('gene SVs: ', geneSV)
						
					#If mutually exclusive mode, skip genes that also have coding SVs in the same sample. 
					if settings.general['nonCoding'] == True and settings.general['coding'] == False:
						
						splitSV = sv.split("_")
						
						#check properly. If DUP or INV, do not exclude the
						#genes within that SV, because these are affected.
						

						#if splitSV[len(splitSV)-2] in geneSVSamples:
						#	continue
					
						#do not exclude inv or dup genes. 
						# if splitSV[7] != 'INV' and splitSV[7] != 'DUP':
						# 	
						# 	if sv in geneSVSamples:
						# 		continue
						# 
					#check if we need to ignore a pair if it has an SNV in that patient.  
					# if settings.general['nonCoding'] == True and settings.general['snvs'] == False: 
					# 	
					# 	splitSV = sv.split("_")
					# 	patientID = splitSV[6]
					# 	if patientID in gene.SNVs:
					# 		print('skipping, SNV')
					# 		continue
					# 	
					# if settings.general['nonCoding'] == True and settings.general['cnvs'] == False: 
					# 	
					# 	splitSV = sv.split("_")
					# 	patientID = splitSV[6]
					# 	if patientID in gene.CNVs:
					# 		print('skipping, cNV')
					# 		continue	
					# 	
					if pairId not in svGeneMap: #Check if we already used this index for a different feature
						svInd = len(svGeneIndices)
						svGeneMap[pairId] = svInd
						svGeneIndices.append(pairId)
					else:
						svInd = svGeneMap[pairId]
				
					gain = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.gainedElementsSVs[sv]:
						if element == elementType:
							gain = True
					if gain == True:
						pairScores[pairId] = 1 #assume that each SV can disrupt a gene only once
						#pairScores[pairId] = len(gene.gainedElementsSVs[sv][element])
					
					
		
		return pairScores, svGeneMap, svGeneIndices	

	def scoreByElementGainsStrengthsSVs(self, genes, elementType):
		
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			if len(gene.gainedElementsStrengthsSVs) > 0:
				
				for sv in gene.gainedElementsStrengthsSVs:
					pairId = gene.name + "_" + sv
					
					for element in gene.gainedElementsStrengthsSVs[sv]:
						if element == elementType:
							
							pairScores[pairId] = gene.gainedElementsStrengthsSVs[sv][element]
		return pairScores 
		
	def scoreByElementLossesStrengthsSVs(self, genes, elementType):
		
		pairScores = dict() #Keep a dictionary because we do not know how large the final matrix will be across all features
		for geneInd in range(0, len(genes)):
			gene = genes[geneInd]
			
			if len(gene.lostElementsStrengthsSVs) > 0:
				
				for sv in gene.lostElementsStrengthsSVs:
					pairId = gene.name + "_" + sv
					
					for element in gene.lostElementsStrengthsSVs[sv]:
						if element == elementType:
							
							pairScores[pairId] = gene.lostElementsStrengthsSVs[sv][element]
		return pairScores 	
