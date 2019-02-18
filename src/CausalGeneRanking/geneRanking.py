import numpy as np
import settings

class GeneRanking:
	"""
		Class responsible for ranking genes by their causality given the SVs in their neighborhood.
		
		Applying the following rules:
		
		For each cancer type, gather all samples in which SVs have been detected.
		
		- All scores are stored in a SV x gene matrix, to allow SVs to affect multiple genes. 
		- If a gene is directly affected by an SV, the score is 1.
		- If the boundary of the left or right TAD of a gene is directly affected by an SV, the score is 1
		- For each sample of a cancer type where we compute these scores, we take the sum across the samples.
		- The genes are then ranked by their final score. The gene with the highest score is most likely causal.
	
		TO DO:
		- Because SNVs were added later, many variables still refer to SVs, but this is actually 'variants' in general and can also contain SNVs. Fix this!
		- Split into more functions. Also large parts of the code can be re-used
		- Taking the distance to TADs into account, perhaps incorporate CTCF sites
		- Interactions/heat diffusion is broken
	
	"""
	
	#This part of the code is still very sloppy, it can be distributed into different functions much better, but for now it is ok just to test
	def __init__(self, genes, svData, mode):
		"""
			TO DO:
			- Document
			- Move code out of constructor and into proper functions
		"""
		
		sampleMap = dict() #samples and their index in the final scoring matrix. 
		geneMap = dict() #genes adn their index
		reverseGeneMap = dict() #also keep a map where we can search by index to later obtain back the gene from the matrix. 
		scores = dict()
		#1. Get all unique cancer types and map the gene objects to the right cancer type
		cancerTypes = np.unique(svData[:,6])
		#2. Make the sample map
		samples = np.unique(svData[:,7])
		for sampleInd in range(0, len(samples)):
			sampleMap[samples[sampleInd]] = sampleInd
		
		#Make the gene maps
		geneIndex = 0
		for gene in genes:
			if gene not in geneMap:
				geneMap[gene] = geneIndex
				geneIndex += 1
		for gene in geneMap:
			index = geneMap[gene]
			reverseGeneMap[index] = gene
		
		#For each gene, get the SVs and SNVs.
		#Only get the variants of the right cancer type (can be different per gene)
		
		#Then do the scoring for the variants in each data type in the neighborhood individually
		
		print "doing the scoring"
		
		for cancerType in cancerTypes:
			print "current cancer type: ", cancerType

			
			
			#Do the scoring of the genes
			#We make a scoring matrix of patients x genes. Each gene has a score in each patient of if an SV overlaps with that element in the neighborhood of the gene yes/no.
			#To get the total score for a gene, we can sum across all patients. And then we need to filter the matrix to make sure that we only sum the scores for the genes if
			#the element in the neighborhood is disrupted by a variant, but the gene itself is not. 
			print "scoring genes:"
			geneScoringMatrix = self.scoreBySVsInGenes(genes, sampleMap, geneMap)
			print "scoring eQTLs"
			#eQTLScoringMatrix = self.scoreBySVsInEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			#eQTLScoringMatrix = self.scoreByEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			eQTLGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, cancerType, "eQTL")
			eQTLLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, cancerType, "eQTL")
			enhancerGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, cancerType, "enhancer")
			enhancerLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, cancerType, "enhancer")
			
			#Scoring only by lost interactions
			scoringMatrix = enhancerGainsScoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.

			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort by highest final score and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]


			#Now map the indices of the scoring matrix back to the actual genes, and report the scores in the different layers per gene. 
			geneScores = []
			for geneInd in sortedGenesInd:
				
				gene = geneMap.keys()[geneMap.values().index(geneInd)] #Isn't this the part that is going wrong? The best index is probably the index in the matrix? 
				gene = reverseGeneMap[geneInd]
				
				
				geneScores.append([gene, np.sum(geneScoringMatrix[:,geneInd]), np.sum(eQTLGainsScoringMatrix[:,geneInd]), np.sum(eQTLLossesScoringMatrix[:,geneInd]), np.sum(enhancerGainsScoringMatrix[:,geneInd]), np.sum(enhancerLossesScoringMatrix[:,geneInd])])
			
			
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
			
			
			
		self.scores = scores #make it global for now because I of course can't return from here in the constructor.... When everything here is a proper function, this should be fixed. 
			
			
			#It would be nice to do the scoring steps above in a function as well
			
	
	
	#First compute how many SVs overlap with each gene itself. These will later be used to filter for elements that have SVs that do not overlap with the gene itself.		
	def scoreBySVsInGenes(self, cancerTypeSVs, sampleMap, geneMap):
		#In this function, loop through the genes and SVs and give the right score.
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		checkedSampleInds = []
		
		for geneInd in range(0, len(cancerTypeSVs)): #first loop over the genes, then get their SVs
			gene = cancerTypeSVs[geneInd]
			
			
			matrixGeneInd = geneMap[gene] #index of the gene in the scoring matrix
			
			#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
			
			#Perform additional check to see if the gene has SVs at all
			
			geneSVs = gene.SVs
			if geneSVs is None:
				continue
			
			geneUniqueSamples = dict()
			for svInd in range(0, len(geneSVs)):
			
				
				#check the sample of this sv and get the right position of the samples in the scoring matrix
				sampleName = geneSVs[svInd][6]
				sampleInd = sampleMap[sampleName]
				
				
				scoringMatrix[sampleInd][matrixGeneInd] = 1
	
		return scoringMatrix
	
	def scoreByElementGains(self, genes, sampleMap, geneMap, cancerType, elementType):
		"""
			For every gene, add a score of 1 if an eQTL is either gained or lost. Later separate losses from gains.
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]
			
	
			
			matrixGeneInd = geneMap[gene]
			
			if len(gene.gainedElements) > 0:
				for sample in gene.gainedElements:
					gain = False
					for element in gene.gainedElements[sample]:
						if element == elementType:
							gain = True
					if gain == True:	#Make sure that we only count every sample once	
						sampleInd = sampleMap[sample]
						scoringMatrix[sampleInd][matrixGeneInd] += 1
				
		
		return scoringMatrix
	
	def scoreByElementLosses(self, genes, sampleMap, geneMap, cancerType, elementType):
		"""
			For every gene, add a score of 1 if an eQTL is either gained or lost. Later separate losses from gains.
			
			Instead of counting the number of losses, normalize for the number of losses compared to the total. 
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			
			if len(gene.lostElements) > 0:
				for sample in gene.lostElements:
					
					loss = False
					for element in gene.lostElements[sample]:
						if element == elementType:
							loss = True
					if loss == True:	#Make sure that we only count every sample once	
						sampleInd = sampleMap[sample]
						scoringMatrix[sampleInd][matrixGeneInd] += 1
					
		
		return scoringMatrix
