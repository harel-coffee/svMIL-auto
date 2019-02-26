import numpy as np
import settings

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
	
	def __init__(self, genes, svData, mode):
		"""
			genes: (numpy array) array with the genes and their information. chr	start	end	Gene (object)
			
			TO DO:
			- svData and mode are currently not used, need to be added later if we include SNVs. 
			
		"""
		self.scoreGenes(genes, svData)
	
			
	def scoreGenes(self, genes, svData):
		"""
			Score the provided genes. Currently we score by:
				- SVs in genes
				- Gains/losses of eQTLs
				- Gains/losses of enhancers
				
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject
		"""
		
		sampleMap = dict() #samples and their index in the final scoring matrix. 
		geneMap = dict() #genes and their index in the final scoring matrix
		reverseGeneMap = dict() #also keep a map where we can search by index to later obtain back the gene from the matrix. 
		scores = dict() #final scores per gene
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
		
		
		print "doing the scoring"
		#Score per cancer type individually
		for cancerType in cancerTypes:
			print "current cancer type: ", cancerType

			#Do the scoring of the genes
			#We make a scoring matrix of patients x genes. Each gene has a score in each patient of if an SV overlaps with that element in the neighborhood of the gene yes/no.
			#To get the total score for a gene, we can sum across all patients. 
			print "scoring genes:"
			geneScoringMatrix = self.scoreBySVsInGenes(genes, sampleMap, geneMap)
			print "scoring eQTL: "
			eQTLGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "eQTL")
			eQTLLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "eQTL")
			print "scoring enhancers: "
			enhancerGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "enhancer")
			enhancerLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "enhancer")
			print "scoring promoters: "
			promoterGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "promoter")
			promoterLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "promoter")
			print "scoring cpg: "
			cpgGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "cpg")
			cpgLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "cpg")
			print "scoring tfs: "
			tfGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "tf")
			tfLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "tf")
			print "scoring hic: "
			hicGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "hic")
			hicLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "hic")
			
			print "scoring histone marks: "
			h3k9me3GainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k9me3")
			h3k9me3LossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k9me3")
			h3k4me3GainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k4me3")
			h3k4me3LossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k4me3")
			h3k27acGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k27ac")
			h3k27acLossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k27ac")
			h3k27me3GainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k27me3")
			h3k27me3LossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k27me3")
			h3k4me1GainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k4me1")
			h3k4me1LossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k4me1")
			h3k36me3GainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "h3k36me3")
			h3k36me3LossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "h3k36me3")
			
			print "scoring tfs: "
			dnaseIGainsScoringMatrix = self.scoreByElementGains(genes, sampleMap, geneMap, "dnaseI")
			dnaseILossesScoringMatrix = self.scoreByElementLosses(genes, sampleMap, geneMap, "dnaseI")
			
			
			#The final scoring matrix selected here determines how the genes will be ranked. For testing, this varies. Eventually, this will be a combined score that we rank by. 
			scoringMatrix = promoterGainsScoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.

			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort by highest final score and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]

			#Now map the indices of the scoring matrix back to the actual genes, and report the scores in the different layers per gene. 
			geneScores = []
			for geneInd in sortedGenesInd:
				gene = reverseGeneMap[geneInd] #Get the gene back from the scoring matrix by index

				geneScores.append([gene, np.sum(geneScoringMatrix[:,geneInd]), np.sum(eQTLGainsScoringMatrix[:,geneInd]), np.sum(eQTLLossesScoringMatrix[:,geneInd]),
								   np.sum(enhancerGainsScoringMatrix[:,geneInd]), np.sum(enhancerLossesScoringMatrix[:,geneInd]), np.sum(promoterGainsScoringMatrix[:,geneInd]), np.sum(promoterLossesScoringMatrix[:,geneInd]),
								   np.sum(cpgGainsScoringMatrix[:,geneInd]), np.sum(cpgLossesScoringMatrix[:,geneInd]), np.sum(tfGainsScoringMatrix[:,geneInd]), np.sum(tfLossesScoringMatrix[:,geneInd]),
								   np.sum(hicGainsScoringMatrix[:,geneInd]), np.sum(hicLossesScoringMatrix[:,geneInd]),
								   np.sum(h3k9me3GainsScoringMatrix[:,geneInd]), np.sum(h3k9me3LossesScoringMatrix[:,geneInd]), np.sum(h3k4me3GainsScoringMatrix[:,geneInd]), np.sum(h3k4me3LossesScoringMatrix[:,geneInd]),
								   np.sum(h3k27acGainsScoringMatrix[:,geneInd]), np.sum(h3k27acLossesScoringMatrix[:,geneInd]), np.sum(h3k27me3GainsScoringMatrix[:,geneInd]), np.sum(h3k27me3LossesScoringMatrix[:,geneInd]),
								   np.sum(h3k4me1GainsScoringMatrix[:,geneInd]), np.sum(h3k4me1LossesScoringMatrix[:,geneInd]), np.sum(h3k36me3GainsScoringMatrix[:,geneInd]), np.sum(h3k36me3LossesScoringMatrix[:,geneInd]),
								   np.sum(dnaseIGainsScoringMatrix[:,geneInd]), np.sum(dnaseILossesScoringMatrix[:,geneInd])])
			
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
	
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
				sampleName = sv
				sampleInd = sampleMap[sampleName]

				scoringMatrix[sampleInd][matrixGeneInd] = 1 #set the score to 1 to ensure that per sample we count only 1 SV. 
	
		return scoringMatrix
	
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
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]

			matrixGeneInd = geneMap[gene]
			
			if len(gene.gainedElements) > 0: #if the gene has gained elements
				for sample in gene.gainedElements:
					gain = False #Make sure that we only focus on gained elements of the provided type. 
					for element in gene.gainedElements[sample]:
						if element == elementType:
							gain = True
					if gain == True: 
						sampleInd = sampleMap[sample]
						scoringMatrix[sampleInd][matrixGeneInd] = 1
				
		
		return scoringMatrix
	
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
		
		for geneInd in range(0, len(genes)):
			
			gene = genes[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			if len(gene.lostElements) > 0:
				for sample in gene.lostElements:
					
					
					loss = False #Make sure that we only focus on lost elements of the provided type. 
					for element in gene.lostElements[sample]:
						if element == elementType:
							loss = True
					if loss == True:
						sampleInd = sampleMap[sample]
						scoringMatrix[sampleInd][matrixGeneInd] = 1
					
		
		return scoringMatrix
