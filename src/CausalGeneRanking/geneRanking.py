import numpy as np

class GeneRanking:
	"""
		Class responsible for ranking genes by their causality given the SVs in their neighborhood.
		
		Applying the following rules:
		
		For each cancer type, gather all samples in which SVs have been detected.
		
		- All scores are stored in a SV x gene matrix, to allow SVs to affect multiple genes. 
		- If a gene is directly affected by an SV, the score is 1.
		- If the boundary of the left or right TAD of a gene is directly affected by an SV, the score is 1 * 1/the distance to the TAD.
		- For each sample of a cancer type where we compute these scores, we take the sum across the samples.
		- The genes are then ranked by their final score. The gene with the highest score is most likely causal.
	
	"""
	

	def __init__(self, genes):
		
		
		
		
		#1. Get all unique cancer types and map the gene objects to the right cancer type
		cancerTypes = dict()
		
		for gene in genes:
			
			if gene.SVs is not None:
				for sv in gene.SVs:
					cancerType = sv.cancerType
					if cancerType not in cancerTypes.keys():
						cancerTypes[cancerType] = []
					
					cancerTypes[cancerType].append(sv)
		
		for cancerType in cancerTypes:
			cancerTypeSVs = cancerTypes[cancerType] #Use these SVs to map to the right position in the scoring matrix
			print cancerTypeSVs
			#For each cancer type, loop through the genes.
			#Define the scoring matrix
			
			scoringMatrix = np.empty([len(cancerTypeSVs), len(genes)])
			
			for geneInd in range(0, len(genes)):
				
				gene = genes[geneInd]
				
				#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
				
				#Perform additional check to see if the gene has SVs at all
				
				geneSVs = gene.SVs
				
				for sv in geneSVs:
					svInd = cancerTypeSVs.index(sv) ##Somehow the SVs here are not in the list, I don't understand why. 
					
					scoringMatrix[svInd][geneInd] = 1
					
					
			print scoringMatrix
					