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
	
		cancerTypeTotalSVs = dict() #save a number of how many SVs are there in total for each cancer type to make the final scoring matrix. 
		for gene in genes:
			
			if gene.SVs is not None:
				for sv in gene.SVs:
					cancerType = sv.cancerType
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					#If the gene is not yet in the list of cancer types, add it
					if gene not in cancerTypes[cancerType]:
						cancerTypes[cancerType][gene] = []
					#Also add all SVs to that gene that are relevant for this cancer type.
					
					cancerTypes[cancerType][gene].append(sv)
					cancerTypeTotalSVs[cancerType] += 1

		#For each gene, get the SVs.
		#Only get the SVs of the right cancer type (can be different per gene)
					
		for cancerType in cancerTypes:
			print "cancer type: ", cancerType
			cancerTypeSVs = cancerTypes[cancerType] #Use these SVs to map to the right position in the scoring matrix
			#print cancerTypeSVs
			#For each cancer type, loop through the genes.
			#Define the scoring matrix
			
			scoringMatrix = np.zeros([cancerTypeTotalSVs[cancerType], len(cancerTypeSVs)])
			genesInOrder = [] #keep these lists to ensure that after we rank we can get back to the genes that are most interesting. 
			svsInOrder = []
			for geneInd in range(0, len(cancerTypeSVs)): #first loop over the genes, then get their SVs
				
				gene = cancerTypeSVs.keys()[geneInd]
				genesInOrder.append(gene)
				#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
				
				#Perform additional check to see if the gene has SVs at all
				
				geneSVs = cancerTypeSVs[gene]
				
				for svInd in range(0, len(geneSVs)):
					
					#If the SV is not in the list of SVs, it is from a different cancer type, do not score it. 

					scoringMatrix[svInd][geneInd] = 1
					svsInOrder.append(geneSVs[svInd])
		
			#print scoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.
			
			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]
		
			
			
			for geneInd in sortedGenesInd:
				
				
				if geneScoresSummed[geneInd] > 20:
					print genesInOrder[geneInd].name, genesInOrder[geneInd].chromosome, genesInOrder[geneInd].start, ": ", geneScoresSummed[geneInd]	
		
			#exit()
			
				
					