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
					else:
						cancerTypes[cancerType].append(gene)
		
		
		
		for cancerType in cancerTypes:
			
			cancerTypeGenes = cancerTypes[cancerType] 
			
			#First check to see if the SVs are really set for all elements
			print "gene SVs:"
			for sv in gene.SVs:
				print sv.chr1, sv.s1, sv.e1
				
			print "right TAD SVs:"
			for sv in gene.rightTAD.SVs:
				print sv.chr1, sv.s1, sv.e1
				
			print "left TAD SVs:"
			for sv in gene.leftTAD.SVs:
				print sv.chr1, sv.s1, sv.e1
			
			
			#For the causal genes, first:
			#1. Make a subset of SVs in a specific cancer type
			
			#2. For each sample, if there are SVs overlapping the gene, the score is 1 for that gene and that SV
			
			
			
			