import numpy as np

class GeneRankingSampleBased:
	"""
		Class responsible for ranking genes by their causality given the SVs in their neighborhood, incorporating mutual exclusivity (SVs in TADs are more relevant if it does not overlap the gene as well)
		
		- Get all samples of a specific cancer type
		- FOr each sample, get all causal genes. Make a gene by sample matrix (1 gene, multiple samples)
		- For each causal gene, assign a score of 1 if it is affected
		- Repeat for the left and right neighboring TAD
		- Then use XOR to only have the SVs that overlap either the gene or the TAD
		- Sum to get a total score for the gene.
		- Rank across all genes, which ones have the highest score?
		
		
	"""
	
	def __init__(self, genes):
		
		
		
	