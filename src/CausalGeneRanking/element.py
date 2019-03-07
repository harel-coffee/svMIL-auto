class Element:
	
	"""
		Generic class for genomic elements that can be described with a chromosome, start and end position. E.g. enhancers, eQTLs, promoters, DNAse I sites, TFs .... 
		We assume that these elements can be linked to genes, but that is not necessary. 
		
	"""
	
	
	def __init__(self, chromosome, start, end):
		
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.SVs = []
		self.SNVs = []
		self.genes = [] #Make sure that we also know which gene(s) the element associates with
		self.type = "element" #Default type. This is set depending on which type of genomic element we are looking at. 
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def addGene(self, gene):
		self.genes.append(gene)
