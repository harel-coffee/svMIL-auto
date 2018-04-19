class Gene:
	"""
		Class to describe a gene. Will hold all other information related to the neighborhood of the gene as well, like the TADs, eQTLs and SVs. 
		
	"""
	
	def __init__(self, name, chromosome, start, end):
		
		self.name = name
		self.chromosome = chromosome
		self.start = start
		self.end = end
		
		
	def setTADs(self, leftTAD, rightTAD):
		
		self.leftTAD = leftTAD
		self.rightTAD = rightTAD
		
	def setEQTLS(self, eQTLS):
		
		self.eQTLS = eQTLS
		
	def setSVs(self, SVs):
		
		self.SVs = SVs