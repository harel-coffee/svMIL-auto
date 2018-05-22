class EQTL:
	
	"""
		Class to descirbe eQTLs and also information such as if these are annotated as enhancer, or if SVs overlap with this eQTL.
		
		Information about which gene this eQTL affects is part of the Gene class. 
		
	"""
	
	
	def __init__(self, chromosome, start, end):
		
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.enhancerStatus = "N" #y default not an enhancer
		self.SVs = []
		self.SNVs = []
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def setEnhancerStatus(self, enhancerStatus): #This funcion can later be replaced to include other regulatory element types as well
		self.enhancerStatus = enhancerStatus
