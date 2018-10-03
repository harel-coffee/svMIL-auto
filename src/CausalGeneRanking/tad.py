class TAD:
	
	"""
		Class to describe the location of a TAD and potentially SVs overlapping it.
	"""
	
	def __init__(self, chromosome, start, end):
		
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.SVs = None
		self.SNVs = None
		self.interactions = None
		self.eQTLInteractions = None
		
	def setSVs(self, SVs): #All SVs that overlap with this TAD boundary (left or right depending on the gene. Is that safe?)
		
		self.SVs = SVs
	
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def setInteractions(self, interactions): #All interactions that take place within this TAD
		self.interactions = interactions
		
	def setEQTLInteractions(self, eQTLInteractions): #All eQTL interactions that take place within the TAD. 
		self.eQTLInteractions = eQTLInteractions 
	