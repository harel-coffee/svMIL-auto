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
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
	
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
	