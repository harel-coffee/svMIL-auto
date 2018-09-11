class Interaction:
	
	"""
		Class to describe the interactions between regions and the SVs/SNVs overlapping it. 
		In the current implementation, region 1 is always the non-coding region, and region 2 is the region containing a gene. 
	"""
	
	def __init__(self, chromosome1, start1, end1, chromosome2, start2, end2, region2Gene):
		
		self.chromosome1 = chromosome1
		self.start1 = start1
		self.end1 = end1
		self.chromosome2 = chromosome2
		self.start2 = start2
		self.end2 = end2
		self.region2Gene = region2Gene #This is currently not a gene object, but only the name
		self.SVs = None
		self.SNVs = None
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
	
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
	