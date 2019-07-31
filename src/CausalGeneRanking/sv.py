class SV:
	"""
		Class describing SV objects. 
	"""
	
	def __init__(self, chr1, s1, e1, o1, chr2, s2, e2, o2, sampleName, cancerType, svType):
		
		self.chr1 = chr1
		self.s1 = s1
		self.e1 = e1
		self.s2 = s2
		self.e2 = e2
		self.chr2 = chr2
		self.sampleName = sampleName
		self.cancerType = cancerType
		self.svType = svType
		self.o1 = o1
		self.o2 = o2
	