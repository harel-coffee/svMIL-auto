import numpy as np

from tad import TAD
from sv import SV
from gene import Gene
from eQTL import EQTL


class NeighborhoodDefiner:
	"""
		Class responsible for defining the neighborhood of causal genes.
		
		Currently, the neighborhood consists of:
		
		- nearest TADs on the left and right of the gene
		- all eQTLs mapped to the gene (and if these are enhancers or not)
		- SVs overlapping either the gene directly, or other elements in the neighborhood
	
	"""

	
	def __init__(self, genes, variantData, mode):
		
		#1. Map genes to TADs
		
		#Make these pats a setting!
		tadFile = "../../data/tads/tads.csv"
		print "Getting TADs"
		tadData = self.getTADsFromFile(tadFile)
		print "mapping TADs to genes"
		self.mapTADsToGenes(genes[:,3], tadData) #only pass the gene objects will suffice
		
		#2. Map genes to eQTLs
		
		eQTLFile = "../../data/eQTLs/eQTLsFilteredForCausalGenes.txt"
		print "getting eQTLs"
		eQTLData = self.getEQTLsFromFile(eQTLFile, genes[:,3])
		
		#3. Map SVs to all neighborhood elements
		if mode == "SV":
			print "Mapping SVs to the neighborhood"
			self.mapSVsToNeighborhood(genes, variantData)
		if mode == "SNV":
			print "Mapping SNVs to the neighborhood"
			self.mapSNVsToNeighborhood(genes, variantData)
		
	def getTADsFromFile(self, tadFile):
		"""
			Read the TADs into NumPy format. I don't read the TADs into object immediately, because the NumPy matrices work very well for
			quick overlapping. I add an object reference to the matrix so that we can later add the right TAD object to the genes. 
		"""
		
		#Read the gene list data into a list
		tadData = []
		with open(tadFile, "r") as f:
			lineCount = 0
			for line in f:
				if lineCount < 2: #skip header
					lineCount += 1
					continue
				line = line.strip()
				splitLine = line.split("\t")
				
				
				#Convert the numbers to integers for quicker comparison. 0 is chromosome, 1 is start, 2 is end. Celltypes are not used for now. 
				splitLine[1] = int(splitLine[1])
				splitLine[2] = int(splitLine[2])
				
				TADObject = TAD(splitLine[0], splitLine[1], splitLine[2])
				
				#chr, start, end
				tadData.append([splitLine[0], splitLine[1], splitLine[2], TADObject])
		
		#Also convert the other dataset to numpy
		tadData = np.array(tadData, dtype='object')
		
		return tadData
		
		
	def getEQTLsFromFile(self, eQTLFile, genes):
		
		eQTLs = []
		with open(eQTLFile, 'rb') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				eQTLObject = EQTL(splitLine[0], int(splitLine[1]), int(splitLine[2])) #chr, start, end
				
				#The mapping information is in the file, so we can already do it here
				self.mapEQTLsToGenes(eQTLObject, genes, splitLine[3])
						
				eQTLs.append(eQTLObject)				
		
		return eQTLs
			
		
	def mapTADsToGenes(self, genes, tadData):
		"""
			Adds the left and right TAD to each gene object.
		"""
		print tadData
		#For each gene, check which TADs are on the correct chromosome.
		
		#Then see which ones are directly on the left and directly on the right of the gene.
		previousChr = None
		for gene in genes:
			
			#1. Make a subset of TADs on the right chromosome. There should be only 1 chromosome
			
			if "chr" + str(gene.chromosome) != previousChr:
				
				#Find the two subsets that match on both chromosomes. 
				matchingChrInd = tadData[:,0] == "chr" + str(gene.chromosome)
				
				#It is not even necessary to make 2 lists if both chromosomes are the same, we could use a reference without re-allocating
				chrSubset = tadData[np.where(matchingChrInd)]
				
				#Make sure to update the previous chromosome when it changes
				previousChr = "chr" + str(gene.chromosome)
				
			if np.size(chrSubset) < 1:
				continue #no need to compute the distance, there are no genes on these chromosomes
			
			#Within this subset, check which TADs are on the right and left of the current gene
			
			#TADs on the left have an end coordinate smaller than the gene start coordinate. Similar for the right TADs
			
			leftTADs = chrSubset[np.where(chrSubset[:,2] <= gene.start)]
			rightTADs = chrSubset[np.where(chrSubset[:,1] > gene.end)]
			
			#Compute the distance to each of these TADs and take the TADs with the minimum distance
			if leftTADs.shape[0] > 0:
				leftTADsDistance = np.abs(leftTADs[:,2] - gene.start)
				nearestLeftTAD = leftTADs[np.argmin(leftTADsDistance),3]
				gene.setLeftTAD(nearestLeftTAD)
			else:
				gene.setLeftTAD(None)
			if rightTADs.shape[0] > 0:
				rightTADsDistance = np.abs(rightTADs[:,1] - gene.end)
				nearestRightTAD = rightTADs[np.argmin(rightTADsDistance),3]
				gene.setRightTAD(nearestRightTAD)
			else:
				gene.setRightTAD(None)

			
			
		
	def mapEQTLsToGenes(self, eQTL, genes, geneSymbol):
		"""
			Map the right gene object to the eQTL object. 
		
		"""
		for gene in genes:
			
			if gene.name == geneSymbol:
				
				gene.addEQTL(eQTL)
			
		
	def mapSVsToNeighborhood(self, genes, svData):
		"""
			Take as input gene objects with the neighborhood pre-set, and search through the SVs to find which SVs overlap the genes, TADs and eQTLs in the gene neighborhood
		
			This function will need to be properly split into multiple functions. 
		
		"""
		
		previousChr = None
		for gene in genes[:,3]:
			
			#Make a subset of SVs that either overlap on their chromosome 1 or chromosome 2
			
			if gene.chromosome != previousChr:
				
				matchingSvChr1Ind = svData[:,0] == str(gene.chromosome)
				matchingSvChr2Ind = svData[:,3] == str(gene.chromosome)
				
				#Intra and chr1 and chr2 will overlap if we don't exclude the positions where both chr1 and chr2 are the same. 
		
				notChr1Matches = svData[:,3] != str(gene.chromosome)
				chr1OnlyInd = matchingSvChr1Ind * notChr1Matches
				
				notChr2Matches = svData[:,0] != str(gene.chromosome)
				chr2OnlyInd = matchingSvChr2Ind * notChr2Matches
				
				#intraSubset: chr1 and chr2 both match
				matchingChr1AndChr2Ind = matchingSvChr1Ind * matchingSvChr2Ind
				intraSubset = svData[matchingChr1AndChr2Ind]
				
				#interChr1Subset: only chr1 matches
				interChr1Subset = svData[chr1OnlyInd]
				
				#interChr2Subset: only chr2 matches
				interChr2Subset = svData[chr2OnlyInd]
				
		
				#Now concatenate them into one set, but keep the formatting the same as: chr, start, end
				
				svChr1Subset = np.empty([interChr1Subset.shape[0],11],dtype='object')
				svChr1Subset[:,0] = interChr1Subset[:,0] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,1] = interChr1Subset[:,1] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,2] = interChr1Subset[:,2] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,3] = interChr1Subset[:,7] #Also keep the sample name
				
				#Keep other regional info as well just to check
				svChr1Subset[:,4] = interChr1Subset[:,0] #Also keep the sample name
				svChr1Subset[:,5] = interChr1Subset[:,1] #Also keep the sample name
				svChr1Subset[:,6] = interChr1Subset[:,2] #Also keep the sample name
				svChr1Subset[:,7] = interChr1Subset[:,3] #Also keep the sample name
				svChr1Subset[:,8] = interChr1Subset[:,4] #Also keep the sample name
				svChr1Subset[:,9] = interChr1Subset[:,5] #Also keep the sample name
				svChr1Subset[:,10] = interChr1Subset[:,8]
				
				#Make the subset for the chr2 matches
				svChr2Subset = np.empty([interChr2Subset.shape[0],11], dtype='object')
				svChr2Subset[:,0] = interChr2Subset[:,0] #For chromosome 2, we use just the second chromosome, s2 and e2.
				svChr2Subset[:,1] = interChr2Subset[:,4] 
				svChr2Subset[:,2] = interChr2Subset[:,5] 
				svChr2Subset[:,3] = interChr2Subset[:,7] #Also keep the sample name
				
				#Keep other regional info as well just to check
				svChr2Subset[:,4] = interChr2Subset[:,0] #Also keep the sample name
				svChr2Subset[:,5] = interChr2Subset[:,1] #Also keep the sample name
				svChr2Subset[:,6] = interChr2Subset[:,2] #Also keep the sample name
				svChr2Subset[:,7] = interChr2Subset[:,3] #Also keep the sample name
				svChr2Subset[:,8] = interChr2Subset[:,4] #Also keep the sample name
				svChr2Subset[:,9] = interChr2Subset[:,5] #Also keep the sample name
				svChr2Subset[:,10] = interChr2Subset[:,8]
				
				#For the intra subset, we need to use s1 and e2.
				svIntraSubset = np.empty([intraSubset.shape[0],11], dtype='object')
				svIntraSubset[:,0] = intraSubset[:,0] #For chromosome 2, we use chromosome 1, s1 and e2.
				svIntraSubset[:,1] = intraSubset[:,1] 
				svIntraSubset[:,2] = intraSubset[:,5] 
				svIntraSubset[:,3] = intraSubset[:,7] #Also keep the sample name
				
				#Keep other regional info as well just to check
				svIntraSubset[:,4] = intraSubset[:,0] #Also keep the sample name
				svIntraSubset[:,5] = intraSubset[:,1] #Also keep the sample name
				svIntraSubset[:,6] = intraSubset[:,2] #Also keep the sample name
				svIntraSubset[:,7] = intraSubset[:,3] #Also keep the sample name
				svIntraSubset[:,8] = intraSubset[:,4] #Also keep the sample name
				svIntraSubset[:,9] = intraSubset[:,5] #Also keep the sample name
				svIntraSubset[:,10] = intraSubset[:,8]
				
				#Now concatenate the arrays
				svSubset = np.concatenate((svChr1Subset, svChr2Subset, svIntraSubset))
			
				previousChr = gene.chromosome
			
			if np.size(svSubset) < 1:
				continue #no need to compute the distance, there are no TADs on this chromosome
			
			#Check which of these SVs overlap with the gene itself
			
			geneStartMatches = gene.start <= svSubset[:,2]
			geneEndMatches = gene.end >= svSubset[:,1]
			
			geneMatches = geneStartMatches * geneEndMatches #both should be true, use AND for concatenating
		
			svsOverlappingGenes = svSubset[geneMatches,10]
			
			
			#Get the SV objects and link them to the gene
			gene.setSVs(svsOverlappingGenes)
			
			#Check which SVs overlap with the right/left TAD
			
			if gene.leftTAD != None:
				
				leftTADStartMatches = gene.leftTAD.start <= svSubset[:,2]
				leftTADEndMatches = gene.leftTAD.end >= svSubset[:,1]
				
				leftTADMatches = leftTADStartMatches * leftTADEndMatches
				
				svsOverlappingLeftTAD = svSubset[leftTADMatches, 10]
				gene.leftTAD.setSVs(svsOverlappingLeftTAD)
			
			if gene.rightTAD != None:
				
				rightTADStartMatches = gene.rightTAD.start <= svSubset[:,2]
				rightTADEndMatches = gene.rightTAD.end >= svSubset[:,1]
				
				rightTADMatches = rightTADStartMatches * rightTADEndMatches
			
				svsOverlappingRightTAD = svSubset[rightTADMatches, 10]
				gene.rightTAD.setSVs(svsOverlappingRightTAD)
			
			#Check which SVs overlap with the eQTLs
			
			#Repeat for eQTLs. Is the gene on the same chromosome as the eQTL? Then use the above chromosome subset.
			
			geneEQTLs = gene.eQTLs
			
			for eQTL in geneEQTLs: #only if the gene has eQTLs
				
				
				
				startMatches = eQTL.start <= svSubset[:,2]
				endMatches = eQTL.end >= svSubset[:,1]
				
				allMatches = startMatches * endMatches
				
				
				
				svsOverlappingEQTL = svSubset[allMatches, 10]

				
				eQTL.setSVs(svsOverlappingEQTL)

		
		
		
	def mapSNVsToNeighborhood(self, genes, snvData):
		"""
			Same as the function for mapping SVs to the neighborhood, but then for SNVs.
			Will also need to be properly split into multiple functions, many pieces of code can be re-used. 
		
		"""
		
		previousChr = None
		for gene in genes[:,3]:
			
			#Make a subset of SNVs on the right chromosome
			
			if gene.chromosome != previousChr:
				
				matchingChrInd = snvData[:,0] == str(gene.chromosome)

				snvSubset = snvData[matchingChrInd]
				
				previousChr = gene.chromosome
			
			if np.size(snvSubset) < 1:
				continue #no need to compute the distance, there are no TADs on this chromosome
			
			#Check which of these SNVs overlap with the gene itself
			
			geneStartMatches = gene.start <= snvSubset[:,2]
			geneEndMatches = gene.end >= snvSubset[:,1]
			
			geneMatches = geneStartMatches * geneEndMatches #both should be true, use AND for concatenating
		
			#Get the SNV objects of the overlapping SNVs
			
			snvsOverlappingGenes = []
			for match in geneMatches:
				snvObject = SNV(match[0], match[1], match[2], match[3], match[4])
				snvsOverlappingGenes.append(snvObject)
			
			#snvsOverlappingGenes = snvSubset[geneMatches,5]
			
			
			#Get the SV objects and link them to the gene
			gene.setSNVs(snvsOverlappingGenes)
			
			#Check which SVs overlap with the right/left TAD
			
			if gene.leftTAD != None:
				
				leftTADStartMatches = gene.leftTAD.start <= snvSubset[:,2]
				leftTADEndMatches = gene.leftTAD.end >= snvSubset[:,1]
				
				leftTADMatches = leftTADStartMatches * leftTADEndMatches
				
				snvsOverlappingLeftTAD = []
				gene.leftTAD.setSNVs(snvsOverlappingLeftTAD)
			
			if gene.rightTAD != None:
				
				rightTADStartMatches = gene.rightTAD.start <= snvSubset[:,2]
				rightTADEndMatches = gene.rightTAD.end >= snvSubset[:,1]
				
				rightTADMatches = rightTADStartMatches * rightTADEndMatches
			
				snvsOverlappingRightTAD = []
				gene.rightTAD.setSNVs(snvsOverlappingRightTAD)
			
			#Check which SVs overlap with the eQTLs
			
			#Repeat for eQTLs. Is the gene on the same chromosome as the eQTL? Then use the above chromosome subset.
			
			geneEQTLs = gene.eQTLs
			
			for eQTL in geneEQTLs: #only if the gene has eQTLs
				
				
				
				startMatches = eQTL.start <= snvSubset[:,2]
				endMatches = eQTL.end >= snvSubset[:,1]
				
				allMatches = startMatches * endMatches
				
				
				
				snvsOverlappingEQTL = []

				
				eQTL.setSNVs(snvsOverlappingEQTL)

		
		
		
		
		
		
		
		
	