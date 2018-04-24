import numpy as np

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
	

	#This part of the code is still very sloppy, it can be distributed into different functions much better, but for now it is ok just to test
	def __init__(self, genes):
		
		#1. Get all unique cancer types and map the gene objects to the right cancer type
		cancerTypes = dict()
	
		cancerTypeTotalSVs = dict() #save a number of how many SVs are there in total for each cancer type to make the final scoring matrix. 
		svMap = dict() #SVs and their index in the final scoring matrix. 
		geneMap = dict() #genes adn their index
		geneIndex = 0 #Keep index for where the genes will be stored in the scoring matrix. 
		svIndex = 0
		
		for gene in genes:
			if gene not in geneMap:
				geneMap[gene] = geneIndex
				geneIndex += 1
				
			
			
			if gene.SVs is not None:
				for sv in gene.SVs:
					
					if sv not in svMap:
						svMap[sv] = svIndex
						svIndex += 1
					
					cancerType = sv.cancerType
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					#If the gene is not yet in the list of cancer types, add it
					if gene not in cancerTypes[cancerType]:
						cancerTypes[cancerType]["genes"][gene] = []	
					
					#Also add all SVs to that gene that are relevant for this cancer type.
					
					cancerTypes[cancerType]["genes"][gene].append(sv)
					cancerTypeTotalSVs[cancerType] += 1
			
			#Also add the SVs affecting TADs to the total list of SVs and the SV map
			
			#Get the left TAD
			leftTAD = gene.leftTAD
			
			
			if leftTAD is not None:
				for sv in leftTAD.SVs:
					if sv not in svMap:
						svMap[sv] = svIndex
						svIndex += 1
					cancerType = sv.cancerType
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if leftTAD not in cancerTypes[cancerType]["TADs"]:
						cancerTypes[cancerType]["TADs"][leftTAD] = []	
					
					cancerTypes[cancerType]["TADs"][leftTAD].append(sv)
					cancerTypeTotalSVs[cancerType] += 1
				
			rightTAD = gene.rightTAD
			
			if rightTAD is not None:
				
				for sv in rightTAD.SVs:
					if sv not in svMap:
						svMap[sv] = svIndex
						svIndex += 1
					cancerType = sv.cancerType
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if rightTAD not in cancerTypes[cancerType]["TADs"]:
						cancerTypes[cancerType]["TADs"][rightTAD] = []	
					
					cancerTypes[cancerType]["TADs"][rightTAD].append(sv)
					cancerTypeTotalSVs[cancerType] += 1
					
		
		#For each gene, get the SVs.
		#Only get the SVs of the right cancer type (can be different per gene)
					
		for cancerType in cancerTypes:
			print "cancer type: ", cancerType
			cancerTypeSVs = cancerTypes[cancerType] #Use these SVs to map to the right position in the scoring matrix
			#print cancerTypeSVs
			#For each cancer type, loop through the genes.
			#Define the scoring matrix
			
			scoringMatrix = np.zeros([len(svMap), len(geneMap)])
			geneScoringMatrix = self.scoreBySVsInGenes(scoringMatrix, cancerTypeSVs["genes"], svMap, geneMap)
			tadScoringMatrix = self.scoreBySVsInTADs(scoringMatrix, cancerTypeSVs, svMap, geneMap)
			#Combine the scoring matrices
			
			#First use xor to remove entries where both the TADs and genes are affected by the same SVs
			xorMatrix = np.logical_xor(geneScoringMatrix, tadScoringMatrix).astype(int)
			affectedPosGenes =  np.where(geneScoringMatrix > 0)
			affectedPosTads = np.where(tadScoringMatrix > 0)
			
			# XOR would be good if genes and TADs are always affected by the same SV. If the gene is already disrupted, it does not really matter anymore that the TAD is gone as well. Now we are counting this double,
			# which may not be necessary. 
			# print affectedPosGenes
			# print affectedPosTads
			# 
			# print np.setdiff1d(affectedPosGenes, affectedPosTads)
			# 
			# print xorMatrix
			# exit()
			#geneXorMatrix = geneScoringMatrix * xorMatrix
			#tadXorMatrix = tadScoringMatrix * xorMatrix
			
			#scoringMatrix = geneXorMatrix + tadXorMatrix
			scoringMatrix = geneScoringMatrix + tadScoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.
			
			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]

			for geneInd in sortedGenesInd:

				gene = geneMap.keys()[geneMap.values().index(geneInd)]
			
				if geneScoresSummed[geneInd] > 0:
					print gene.name, gene.chromosome, gene.start, ": ", geneScoresSummed[geneInd]	
		
			exit()
	
	
	#The issue that we run into here is that there are different SVs, but mapping between indices of all SVs and genes is rather tough. 		
	def scoreBySVsInGenes(self, scoringMatrix, cancerTypeSVs, svMap, geneMap):
		#In this function, loop through the genes and SVs and give the right score.  

		for geneInd in range(0, len(cancerTypeSVs)): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs.keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
			
			#Perform additional check to see if the gene has SVs at all
			
			geneSVs = cancerTypeSVs[gene]
			
			for svInd in range(0, len(geneSVs)):
				
				#If the SV is not in the list of SVs, it is from a different cancer type, do not score it. 
				matrixSvInd = svMap[geneSVs[svInd]]
				
				scoringMatrix[matrixSvInd][matrixGeneInd] = 1
		
		return scoringMatrix
		
	#Make sure that the scoring matrix has the same dimensions
	#However here we need to make sure that the filtered set of SVs is also applicable to the TADs. 
	def scoreBySVsInTADs(self, scoringMatrix, cancerTypeSVs, svMap, geneMap):
		
		#Get the left and right TAD for each gene.
		#The score is 1 for each TAD boundary affected.
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#1. Check the left and right TAD
			
			leftTAD = gene.leftTAD
			rightTAD = gene.rightTAD
			
			#If the TAD does not have any SVs in that cancer type, skip it.
			if leftTAD in cancerTypeSVs["TADs"]:
				leftTADSVs = cancerTypeSVs["TADs"][leftTAD]
				
				#2. Add a score of 1 for every affected TAD.
			
				for sv in leftTADSVs:
					
					matrixSvInd = svMap[sv]
					
					scoringMatrix[matrixSvInd][geneInd] = 1
			
			if rightTAD in cancerTypeSVs["TADs"]:
				rightTADSVs = cancerTypeSVs["TADs"][rightTAD]

				for sv in rightTADSVs:
					
					matrixSvInd = svMap[sv]
					
					scoringMatrix[matrixSvInd][geneInd] = 1
			
		return scoringMatrix
		
		
		
		
		
		
		1+1
		
					