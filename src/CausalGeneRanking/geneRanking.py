import numpy as np
import settings

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
	def __init__(self, genes, mode):
		
		#Check if the genes are always annotated int he same way.Where is the random factor?
		
		#1. Get all unique cancer types and map the gene objects to the right cancer type
		cancerTypes = dict()
	
		cancerTypeTotalSVs = dict() #save a number of how many SVs are there in total for each cancer type to make the final scoring matrix. 
		sampleMap = dict() #samples and their index in the final scoring matrix. 
		geneMap = dict() #genes adn their index
		geneIndex = 0 #Keep index for where the genes will be stored in the scoring matrix. 
		sampleIndex = 0 #Index for the samples in the scoring matrix
		reverseGeneMap = dict() #also keep a map where we can search by index to later obtain back the gene from the matrix. 
		
		scores = dict() #save the final scores per gene per cancer type.
		
		
		print "ordering genes by cancer types"
		for gene in genes:
			samplesAndSVCounts = dict()
			if gene not in geneMap:
				geneMap[gene] = geneIndex
				geneIndex += 1
				
			
			#First collect all variants of the gene, both SVs and SNVs
			geneVariants = [] #have one big list of objects
			if gene.SVs is not None:
				geneVariants += list(gene.SVs)
			if gene.SNVs is not None:
				geneVariants += list(gene.SNVs)
				
			
				
			if len(geneVariants) > 0:
				
				for variant in geneVariants:
					
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
						
					if variant[6] not in samplesAndSVCounts:
						samplesAndSVCounts[variant[6]] = 0
					samplesAndSVCounts[variant[6]] += 1
					
					cancerType = variant[7]
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
						cancerTypes[cancerType]["interactions"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					#If the gene is not yet in the list of cancer types, add it
					if gene not in cancerTypes[cancerType]["genes"]:
						cancerTypes[cancerType]["genes"][gene] = []	
					
					#Also add all SVs to that gene that are relevant for this cancer type.
					
					cancerTypes[cancerType]["genes"][gene].append(variant)
					cancerTypeTotalSVs[cancerType] += 1
				
			#Also add the SVs affecting TADs to the total list of SVs and the SV map
			
			#print "mapping left TAD SVs to cancer type"
			
			#Get the left TAD
			leftTAD = gene.leftTAD
			
			leftTADVariants = [] #have one big list of objects
			if leftTAD is not None:
				if leftTAD.SVs is not None:
					leftTADVariants += list(leftTAD.SVs)
				if leftTAD.SNVs is not None:
					leftTADVariants += list(leftTAD.SNVs)
				
			
			if len(leftTADVariants) > 0:
				for variant in leftTADVariants:
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
					cancerType = variant[7]
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
						cancerTypes[cancerType]["interactions"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if leftTAD not in cancerTypes[cancerType]["TADs"]:
						cancerTypes[cancerType]["TADs"][leftTAD] = []	
					
					cancerTypes[cancerType]["TADs"][leftTAD].append(variant)
					cancerTypeTotalSVs[cancerType] += 1
			
			#print "mapping right tad svs to cancer types"
				
			rightTAD = gene.rightTAD
			rightTADVariants = [] #have one big list of objects
			if rightTAD is not None:
				
				if rightTAD.SVs is not None:
					rightTADVariants += list(rightTAD.SVs)
				if rightTAD.SNVs is not None:
					rightTADVariants += list(rightTAD.SNVs)
				
			if len(rightTADVariants) > 0:
				
				for variant in rightTADVariants:
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
					cancerType = variant[7]
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
						cancerTypes[cancerType]["interactions"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if rightTAD not in cancerTypes[cancerType]["TADs"]:
						cancerTypes[cancerType]["TADs"][rightTAD] = []	
					
					cancerTypes[cancerType]["TADs"][rightTAD].append(variant)
					cancerTypeTotalSVs[cancerType] += 1
		
			#Get the eQTLs and add the SVs to the right cancer type
			
			#print "mapping eQTLs to cancer type"
			
			eQTLs = gene.eQTLs
			
			
			for eQTL in eQTLs:
				
				#Collect all variants of this eQTL	
				eQTLVariants = [] #have one big list of objects
				if eQTL.SVs is not None:
					eQTLVariants += list(eQTL.SVs)
				if eQTL.SNVs is not None:
					eQTLVariants += list(eQTL.SNVs)
					
				
				for variant in eQTLVariants:
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
						
					cancerType = variant[7]
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
						cancerTypes[cancerType]["interactions"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if eQTL not in cancerTypes[cancerType]["eQTLs"]:
						cancerTypes[cancerType]["eQTLs"][eQTL] = []	
					
					cancerTypes[cancerType]["eQTLs"][eQTL].append(variant)
					cancerTypeTotalSVs[cancerType] += 1

			print "mapping interactions to cancer types"
			
			for interaction in gene.interactions:
				
				#Collect all variants of this eQTL	
				interactionVariants = [] #have one big list of objects
				if interaction.SVs is not None:
					interactionVariants += list(interaction.SVs)
				if interaction.SNVs is not None:
					interactionVariants += list(interaction.SNVs)
					
				
				for variant in interactionVariants:
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
						
					cancerType = variant[7]
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["TADs"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
						cancerTypes[cancerType]["interactions"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if interaction not in cancerTypes[cancerType]["interactions"]:
						cancerTypes[cancerType]["interactions"][interaction] = []	
					
					cancerTypes[cancerType]["interactions"][interaction].append(variant)
					cancerTypeTotalSVs[cancerType] += 1
		
		for gene in geneMap:
			index = geneMap[gene]
			reverseGeneMap[index] = gene
		
		#For each gene, get the SVs.
		#Only get the SVs of the right cancer type (can be different per gene)
		
		print "doing the scoring"
		print cancerTypes.keys()
		for cancerType in cancerTypes:
			print "current cancer type: ", cancerType
			# if mode == "SV":
			# 	if cancerType != "breast/gastric": #restrict to one cancer type for now
			# 		continue
			# if mode == "SNV":
			# 	if cancerType != "breast": #Use a different cancer type to filter for in the SNVs, these come form different patients. 
			# 		continue
			#
			
			if cancerType != "breast": #focus on one cancer type for now
				continue
			
			print "cancer type: ", cancerType
			cancerTypeSVs = cancerTypes[cancerType] #Use these SVs to map to the right position in the scoring matrix
			
			#For each cancer type, loop through the genes.
			#Define the scoring matrix

			geneScoringMatrix = self.scoreBySVsInGenes(cancerTypeSVs["genes"], sampleMap, geneMap)
			eQTLScoringMatrix = self.scoreBySVsInEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			tadScoringMatrix = self.scoreBySVsInTADs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			interactionScoringMatrix = self.scoreBySVsInInteractions(cancerTypeSVs, sampleMap, geneMap, cancerType)

			#? Does this part make sense, to first xor with eQtls and then again later?
			xorMatrix = np.logical_xor(geneScoringMatrix, eQTLScoringMatrix).astype(int)
			geneXorMatrix = geneScoringMatrix * xorMatrix
			eQTLXorMatrix = eQTLScoringMatrix * xorMatrix
			tadXorMatrix = tadScoringMatrix * xorMatrix
			interactionXorMatrix = interactionScoringMatrix * xorMatrix

			scoringMatrix = geneXorMatrix + eQTLXorMatrix + tadXorMatrix  + interactionXorMatrix
			
			#exit()
			#Sum the total score per gene and report the genes by which ones are most likely causal.

			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]

			geneCount = 0
			geneScores = []
			for geneInd in sortedGenesInd:
				
				gene = geneMap.keys()[geneMap.values().index(geneInd)] #Isn't this the part that is going wrong? The best index is probably the index in the matrix? 
				gene = reverseGeneMap[geneInd]
				
				
				if geneScoresSummed[geneInd] > 0:
					geneCount += 1
					#print gene.name, gene.chromosome, gene.start, ": ", geneScoresSummed[geneInd], " gene score: ", np.sum(geneXorMatrix[:,geneInd]), " eQTL score: ", np.sum(eQTLXorMatrix[:,geneInd]), " TAD score: ", np.sum(tadXorMatrix[:,geneInd])	
			
				
				geneScores.append([gene, np.sum(geneXorMatrix[:,geneInd]), np.sum(eQTLXorMatrix[:,geneInd]), np.sum(tadXorMatrix[:,geneInd]), np.sum(interactionXorMatrix[:,geneInd])])
			
			
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
			#print "total genes: ", geneCount
			
			
			
		self.scores = scores #make it global for now because I can't return from here. When everything here is a proper function, this should be fixed. 
			
			
			#It would be nice to do the scoring steps above in a function as well
			
	
	
	#The issue that we run into here is that there are different SVs, but mapping between indices of all SVs and genes is rather tough. 		
	def scoreBySVsInGenes(self, cancerTypeSVs, sampleMap, geneMap):
		#In this function, loop through the genes and SVs and give the right score.
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		checkedSampleInds = []
		
		for geneInd in range(0, len(cancerTypeSVs)): #first loop over the genes, then get their SVs
			gene = cancerTypeSVs.keys()[geneInd]
			
			
			matrixGeneInd = geneMap[gene]
			
			#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
			
			#Perform additional check to see if the gene has SVs at all
			
			geneSVs = cancerTypeSVs[gene]
		
			
			geneUniqueSamples = dict()
			for svInd in range(0, len(geneSVs)):
			
				
				#check the sample of this sv
				sampleName = geneSVs[svInd][6]
				sampleInd = sampleMap[sampleName]
				
				
					
				#If the SV is not in the list of SVs, it is from a different cancer type, do not score it. 
				#matrixSvInd = svMap[geneSVs[svInd]]
				
				scoringMatrix[sampleInd][matrixGeneInd] += 1
	
		return scoringMatrix
		
	#Make sure that the scoring matrix has the same dimensions
	#However here we need to make sure that the filtered set of SVs is also applicable to the TADs. 
	def scoreBySVsInTADs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		#Get the left and right TAD for each gene.
		#The score is 1 for each TAD boundary affected.
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			
			matrixGeneInd = geneMap[gene]
			
			
			#1. Check the left and right TAD
			
			leftTAD = gene.leftTAD
			rightTAD = gene.rightTAD
			
			# if gene.name == "MUC1":
			# 	print leftTAD.chromosome, leftTAD.start, leftTAD.end
			# 	print rightTAD.chromosome, rightTAD.start, rightTAD.end
			# 
			
			#If the TAD does not have any SVs in that cancer type, skip it.
			if leftTAD in cancerTypeSVs["TADs"]:
				leftTADSVs = cancerTypeSVs["TADs"][leftTAD]
				
				#2. Add a score of 1 for every affected TAD.
				
				# if gene.name == "MUC1":
				# 	print "Number of SVs in left TAD: "
				# 	print len(leftTADSVs)
				# 
				for sv in leftTADSVs:
					sampleInd = sampleMap[sv[6]]
					if sv[7] == cancerType:
						scoringMatrix[sampleInd][matrixGeneInd] += 1
			
			if rightTAD in cancerTypeSVs["TADs"]:
				rightTADSVs = cancerTypeSVs["TADs"][rightTAD]
				# 
				# if gene.name == "MUC1":
				# 	print "Number of SVs in right TAD: "
				# 	print len(rightTADSVs)

					
				for sv in rightTADSVs:
					
					sampleInd = sampleMap[sv[6]]
					
					if sv[7] == cancerType:
					
						scoringMatrix[sampleInd][matrixGeneInd] += 1
		
		return scoringMatrix
							
	def scoreBySVsInEQTLs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		#Get the left and right TAD for each gene.
		#The score is 1 for each TAD boundary affected.
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		foundSamples = []
		sampleCount = dict()
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#Score SVs in the eQTLs
			eQTLs = gene.eQTLs
			
			for eQTL in eQTLs:
				
				
				if eQTL in cancerTypeSVs["eQTLs"]: #only score this eQTL if it has overlapping SVs in this particular cancer type
					
					for sv in eQTL.SVs:
						
						if sv[7] == cancerType:
							
							if gene.name == "CNBD1":
							
								
								if sv[6] not in sampleCount:
									sampleCount[sv[6]] = 0
								sampleCount[sv[7]] += 1
							
									
								foundSamples.append([sv[7], sv[6], sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2])
							
							sampleInd = sampleMap[sv[6]]
							scoringMatrix[sampleInd][matrixGeneInd] += 1
			
		
			
			#Divide the score at each position in the matrix by the total number of eQTLs mapped to this gene. It would be ideal to do this on the final summed matrix, but is more difficult because of the combination.
			for row in range(0, scoringMatrix.shape[0]):
				if scoringMatrix[row][matrixGeneInd] > 0:
					scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd] / len(eQTLs)
	
						
		return scoringMatrix
		
	def scoreBySVsInEQTLs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		#Get the left and right TAD for each gene.
		#The score is 1 for each TAD boundary affected.
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		foundSamples = []
		sampleCount = dict()
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#Score SVs in the eQTLs
			eQTLs = gene.eQTLs
			
			for eQTL in eQTLs:
				
				
				if eQTL in cancerTypeSVs["eQTLs"]: #only score this eQTL if it has overlapping SVs in this particular cancer type
					
					for sv in eQTL.SVs:
						
						if sv[7] == cancerType:
							
							if gene.name == "CNBD1":
							
								
								if sv[6] not in sampleCount:
									sampleCount[sv[6]] = 0
								sampleCount[sv[7]] += 1
							
									
								foundSamples.append([sv[7], sv[6], sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2])
							
							sampleInd = sampleMap[sv[6]]
							scoringMatrix[sampleInd][matrixGeneInd] += 1
			
		
			
			#Divide the score at each position in the matrix by the total number of eQTLs mapped to this gene. It would be ideal to do this on the final summed matrix, but is more difficult because of the combination.
			for row in range(0, scoringMatrix.shape[0]):
				if scoringMatrix[row][matrixGeneInd] > 0:
					scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd] / len(eQTLs)
	
						
		return scoringMatrix
		
	def scoreBySVsInEQTLs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		#Get the left and right TAD for each gene.
		#The score is 1 for each TAD boundary affected.
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		foundSamples = []
		sampleCount = dict()
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#Score SVs in the eQTLs
			eQTLs = gene.eQTLs
			
			for eQTL in eQTLs:
				
				
				if eQTL in cancerTypeSVs["eQTLs"]: #only score this eQTL if it has overlapping SVs in this particular cancer type
					
					for sv in eQTL.SVs:
						
						if sv[7] == cancerType:
							
							if gene.name == "CNBD1":
							
								
								if sv[6] not in sampleCount:
									sampleCount[sv[6]] = 0
								sampleCount[sv[7]] += 1
							
									
								foundSamples.append([sv[7], sv[6], sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2])
							
							sampleInd = sampleMap[sv[6]]
							scoringMatrix[sampleInd][matrixGeneInd] += 1
			
		
			
			#Divide the score at each position in the matrix by the total number of eQTLs mapped to this gene. It would be ideal to do this on the final summed matrix, but is more difficult because of the combination.
			for row in range(0, scoringMatrix.shape[0]):
				if scoringMatrix[row][matrixGeneInd] > 0:
					scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd] / len(eQTLs)
	
						
		return scoringMatrix
		
	#Scoring specific for interactions	
	def scoreBySVsInInteractions(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		#Here it goes wrong because we go through the genes that are directly affected by SVs in some cancer type, but the TADs may not be affected in that cancer type. 
		foundSamples = []
		sampleCount = dict()
		genePositions = []
		for geneInd in range(0, len(cancerTypeSVs["genes"])): #first loop over the genes, then get their SVs
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			genePositions.append(['chr' + gene.chromosome, gene.start, gene.end, gene.name, gene])
			
			matrixGeneInd = geneMap[gene]
			
			#Score SVs in the interactions
			
			for interaction in gene.interactions:
				
				
				if interaction in cancerTypeSVs["interactions"]: #only score this eQTL if it has overlapping SVs in this particular cancer type
					
					for sv in interaction.SVs:
						
						if sv[7] == cancerType:
							
							if gene.name == "CNBD1":
							
								
								if sv[6] not in sampleCount:
									sampleCount[sv[6]] = 0
								sampleCount[sv[7]] += 1
							
									
								foundSamples.append([sv[7], sv[6], sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2])
							
							sampleInd = sampleMap[sv[6]]
							scoringMatrix[sampleInd][matrixGeneInd] += 1

			#Normalize for the total number of interactions
			for row in range(0, scoringMatrix.shape[0]):
				if scoringMatrix[row][matrixGeneInd] > 0:
					scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd] / len(gene.interactions)
	
		#These are the initial heat scores per gene, now we need to overlap these with the diffusion scores
		#The heat * diffusion scores will be the final scores for the genes that we use to rank.
		
		genePositions = np.array(genePositions, dtype='object')
		
		heatDiffusionScoresMatrix = self.computeHeatDiffusionScoresInteractions(cancerTypeSVs, scoringMatrix, geneMap, genePositions)
		
		return heatDiffusionScoresMatrix
	
	def computeHeatDiffusionScoresInteractions(self, cancerTypeSVs, svScoresMatrix, geneMap, genePositions):
		"""
			Obtain the heat diffusion scores from the file for each region.
			Then overlap each gene with the regions, and obtain the correct score.
			If a gene overlaps with multiple regions, take the sum of the scores of those regions.
			Multiply the diffusion score with the SV-based score for that gene.
			Return the scoring matrix with the multiplied scores. 
		"""

		diffusionScores = self.readHiCDiffusionScores(genePositions)
		
		
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#get the position of the gene
			if gene in diffusionScores: #not all genes have overlap with a region for some reason
				diffusionScore = diffusionScores[gene] 
			
				svScoresMatrix[:,matrixGeneInd] = svScoringMatrix[:,matrixGeneInd] * diffusionScore #multiply the entire column for that gene with the diffusion score to get the total heat score
			
		return svScoresMatrix


	def readHiCDiffusionScores(self, genePositions):
		"""
			Read the diffusion scores from a file. We can initialize the scoring matrix here for the genes, containing the diffusion scores at each gene. 
		"""
		
		#get the diffusion scores from the file
		print genePositions
	
		
		heatDiffusionFile = settings.files['heatDiffusionFile']
		diffusionScores = dict()
		#Read the file into a numpy array, split the regions
		with open(heatDiffusionFile, 'r') as inFile:
			
			for line in inFile:
				line = line.strip()
				
				splitLine = line.split("\t")
				region = splitLine[0]

				splitRegion = region.split("_")
				chrom = splitRegion[0]
				
				#Check if this region contains any gene. If true, find which gene it is and get the position in the diffusion matrix. Fill this with the diffusion score.
				
				regionStart = int(splitRegion[1])
				regionEnd = regionStart + 5000 #should be a setting, bin size of HiC data. 
				
				#print regionStart
				#print regionEnd
				
				#match only on the right chromosome
				
				chrSubset = genePositions[np.where(genePositions[:,0] == chrom)]
				
				
				#Start of the gene should be before the end of the Hi-C region, and the end of the gene should be after the start of the region. 
				matchingGenesStart = regionStart < chrSubset[:,2]
				matchingGenesEnd = regionEnd > chrSubset[:,1]
				
				matchingGenesInd = matchingGenesStart * matchingGenesEnd
				
				matchingGenes = chrSubset[matchingGenesInd]
				
				
				
				for matchingGene in matchingGenes: #for the genes in the region, get their position in the diffusion matrix and assign the diffusion score.
					
					diffusionScores[matchingGene[3]] = splitLine[1]
		
		
		return diffusionScores		
		
		
		
		
		