import numpy as np
import settings

class GeneRanking:
	"""
		Class responsible for ranking genes by their causality given the SVs in their neighborhood.
		
		Applying the following rules:
		
		For each cancer type, gather all samples in which SVs have been detected.
		
		- All scores are stored in a SV x gene matrix, to allow SVs to affect multiple genes. 
		- If a gene is directly affected by an SV, the score is 1.
		- If the boundary of the left or right TAD of a gene is directly affected by an SV, the score is 1
		- For each sample of a cancer type where we compute these scores, we take the sum across the samples.
		- The genes are then ranked by their final score. The gene with the highest score is most likely causal.
	
		TO DO:
		- Because SNVs were added later, many variables still refer to SVs, but this is actually 'variants' in general and can also contain SNVs. Fix this!
		- Split into more functions. Also large parts of the code can be re-used
		- Taking the distance to TADs into account, perhaps incorporate CTCF sites
		- Interactions/heat diffusion is broken
	
	"""
	
	#This part of the code is still very sloppy, it can be distributed into different functions much better, but for now it is ok just to test
	def __init__(self, genes, mode):
		"""
			TO DO:
			- Document
			- Move code out of constructor and into proper functions
		"""
		
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
			#This way the scoring can be done at once for all variants, they are encoded in the same way. 
			geneVariants = [] #have one big list of objects
			if gene.SVs is not None:
				geneVariants += list(gene.SVs)
			if gene.SNVs is not None:
				geneVariants += list(gene.SNVs)
				
			
			#Get all the genes and make a dictionary where the variants (SVs and SNVs) are listed per cancer type, and directly under taht per data type.
			#This can be way more efficient to filter by cancer type early on, we don't need to go through all SVs and SNVs anymore in that case. 
			if len(geneVariants) > 0:
				
				for variant in geneVariants:
					
					if variant[6] not in sampleMap:
						sampleMap[variant[6]] = sampleIndex
						sampleIndex += 1
						
					if variant[6] not in samplesAndSVCounts:
						samplesAndSVCounts[variant[6]] = 0
					samplesAndSVCounts[variant[6]] += 1
					
					cancerType = variant[7]
					#Make the different data types for each cancer type, this is a bit of a strange setup
					if cancerType not in cancerTypes:
						cancerTypes[cancerType] = dict()
						cancerTypes[cancerType]["genes"] = dict()
						cancerTypes[cancerType]["eQTLs"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0 #counting of SVs, not used anymore
					
					#If the gene is not yet in the list of cancer types, add it
					if gene not in cancerTypes[cancerType]["genes"]:
						cancerTypes[cancerType]["genes"][gene] = []	
					
					#Also add all SVs to that gene that are relevant for this cancer type.
					
					cancerTypes[cancerType]["genes"][gene].append(variant)
					cancerTypeTotalSVs[cancerType] += 1
				
			#Also add the SVs affecting TADs to the total list of SVs and the SV map
			
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
						cancerTypes[cancerType]["eQTLs"] = dict()
					if cancerType not in cancerTypeTotalSVs:
						cancerTypeTotalSVs[cancerType] = 0
					
					if eQTL not in cancerTypes[cancerType]["eQTLs"]:
						cancerTypes[cancerType]["eQTLs"][eQTL] = []	
					
					cancerTypes[cancerType]["eQTLs"][eQTL].append(variant)
					cancerTypeTotalSVs[cancerType] += 1

		for gene in geneMap:
			index = geneMap[gene]
			reverseGeneMap[index] = gene
		
		#For each gene, get the SVs and SNVs.
		#Only get the variants of the right cancer type (can be different per gene)
		
		#Then do the scoring for the variants in each data type in the neighborhood individually
		
		print "doing the scoring"
		print cancerTypes.keys()
		for cancerType in cancerTypes:
			print "current cancer type: ", cancerType
			
			
			if cancerType != "breast": #focus on one cancer type for now, will later be a setting and we need to make sure to filter the variants early on if we restrict to cancer types. 
				continue
			
			print "cancer type: ", cancerType
			cancerTypeSVs = cancerTypes[cancerType] #Get all variants (not just SVs anymore!!! update naming) found overlapping with a neighborhood element in this cancer type. 
			
			#Do the scoring of the genes
			#We make a scoring matrix of patients x genes. Each gene has a score in each patient of if an SV overlaps with that element in the neighborhood of the gene yes/no.
			#To get the total score for a gene, we can sum across all patients. And then we need to filter the matrix to make sure that we only sum the scores for the genes if
			#the element in the neighborhood is disrupted by a variant, but the gene itself is not. 
			print "scoring genes:"
			geneScoringMatrix = self.scoreBySVsInGenes(cancerTypeSVs["genes"], sampleMap, geneMap)
			print "scoring eQTLs"
			#eQTLScoringMatrix = self.scoreBySVsInEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			#eQTLScoringMatrix = self.scoreByEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			eQTLGainsScoringMatrix = self.scoreByEQTLGains(cancerTypeSVs, sampleMap, geneMap, cancerType)
			eQTLLossesScoringMatrix = self.scoreByEQTLLosses(cancerTypeSVs, sampleMap, geneMap, cancerType)
			
			# 
			# print "scoring gained eQTL interactions:"
			# eQTLGainedInteractionScoringMatrix = self.scoreByGainedEQTLInteractions(cancerTypeSVs["genes"], sampleMap, geneMap)
			
			#interactionScoringMatrix = self.scoreBySVsInInteractions(cancerTypeSVs, sampleMap, geneMap, cancerType)

			#? Does this part make sense, to first xor with eQtls and then again later?
			#I think we need to do this for all data types right? Because all of them need to be exclusive? 
			#The idea is that we end up with one scoring matrix where only the scores remain for the elements in the neighborhood where SVs in those patients do not overlap with the gene itself. 
			
			#Filter out all patients in which the SV also disrupts the gene. 
			eQTLGainsXorMatrix = np.logical_xor(geneScoringMatrix, eQTLGainsScoringMatrix).astype(int)
			
			#Actually we don't want to use XOR, we only want a score of 1 if the losses are 1, not if the genes are 1 but the losses are 0
			eQTLLossesXorMatrix = np.zeros(eQTLLossesScoringMatrix.shape)
			for row in range(0, eQTLLossesScoringMatrix.shape[0]):
				for col in range(0, eQTLLossesScoringMatrix.shape[1]):
					
					if eQTLLossesScoringMatrix[row][col] == 1:
						if geneScoringMatrix[row][col] == 0:
							eQTLLossesXorMatrix[row][col] = 1
					
			#eQTLLossesXorMatrix = np.logical_xor(eQTLLossesScoringMatrix, geneScoringMatrix).astype(int)
			
			
			#Scoring only by lost interactions
			scoringMatrix = eQTLLossesScoringMatrix
			
			#Sum the total score per gene and report the genes by which ones are most likely causal.

			geneScoresSummed = np.sum(scoringMatrix, axis=0)
			
			#Sort by highest final score and report the names of the genes that are involved
			sortedGenesInd = np.argsort(geneScoresSummed)[::-1]


			#Now map the indices of the scoring matrix back to the actual genes, and report the scores in the different layers per gene. 
			geneCount = 0
			geneScores = []
			for geneInd in sortedGenesInd:
				
				gene = geneMap.keys()[geneMap.values().index(geneInd)] #Isn't this the part that is going wrong? The best index is probably the index in the matrix? 
				gene = reverseGeneMap[geneInd]
				
				if gene.name == "BRCA1":
					print eQTLGainsScoringMatrix[:,geneInd]
					print eQTLLossesScoringMatrix[:,geneInd]
					print eQTLLossesXorMatrix[:,geneInd]
					print geneScoringMatrix[:,geneInd]
					print sampleMap
					
				
				
				if geneScoresSummed[geneInd] > 0:
					geneCount += 1
					#print gene.name, gene.chromosome, gene.start, ": ", geneScoresSummed[geneInd], " gene score: ", np.sum(geneXorMatrix[:,geneInd]), " eQTL score: ", np.sum(eQTLXorMatrix[:,geneInd]), " TAD score: ", np.sum(tadXorMatrix[:,geneInd])	
			
				
				geneScores.append([gene, np.sum(geneScoringMatrix[:,geneInd]), np.sum(eQTLGainsScoringMatrix[:,geneInd]), np.sum(eQTLLossesScoringMatrix[:,geneInd])])
			
			
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
			#print "total genes: ", geneCount
			
			
			
		self.scores = scores #make it global for now because I of course can't return from here in the constructor.... When everything here is a proper function, this should be fixed. 
			
			
			#It would be nice to do the scoring steps above in a function as well
			
	
	
	#First compute how many SVs overlap with each gene itself. These will later be used to filter for elements that have SVs that do not overlap with the gene itself.		
	def scoreBySVsInGenes(self, cancerTypeSVs, sampleMap, geneMap):
		#In this function, loop through the genes and SVs and give the right score.
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		checkedSampleInds = []
		
		for geneInd in range(0, len(cancerTypeSVs)): #first loop over the genes, then get their SVs
			gene = cancerTypeSVs.keys()[geneInd]
			
			
			matrixGeneInd = geneMap[gene] #index of the gene in the scoring matrix
			
			#1. Check which genes are directly overlapped by an SV, these get a score of 1 for these SVs.
			
			#Perform additional check to see if the gene has SVs at all
			
			geneSVs = cancerTypeSVs[gene]
		
			
			geneUniqueSamples = dict()
			for svInd in range(0, len(geneSVs)):
			
				
				#check the sample of this sv and get the right position of the samples in the scoring matrix
				sampleName = geneSVs[svInd][6]
				sampleInd = sampleMap[sampleName]
				
				
				scoringMatrix[sampleInd][matrixGeneInd] += 1
	
		return scoringMatrix
	
	def scoreByGainedEQTLInteractions(self, genes, sampleMap, geneMap):
		"""
			For each gene, get how many gains of interactions are in each sample, and make the scoring matrix of genes by samples. 
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])

		for gene in genes:
			
			geneInd = geneMap[gene]
			
			for sample in gene.gainedEQTLs: #From the number of eQTLs alone we do not know the patients. 
				sampleInd = sampleMap[sample]
				scoringMatrix[sampleInd][geneInd] = len(gene.gainedEQTLs[sample])
			
			
			
		
		
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

	def scoreByEQTLs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		"""
			For every gene, add a score of 1 if an eQTL is either gained or lost. Later separate losses from gains.
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			if len(gene.lostEQTLs) > 0:
				for sample in gene.lostEQTLs:
					sampleInd = sampleMap[sample]
					scoringMatrix[sampleInd][matrixGeneInd] += 1
		
			if len(gene.gainedEQTLs) > 0:
				for sample in gene.gainedEQTLs:
					sampleInd = sampleMap[sample]
					scoringMatrix[sampleInd][matrixGeneInd] += 1
				
		
		return scoringMatrix
	
	def scoreByEQTLGains(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		"""
			For every gene, add a score of 1 if an eQTL is either gained or lost. Later separate losses from gains.
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			if gene.name == "PTK6":
				print "Gained eQTLs of PTK6: ", len(gene.gainedEQTLs)
			
			matrixGeneInd = geneMap[gene]
			
			if len(gene.gainedEQTLs) > 0:
				for sample in gene.gainedEQTLs:
					if gene.name == "PTK6":
						print "sample: ", sample
						for eQTL in gene.gainedEQTLs[sample]:
							print eQTL.start
					sampleInd = sampleMap[sample]
					scoringMatrix[sampleInd][matrixGeneInd] += 1
				
		
		return scoringMatrix
	
	def scoreByEQTLLosses(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		"""
			For every gene, add a score of 1 if an eQTL is either gained or lost. Later separate losses from gains.
			
			Instead of counting the number of losses, normalize for the number of losses compared to the total. 
		"""
		
		scoringMatrix = np.zeros([len(sampleMap), len(geneMap)])
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			if gene.name == "PTK6":
				print "Lost eQTLs of PTK6: ", len(gene.lostEQTLs)
				for eQTL in gene.lostEQTLs:
					print eQTL
		
			if len(gene.lostEQTLs) > 0:
				for sample in gene.lostEQTLs:
					delta = len(gene.eQTLs)
					sampleInd = sampleMap[sample]
					for eQTL in gene.lostEQTLs:
						delta -= 1
					#scoringMatrix[sampleInd][matrixGeneInd] = delta
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
					#scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd] / len(eQTLs)
					#Don't normalize, this should be captured by background in the permutations!
					scoringMatrix[row][matrixGeneInd] = scoringMatrix[row][matrixGeneInd]
						
		return scoringMatrix
		
	#Scoring specific for interactions	(Broken, somehow the scores are still high for interactions even when we turn it off, so ignore this part for now)
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
		
		originalRanking = []
		diffusionRanking = [] 
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			gene = cancerTypeSVs["genes"].keys()[geneInd]
			
			matrixGeneInd = geneMap[gene]
			
			#get the position of the gene
			
			if gene.name in diffusionScores: #not all genes have overlap with a region for some reason
				diffusionScore = diffusionScores[gene.name]
				
				svScoresMatrix[np.where(svScoresMatrix[:,matrixGeneInd] == 0),matrixGeneInd] = 1
				
				originalRanking.append([gene, sum(svScoresMatrix[:,matrixGeneInd])])
				diffusionRanking.append([gene, sum(svScoresMatrix[:,matrixGeneInd] * diffusionScore)])
				
				if sum(svScoresMatrix[:,matrixGeneInd] * diffusionScore) > sum(svScoresMatrix[:,matrixGeneInd]):
					print "gene ", gene.name, " has a higher diffusion score ", sum(svScoresMatrix[:,matrixGeneInd] * diffusionScore), " vs ", sum(svScoresMatrix[:,matrixGeneInd])


				#svScoresMatrix[:,matrixGeneInd] = svScoresMatrix[:,matrixGeneInd] * diffusionScore #multiply the entire column for that gene with the diffusion score to get the total heat score
				#temporarily turn off the diffusion scores, to see the influence
				#svScoresMatrix[:,matrixGeneInd] = svScoresMatrix[:,matrixGeneInd]
				svScoresMatrix[:,matrixGeneInd] = 0
		
		originalRanking = np.array(originalRanking)
		diffusionRanking = np.array(diffusionRanking)
		
		sortedOriginalRanking = np.sort(originalRanking)
		sortedDiffusionRanking = np.sort(diffusionRanking)
		
		#print sortedOriginalRanking
		#print sortedDiffusionRanking
		#exit()
			
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
					
					diffusionScores[matchingGene[3]] = float(splitLine[1])
		
		
		return diffusionScores		
		
		
		
		
		