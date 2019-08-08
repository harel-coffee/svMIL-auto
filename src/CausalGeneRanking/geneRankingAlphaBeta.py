from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import settings
import math
import matplotlib.pyplot as plt
from six.moves import range

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
		
		
		print("ordering genes by cancer types")
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
						
						
						
						cancerTypes[cancerType]["genes"][gene] = dict()
					
					#Also add all SVs to that gene that are relevant for this cancer type.
					
					#Make sure that every variant is counted only once.
					variantUniqueId = str(variant[0]) + "_" + str(variant[1]) + "_" + str(variant[2]) + "_" + str(variant[3]) + "_" + str(variant[4]) + "_" + str(variant[5])
					if variantUniqueId not in cancerTypes[cancerType]["genes"][gene]:
						cancerTypes[cancerType]["genes"][gene][variantUniqueId] = variant
					
				
			#Map the eQTLs to the cancer types
			
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
					
					if gene not in cancerTypes[cancerType]["genes"]:
						cancerTypes[cancerType]["genes"][gene] = dict()
					
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
		
		print("doing the scoring")
		print(list(cancerTypes.keys()))
		for cancerType in cancerTypes:
			print("current cancer type: ", cancerType)
			
			
			if cancerType != "breast": #focus on one cancer type for now, will later be a setting and we need to make sure to filter the variants early on if we restrict to cancer types. 
				continue
			
			print("cancer type: ", cancerType)
			cancerTypeSVs = cancerTypes[cancerType] #Get all variants (not just SVs anymore!!! update naming) found overlapping with a neighborhood element in this cancer type. 
			
			#Score the genes.
			#To compute beta, get the number of eQTLs per gene and define a sigmoid.
			#To compute alpha, count the number of SVs in the eQTL layer.
			#Multiply the alpha and beta per gene.
			
			print("scoring eQTLs")
			eQTLBasedGeneScores = self.scoreBySVsInEQTLs(cancerTypeSVs, sampleMap, geneMap, cancerType)
			
			#Sort by highest final score and report the names of the genes that are involved
			print("sorting genes by eQTL scores")
			sortedGenesInd = np.argsort(eQTLBasedGeneScores[:,0])[::-1]


			#Now map the indices of the scoring matrix back to the actual genes, and report the scores in the different layers per gene. 
			geneCount = 0
			geneScores = []
			print("obtaining scores")
			for geneInd in sortedGenesInd:
				
				gene = list(geneMap.keys())[list(geneMap.values()).index(geneInd)] #Isn't this the part that is going wrong? The best index is probably the index in the matrix? 
				gene = reverseGeneMap[geneInd]
				
				
				geneScores.append([gene, eQTLBasedGeneScores[geneInd][0], eQTLBasedGeneScores[geneInd][1], eQTLBasedGeneScores[geneInd][2], eQTLBasedGeneScores[geneInd][3]])
			
			
			geneScores = np.array(geneScores, dtype="object")
			scores[cancerType] = geneScores
			#print "total genes: ", geneCount
			
			
		print("done scoring")	
		self.scores = scores #make it global for now because I of course can't return from here in the constructor.... When everything here is a proper function, this should be fixed. 
			
			
			#It would be nice to do the scoring steps above in a function as well

	

	def sigmoid(self, x):
		#return logistic.cdf(x)
		return 1 / -(1 + math.exp(-x)) + 1 #1 / (1 + math.exp(-x)) is the normal sigmoid, do - the denominator to reverse the sigmoid, +1 to scale it back to the 0-1 y range
		
		#y = y0 + (y1 - y0)
		y = 0.2 + (1 - 0.2)
		
		#return 1 / float(1 + np.exp(-x))
	

	def computeEQTLBetaSigmoid(self, genes):
		
		#Define the sigmoid to use for beta.
		#1. Count the number of eQTLs that each gene has
		#2. Set a value of 0.2 at the maximum (how?????)
		#3. Between the minimum and maximum should be the sigmoid from 0.2 to 1
		
		eQTLCounts = []
		
		for gene in genes:
			
			#if len(gene.eQTLs) > 0:
			eQTLCounts.append(len(gene.eQTLs))
		
		#Plot the raw eQTL counts
		#plt.hist(eQTLCounts)
		#plt.show()
		
		#scaledEQTLCounts = np.interp(eQTLCounts, (np.min(eQTLCounts), np.max(eQTLCounts)), (-10, +10))
		scaledEQTLCounts = np.interp(eQTLCounts, (np.min(eQTLCounts), np.max(eQTLCounts)), (-1, +1))


		#Also plot the scaled eQTL counts
				
		testX = list(range(-10,10))
		sigmoidValues = []
		for x in testX:
			sigmoidValues.append(self.sigmoid(x))
		
		
		#Make a plot
		# import matplotlib.pyplot as plt
		#
		# plt.clf()
		# plt.plot(testX, sigmoidValues)
		# plt.show()
		# exit()

		return scaledEQTLCounts

	def scoreBySVsInEQTLs(self, cancerTypeSVs, sampleMap, geneMap, cancerType):
		
		#Temporarily do the test with the SOMs here, not sure if there should be another class for this yet
		
		#We want a matrix of the distance from each eQTL to all other eQTLs, also to the gene.
		
		
		selectedGene = None
		
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			gene = list(cancerTypeSVs["genes"].keys())[geneInd]
			
			if len(gene.eQTLs) > 0:
				selectedGene = list(cancerTypeSVs["genes"].keys())[geneInd] #take a random gene for now with eQTLs
				break
			
		
		
		
		distMat = np.zeros([len(selectedGene.eQTLs)+1, len(selectedGene.eQTLs) + 1])
		
		for eQTLInd in range(0, len(selectedGene.eQTLs)):
			eQTL1 = selectedGene.eQTLs[eQTLInd]
			for eQTLInd2 in range(eQTLInd, len(selectedGene.eQTLs)):
				
				#Compute distance to this eQTL
				dist = abs(eQTL1.start - selectedGene.eQTLs[eQTLInd2].start)
				
				distMat[eQTLInd][eQTLInd2] = dist
				distMat[eQTLInd2][eQTLInd] = dist
				
			#Also compute the distance to the gene
			geneDistStart = abs(eQTL1.start - selectedGene.start)
			geneDistEnd = abs(eQTL1.start - selectedGene.end)
				
			if geneDistStart < geneDistEnd:
				distMat[eQTLInd][len(selectedGene.eQTLs)] = geneDistStart
				distMat[len(selectedGene.eQTLs)][eQTLInd] = geneDistStart
			else:
				distMat[eQTLInd][len(selectedGene.eQTLs)] = geneDistEnd
				distMat[len(selectedGene.eQTLs)][eQTLInd] = geneDistEnd
		
		print(distMat.shape)
		
		#What do the data look like when we simply make a scatterplot?
		
		#Should be in 2D right? 
		
		#plt.scatter(distMat[:,0], distMat[:,1])
		#plt.show()
		
		#Apply the SOM to the data. What do we see? 
		
		# import sompy
		# mapsize = [20,20]
		# som = sompy.SOMFactory.build(distMat, mapsize, mask=None, mapshape='planar', lattice='rect', normalization='var', initialization='pca', neighborhood='gaussian', training='batch', name='sompy')  # this will use the default parameters, but i can change the initialization and neighborhood methods
		# som.train(n_job=1, verbose='info')  # verbose='debug' will print more, and verbose=None wont print anything
		# v = sompy.mapview.View2DPacked(50, 50, 'test')
		# v.show(som, what='codebook', which_dim=[0,1], cmap=None, col_sz=6)
		#
		
		exit()
		
		
		
		#Plot: distance from genes to eQTLs
		distances = dict()
		distancesList = []
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			gene = list(cancerTypeSVs["genes"].keys())[geneInd]
			
			start = gene.start
			end = gene.end
			for eQTL in gene.eQTLs:
				eQTLPos = eQTL.start
				
				#check for the closest distance
				#start to eQTL pos
				#end to eQTL pos
				
				startEQTLDist = abs(eQTLPos - start)
				endEQTLDist = abs(eQTLPos - end)
				if startEQTLDist < endEQTLDist:
					#if startEQTLDist not in distances:
					#	distances[startEQTLDist] = 0
					#distances[startEQTLDist] += 1
					if startEQTLDist < 200000:
						distancesList.append(startEQTLDist)
				else:
					#if endEQTLDist not in distances:
					#	distances[endEQTLDist] = 0
					#distances[endEQTLDist] += 1
					if endEQTLDist < 200000:
						distancesList.append(endEQTLDist)
		
		#histogram of the distances
		distancesList = np.array(distancesList)
		#plt.hist(np.log(distancesList))
		plt.hist(distancesList)
		plt.show()
		exit()
		
		#Define the sigmoid for beta for eQTLs specifically (for now that is just the counts, we can work from there).
		scaledEQTLCounts = self.computeEQTLBetaSigmoid(cancerTypeSVs["genes"])
		
		
		#For every gene, get how many eQTLs it has.
		#Then determine the beta based on the number of eQTLs
		#Finally, determine alpha by counting the number of SVs overlapping the eQTL layer (what about the ones also overlapping the gene? Should we exclude those for now? )
		allBetas = []
		allAlphas = []
		nonZeroScore = 0
		noEQTLs = 0
		zeroBetaSigmoid = 0
		geneScores = np.zeros([len(geneMap),4])
		for geneInd in range(0, len(cancerTypeSVs["genes"])):
			gene = list(cancerTypeSVs["genes"].keys())[geneInd]
			eQTLCount = scaledEQTLCounts[geneInd] #use the counts scaled to a sigmoid range of -10, 10
			
			beta = 0
			alpha = 0
			
			if len(gene.eQTLs) > 0:
				beta = self.sigmoid(eQTLCount)
				if beta == 0:
					zeroBetaSigmoid += 1
				
				allBetas.append(beta)
			else:
				noEQTLs += 1
				#allBetas.append(beta)
				
				continue #for testing
			
			#Determine alpha
			#Only increase alpha if the variant has not been seen before, to ensure that we do not rank large SVs affecting many eQTLs highest
			seenVariants = dict()
			for eQTL in gene.eQTLs:
				if eQTL in cancerTypeSVs["eQTLs"]: #only the eQTLs that have an SV overlapping them
					
					variants = cancerTypeSVs["eQTLs"][eQTL]
					for variant in variants:
						variantUniqueId = str(variant[0]) + "_" + str(variant[1]) + "_" + str(variant[2]) + "_" + str(variant[3]) + "_" + str(variant[4]) + "_" + str(variant[5])
						if variantUniqueId in seenVariants or variantUniqueId in cancerTypeSVs["genes"][gene]: #also skip the variant if it overlaps the gene itself as well. 
							continue
						seenVariants[variantUniqueId] = 1
						
						alpha += 1 #only increase alpha when we have not seen the variant before. We count the number of variants in the eQTL layer, not the number of eQTLs with at least 1 SV
			
			
			geneScore = alpha * beta
			if geneScore > 0:
				nonZeroScore += 1
			geneInd = geneMap[gene]
			geneScores[geneInd] = [geneScore, alpha, beta, len(gene.eQTLs)] #Also add the number of eQTLs
			
			print("alpha: ", alpha)
			print("beta: ", beta)
			print("score: ", geneScore)
		
		# plt.hist(allBetas)
		# plt.show()
		# exit()
		
		
		#Return a multiplied score for each gene (across all patients). 
		
		
		return geneScores
		