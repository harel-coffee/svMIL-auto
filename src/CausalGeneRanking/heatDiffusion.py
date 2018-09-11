# Script for testing heat diffusion even when graphs are large

import numpy as np
import pickle as pkl
from numpy import linalg
from numpy.linalg import inv

#1. Read the data from regions.csv

#All relations should be within a chromosome, so we can look per chromosome individually

#It is better to keep the regions sorted in the adjacency matrix. So, if we keep a numpy matrix, we can sort and at the same time use the index as an index for the region.
#In that case we should not use dictionaries, or convert that to a dictionary after sorting. 

regionsFile = "../Parsers/HiC/Intrachromosomal/regions.csv"

#1. Get all regions on a specific chromosome
#Sort these regions by their position on the chromosome

regions = dict()
relations = dict()
with open(regionsFile, 'r') as inFile:
	lineCount = 0
	for line in inFile:
		if lineCount < 1:
			lineCount += 1
			continue
		line = line.strip()
		splitLine = line.split(",")
		
		splitRegion = splitLine[0].split("_")
		chrom = splitRegion[0]
	
		if chrom not in regions:
			regions[chrom] = []
			
		if chrom not in relations:
			relations[chrom] = dict()
			
		regions[chrom].append([splitLine[0], int(splitRegion[1])])
			

#do the sorting per chromosome and store the sorted indices in a dictionary
sortedRegions = dict()
for chrom in regions:
	
	chromRegions = np.array(regions[chrom], dtype=object)
	#sort this by the 2nd column
	
	chromRegions = chromRegions[chromRegions[:,1].argsort()]
	sortedRegions[chrom] = chromRegions
	
#Now make the lookup based on the sorted indices

regionsLookup = dict()
reverseRegionsLookup = dict()
for chrom in sortedRegions:
	if chrom not in regionsLookup:
		regionsLookup[chrom] = dict()
	if chrom not in reverseRegionsLookup:
		reverseRegionsLookup[chrom] = dict()
	
	currentIndex = 0
	for region in range(0, len(sortedRegions[chrom])):
		regionName = sortedRegions[chrom][region,0]
		
		regionsLookup[chrom][regionName] = currentIndex
		reverseRegionsLookup[chrom][currentIndex] = regionName
		currentIndex += 1	

print "done making lookup"

#1.2 Also read the relationships between the regions

regionsRelFile = "../Parsers/HiC/Intrachromosomal/regions_regions_rel.csv"

with open(regionsRelFile, 'r') as inFile:
	lineCount = 0
	for line in inFile:
		if lineCount < 1:
			lineCount += 1
			continue
		line = line.strip()
		splitLine = line.split(",")
		
		splitRegion = splitLine[0].split("_")
		chrom = splitRegion[0]
		
		regionIndex = regionsLookup[chrom][splitLine[0]]
		relationIndex = regionsLookup[chrom][splitLine[1]]
		
		if regionIndex not in relations[chrom]:
			relations[chrom][regionIndex] = []

		relations[chrom][regionIndex].append(relationIndex)

print "done initializing relations"



#2. Store in an adjacency matrix. How large is this matrix and can we store it in memory?
#delete this file before appending to it, each re-run should be a replacement of the file
scoreOutFile = "diffusionScores.txt"
open(scoreOutFile, 'w').close()

for chrom in relations:
	print "chrom is: ", chrom
	

	print "processing: ", chrom
		
	
	noOfRegions = max(relations[chrom]) + 3 #I'm not sure about this number, will need to be fixed later but now works to prevent out of bounds. 
	
	
	#Make submatrices for the regions.
	#Keep the indices to ensure that we can map back to the final big matrix
	
	blockNum = 20 #how many blocks do we want to divide our data in?
	overlap = 0.5 #how much overlap between blocks?
	blocksToExplore = int(blockNum / overlap)
	
	blockSize = int(noOfRegions / float(blockNum))
	currentIndex = 0

	print "max rel: "
	print max(relations[chrom][max(relations[chrom])])

	mat = np.zeros([noOfRegions, noOfRegions])
	
	for regionId, relationIds in relations[chrom].items():
	#	print "rel:", max(relationIds)
		
		for relationId in relationIds:
			#print 'regionId', regionId
			#print 'relId', relationId
			mat[regionId][relationId] = 1
			mat[relationId][regionId] = 1
	
	print mat.shape

	### Block processing 
	
	#Store the final diffusion values in here
	finalF = np.zeros([noOfRegions, noOfRegions], dtype=float)
	
	#Process the adjacency matrix in blocks
	#We can assume that most important interactions will not be too far away
	#Although one issue here is that these again may operate further away, and thus in steps we would still end up at the end of the chromosome.
	#What is the distance between min and max within blocks?
	
	print "processing in blocks of size: ", blockSize
	
	currentIndex = 0
	for block in range(0, blocksToExplore):
		currentMax = currentIndex + blockSize
		currentBlock = mat[currentIndex:currentMax, currentIndex:currentMax]
		
		#This part should ideally be before the block processing, but crashes on chromosome 1, so I will do block processing for this part in block processing as well for all chromosome
		#First compute the weighted density
		weightedAdj = np.zeros([currentBlock.shape[0], currentBlock.shape[0]], dtype=float)
		
		#Loop over all regions in this current block
		for row in range(0, currentBlock.shape[0]):
			region = reverseRegionsLookup[chrom][row]
			
			#rowDegree = len(relations[chrom][row]) #The number of interactions each region has
			rowDegree = sum(currentBlock[row,:])
			
			if rowDegree > 0:
				weightedAdj[row] = currentBlock[row] / float(rowDegree)
				
		#print "weighted matrix"
		#Do computation of 1-B times weightedAdj
		beta = 0.8
		#wow why are there no good ways to call this thing I don't know what to call it
		#print "multiplying with beta"
		weightedAdjSubmatrix = (1-beta) * weightedAdj
		
		#Look at the minimum and maximum region, how far are these apart within a block?
		# minRegion = float("inf")
		# maxRegion = -float("inf")
		# 
		# 
		# for region in range(currentIndex, currentMax):
		# 	
		# 	regionName = reverseRegionsLookup[chrom][region]
		# 	splitRegionName = regionName.split("_")
		# 	regionNum = int(splitRegionName[1])
		# 	
		# 	if regionNum < minRegion:
		# 		minRegion = regionNum
		# 	if regionNum > maxRegion:
		# 		maxRegion = regionNum
		# 	
		# print minRegion
		# print maxRegion
		
		print "current block size: ", currentBlock.shape
		
		#Do the heat diffusion for each block separately
		
		
		#obtain the submatrix of the weighted adjacency matrix
		#weightedAdjSubmatrix = weightedWeightedAdj[currentIndex:currentMax, currentIndex:currentMax]

		#subtract the weighted submatrix from the adjacency matrix
		subtractedMatrices = currentBlock - weightedAdjSubmatrix
		
		print "inverting the matrix"
		#Compute the inverse on the subtracted submatrix
		invertedMatrix = np.linalg.pinv(subtractedMatrices)

		#combine the blocks by taking the highest value if the indices overlap
		
		#do a max to compare the two regions of the matrices (the full matrix and the current submatrix)
		#which values are true should be filled with the current matrix
		
		with open(scoreOutFile, 'a') as scoreFile:
				
			#not sure if there is a faster way
			subFinalF = finalF[currentIndex:currentMax, currentIndex:currentMax]
			for row in range(0, invertedMatrix.shape[0]):
				for col in range(0, invertedMatrix.shape[1]):
					if invertedMatrix[row][col] > subFinalF[row][col]:
						mappedRow = currentIndex + row
						mappedCol = currentIndex + col
						finalF[mappedRow][mappedCol] = beta * invertedMatrix[row][col]
					
				#Get the value of each region	
				#write this value to a file right away to prevent having to do this loop again.
				region = reverseRegionsLookup[chrom][row]
				
				
				scoreFile.write(region + "\t" + str(finalF[row][row]) + "\n")
				#print region
				#print finalF[row][row]
				#exit()
			
		#higherHeatIndices = np.where(invertedMatrix > finalF[currentIndex:currentMax, currentIndex:currentMax])
		
		currentIndex += int(blockSize * overlap) #make sure that blocks overlap

	#clean up properly
	del finalF
	
	
	
		
	#Compute F by weighting with beta after all submatrices have been explored
	#F = beta * finalF
	#print F
	
	#Link these scores back to the corresponding regions and write this to a file.
	
	
	#
	
	#exit() #break after 1 chromosome
	
#After this, we need to write the scores of the network somewhere so that these can be used in the tool.
#
