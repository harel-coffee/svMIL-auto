"""
	Instead of plotting each TAD individually, plot the pairs of TADs that are disrupted by an SV. Make sure to split this across SV types. 

"""


import sys
import numpy as np
import settings
from inputParser import InputParser
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import matplotlib.pyplot as plt
import glob

bins = 10

def plotExpressionComparison(zScores, allShuffledZScores, tadData, allGenes, filteredSVs):
	
	genes = []
	for gene in allGenes:
		
		genes.append([gene[0], gene[1], gene[2], gene[3].name])
	
	genes = np.array(genes, dtype='object')
	
	#svType = 'DEL'
	svType = sys.argv[1]
	
	
	disruptedZScores = []
	shuffledZScores = []
	
	if svType in ['DEL', 'DUP', 'INV']:
		
		#First find which TAD the gene of the z-score is in.
		
		for pair in zScores:
			
			#get the position of the gene
			
			geneInfo = genes[genes[:,3] == pair[1]][0]
	
			#which TAD is the gene in?
			
			tadChrSubset = tadData[tadData[:,0] == geneInfo[0]]

			matchingTads = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])]
			
			if len(matchingTads) < 1:
					continue
	
			matchingTad = matchingTads[0]

			#Then see which SV disrupts this TAD. 
			
			#first check intrachromosomal SVs
			
			svChr1Subset = filteredSVs[(filteredSVs[:,0] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
			
			
			overlappingSVs = svChr1Subset[(svChr1Subset[:,1] <= matchingTad[2]) * (svChr1Subset[:,2] >= matchingTad[1])]
			
			foundSV = False #make sure that this is added only once if there are multiple SVs
			
			#if there is AN SV overlapping the TAD, add the current z-score.
			for sv in overlappingSVs:
				if sv[8].svType == svType:
					foundSV = True
			if foundSV == True:
				disruptedZScores.append(float(pair[2]))
				
			
	
		
		#repeat but then for the shuffled case.
		for shuffledCase in allShuffledZScores:
			
			for pair in shuffledCase:
				
				#get the position of the gene
				
				geneInfo = genes[genes[:,3] == pair[1]][0]
		
				#which TAD is the gene in?
				
				tadChrSubset = tadData[tadData[:,0] == geneInfo[0]]
				
				if len(tadChrSubset) < 1:
					continue
				
				matchingTads = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])]
	
				if len(matchingTads) < 1:
					continue
	
				matchingTad = matchingTads[0]
	
				#Then see which SV disrupts this TAD. 
				
				#first check intrachromosomal SVs
				
				svChr1Subset = filteredSVs[(filteredSVs[:,0] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
				
				overlappingSVs = svChr1Subset[(svChr1Subset[:,1] <= matchingTad[2]) * (svChr1Subset[:,2] >= matchingTad[1])]
				
				foundSV = False
				
				#if there is AN SV overlapping the TAD, add the current z-score.
				for sv in overlappingSVs:
					
					if sv[8].svType == svType:
						foundSV = True
				if foundSV == True:		
					shuffledZScores.append(float(pair[2]))
		
	elif svType == 'ITX':
		
		#find which SV is in the TAD
		for pair in zScores:
			
			#get the position of the gene
			
			geneInfo = genes[genes[:,3] == pair[1]][0]
	
			#which TAD is the gene in?
			
			tadChrSubset = tadData[tadData[:,0] == geneInfo[0]]

			matchingTads = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])]
			
			if len(matchingTads) < 1:
					continue
	
			matchingTad = matchingTads[0]

			#Then see which SV disrupts this TAD. 
			
			#first check which SVs have chr1 matches			
			svChr1Subset = filteredSVs[(filteredSVs[:,0] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
			
			overlappingSVs = svChr1Subset[(svChr1Subset[:,1] >= matchingTad[1]) * (svChr1Subset[:,1] <= matchingTad[2])]
			
			foundSV = False #make sure that each pair is added only once, even if there are more SVs
			
			for sv in overlappingSVs:
				if sv[8].svType == svType:
					foundSV = True
		
			#repeat for chr2
			svChr2Subset = filteredSVs[(filteredSVs[:,3] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
			
			overlappingSVs = svChr2Subset[(svChr2Subset[:,5] >= matchingTad[1]) * (svChr2Subset[:,5] <= matchingTad[2])]
			
			for sv in overlappingSVs:
				if sv[8].svType == svType:
					foundSV = True
			if foundSV == True:
				disruptedZScores.append(float(pair[2]))


		#repeat for the shuffled case
		for shuffledCase in allShuffledZScores:
			
			for pair in shuffledCase:
				
				geneInfo = genes[genes[:,3] == pair[1]][0]
	
				#which TAD is the gene in?
				
				tadChrSubset = tadData[tadData[:,0] == geneInfo[0]]
	
				matchingTads = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])]
	
				if len(matchingTads) < 1:
						continue
		
				matchingTad = matchingTads[0]
				#Then see which SV disrupts this TAD. 
				
				#first check which SVs have chr1 matches			
				svChr1Subset = filteredSVs[(filteredSVs[:,0] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
				
				overlappingSVs = svChr1Subset[(svChr1Subset[:,1] >= matchingTad[1]) * (svChr1Subset[:,1] <= matchingTad[2])]
				
				foundSV = False #make sure that each pair is added only once, even if there are more SVs
				
				for sv in overlappingSVs:
					if sv[8].svType == svType:
						foundSV = True
			
				#repeat for chr2
				svChr2Subset = filteredSVs[(filteredSVs[:,3] == geneInfo[0]) * (filteredSVs[:,7] == pair[0])]
				
				overlappingSVs = svChr2Subset[(svChr2Subset[:,5] >= matchingTad[1]) * (svChr2Subset[:,5] <= matchingTad[2])]
				
				for sv in overlappingSVs:
					if sv[8].svType == svType:
						foundSV = True
				if foundSV == True:
					shuffledZScores.append(float(pair[2]))
	
	else:
		
		#if we do not check for SVs, just plot the z-scores normally. 
		
		#make a boxplot of the normal z-scores compared to the shuffled z-scores.
		
		for zScore in zScores:
			disruptedZScores.append(float(zScore[2]))

		for shuffledCase in allShuffledZScores:
			for zScore in shuffledCase:
	
				shuffledZScores.append(float(zScore[2]))

	
	#to do the log transform, add a constant to remove negative values.
	#this is the smallest value in the shuffled and non-shuffled set.
	
	# minValue = 0
	# if np.min(disruptedZScores) < np.min(shuffledZScores):
	# 	minValue = np.abs(np.min(disruptedZScores))
	# else:
	# 	minValue = np.abs(np.min(shuffledZScores))
	# 	
	# constant = minValue + 0.0001
	# disruptedZScores = np.array(disruptedZScores)
	# shuffledZScores = np.array(shuffledZScores)
	# 
	# 
	# 
	# disruptedZScores += constant
	# shuffledZScores += constant
	# 
	# disruptedZScores = np.log(disruptedZScores)
	# shuffledZScores = np.log(shuffledZScores)
	
	disruptedZScores = np.array(disruptedZScores)
	shuffledZScores = np.array(shuffledZScores)
	
	
	filteredDisruptedZScores = disruptedZScores[disruptedZScores < np.percentile(disruptedZScores,95)]
	filteredShuffledZScores = shuffledZScores[shuffledZScores < np.percentile(shuffledZScores,95)]
	
	print('plotting')
	
	allData = [filteredDisruptedZScores, filteredShuffledZScores]
	
	plt.boxplot(allData)
	plt.ylim([-10,30])
	plt.savefig('zScores_' + svType + '.svg')
	#plt.show()
	plt.clf()
	allData = [disruptedZScores, shuffledZScores]
	
	plt.boxplot(allData)
	plt.savefig('zScores_' + svType + '_unfiltered.svg')
	#plt.show()

	#do a quick t-test
	z = (np.mean(filteredDisruptedZScores) - np.mean(filteredShuffledZScores)) / float(np.std(filteredShuffledZScores))
	pValue = stats.norm.sf(abs(z))*2
	
	print(pValue)

	return 0
	

def getBinScores(zScores, randomPValuesDir):
	
	splitZScores = []
	for zScore in zScores:
		splitScore = zScore[0].split("_")

		#if zScore[3] != 'False': #only look at the ones which we actually mapped. 
		splitZScores.append([splitScore[0], splitScore[1], float(zScore[5])])
	
	zScores = np.array(splitZScores, dtype='object')
	
	#filteredZScores = zScores[zScores[:,2] < np.percentile(zScores[:,2],95)]
	#zScores = filteredZScores

	#also load all the zScores of the random shuffles.
	shuffledZScoreFiles = glob.glob(randomPValuesDir + '/pValues*')
	
	allShuffledZScores = []
	shuffleCount = 0
	genePatientShuffledZScores = dict()
	for shuffledFile in shuffledZScoreFiles:
		#if shuffleCount > 0: #check only 1 shuffle for now. 
		#	continue

		#if shuffledFile != 'tadDisr/pValues_shuffled_1.txt':
		#	continue
		
		print(shuffledFile)

		shuffledZScores = np.loadtxt(shuffledFile, dtype='object')
		splitShuffledZScores = []
		for zScore in shuffledZScores:
			splitScore = zScore[0].split("_")
			
			if splitScore[0] not in genePatientShuffledZScores:
				genePatientShuffledZScores[splitScore[0]] = dict()
			if splitScore[1] not in genePatientShuffledZScores[splitScore[0]]:
				genePatientShuffledZScores[splitScore[0]][splitScore[1]] = []
			
			#if zScore[3] != 'False':
			
			splitShuffledZScores.append([splitScore[0], splitScore[1], float(zScore[5])])
			genePatientShuffledZScores[splitScore[0]][splitScore[1]].append(float(zScore[5]))
			
		splitShuffledZScores = np.array(splitShuffledZScores, dtype='object')
		
		#filteredZScores = splitShuffledZScores[splitShuffledZScores[:,2] < np.percentile(splitShuffledZScores[:,2],95)]
		#splitShuffledZScores = filteredZScores

		allShuffledZScores.append(splitShuffledZScores)
		shuffleCount += 1
	
	allShuffledZScores = np.array(allShuffledZScores, dtype='object')

	#compute z-score of z-scores
	zScoresOfZScores = []
	for zScore in zScores:
		
		if zScore[0] not in genePatientShuffledZScores:
			continue
		if zScore[1] not in genePatientShuffledZScores[zScore[0]]:
			continue

		negScores = genePatientShuffledZScores[zScore[0]][zScore[1]]
		z = (zScore[2] - np.mean(negScores)) / np.std(negScores)

		zScoresOfZScores.append([zScore[0], zScore[1], z])

	zScoresOfZScores = np.array(zScoresOfZScores, dtype='object')
	print(zScores.shape)
	print(zScoresOfZScores.shape)
	zScores = zScoresOfZScores

	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
			

	#then go through the TADs that are disrupted by a non-coding SV. 
	
	#Get all SVs
	svDir = settings.files['svDir']
	svData = InputParser().getSVsFromFile_hmf(svDir)
	
	#Filter out the coding effect SVs, we want to focus on non-coding SVs. 
	excludedSVs = np.loadtxt(settings.files['excludedSVs'], dtype='object')
	
	svType = 'ITX'
	
	filteredSVs = []
	types = []
	for sv in svData:
	
		#if sv[8].svType != svType:
		#	continue
		
		svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName
		if svEntry not in excludedSVs:
			filteredSVs.append(sv)
			
			if sv[8].svType not in types:
				types.append(sv[8].svType)
	
	print(types)
	filteredSVs = np.array(filteredSVs, dtype='object')
	
	np.save('filteredSVs_tadPairs.npy', filteredSVs)
	filteredSVs = np.load('filteredSVs_tadPairs.npy', allow_pickle=True, encoding='latin1')
	
	#For each SV, determine which TAD it starts and ends in.
	#Keep this as a TAD pair.
	tadFile = settings.files['tadFile']
	tadData = InputParser().getTADsFromFile(tadFile)
	
	#plot the expression between disrupted and non-disrupted TADs
	#plotExpressionComparison(zScores, allShuffledZScores, tadData, allGenes, filteredSVs)
	#exit()
	
	
	tadPairs = dict() #keep the pair as name, and the patients as value. 
	for sv in filteredSVs:
		
		#get the left and rightmost TAD.
		
		#if intrachromosomal, check overlap
		if sv[0] == sv[3]:
			tadChrSubsetInd = sv[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			
			#If the SV start is before the end of the TAD, and the SV end after the start of the TAD, the TAD is overlapped.
			startMatches = sv[1] <= tadChrSubset[:,2]
			endMatches = sv[5] >= tadChrSubset[:,1]
			
			tadMatches = tadChrSubset[startMatches * endMatches]
			
			if tadMatches.shape[0] < 2: #no matches, or overlapping just 1 TAD. 
				continue
			
			#Get the leftmost and rightmost TADs
			farLeftTad = tadMatches[0] #This list is sorted
			farRightTad = tadMatches[tadMatches.shape[0]-1]
		
			tadPair = farLeftTad[0] + '_' + str(farLeftTad[1]) + '_' + str(farLeftTad[2]) + '_' + farRightTad[0] + '_' + str(farRightTad[1]) + '_' + str(farRightTad[2])
			
			if tadPair not in tadPairs:
				tadPairs[tadPair] = []
			tadPairs[tadPair].append(sv[7])
		
		else: #if interchromosomal, determine the TAD based on breakpoints on either chromosome.

			tadChr1SubsetInd = sv[0] == tadData[:,0]
			tadChr1Subset = tadData[tadChr1SubsetInd]
			
			#If the SV start is before the end of the TAD, and the SV end after the start of the TAD, the TAD is overlapped.
			startMatches = sv[1] <= tadChr1Subset[:,2]
			endMatches = sv[5] >= tadChr1Subset[:,1]
			
			tadMatches = tadChr1Subset[startMatches * endMatches]
			
			if tadMatches.shape[0] < 1: #no matches
				continue
			
			#Get the leftmost and rightmost TADs
			farLeftTad = tadMatches[0] #This list is sorted
			
			#repeat for right TAD
			tadChr2SubsetInd = sv[0] == tadData[:,0]
			tadChr2Subset = tadData[tadChr2SubsetInd]
			
			#If the SV start is before the end of the TAD, and the SV end after the start of the TAD, the TAD is overlapped.
			startMatches = sv[1] <= tadChr2Subset[:,2]
			endMatches = sv[5] >= tadChr2Subset[:,1]
			
			tadMatches = tadChr2Subset[startMatches * endMatches]
			
			if tadMatches.shape[0] < 1: #no matches
				continue
			
			farRightTad = tadMatches[0]
			
			tadPair = farLeftTad[0] + '_' + str(farLeftTad[1]) + '_' + str(farLeftTad[2]) + '_' + farRightTad[0] + '_' + str(farRightTad[1]) + '_' + str(farRightTad[2])
			
			if tadPair not in tadPairs:
				tadPairs[tadPair] = []
			tadPairs[tadPair].append(sv[7])
			
	
	#have an additional filter here for the TADs; if there is one TAD pair where we also see the same TAD boundary disrupted again in the same patient, but on another side, we should ignore it for now. 
	
	#if the start of the left TAD is also the end of another pair, or te end of the right TAD is the start of another pair, then we should remove this pair.
	#how often does this happen?
	splitPairs = []
	for pair in tadPairs:
		splitPair = pair.split('_')
		splitPairs.append([splitPair[0], int(splitPair[1]), int(splitPair[2]), splitPair[3], int(splitPair[4]), int(splitPair[5])])
		
	splitPairs = np.array(splitPairs, dtype='object')
	
	
	tadPairsFiltered = dict()
	for pair in splitPairs:
		
		pairChrSubset = splitPairs[splitPairs[:,3] == pair[0]]
		pairStr = '_'.join([str(i) for i in pair])
	
		pairPatients = tadPairs[pairStr]
		
		matched = False
		
		if pair[1] in pairChrSubset[:,5]:
			matchingPairs = pairChrSubset[pairChrSubset[:,5] == pair[1]]
			
			#for these matches, check if they are also disrupted in the same patient.
			for matchedPair in matchingPairs:
				matchedPairStr = '_'.join([str(i) for i in matchedPair])
				
				matchedPairPatients = tadPairs[matchedPairStr]
				
				for patient in matchedPairPatients:
					if patient in pairPatients:
						#print(pair, ' has match in : ', matchedPairStr, ' patient: ', patient)
						matched = True
						
		if pair[5] in pairChrSubset[:,1]:
			matchingPairs = pairChrSubset[pairChrSubset[:,1] == pair[5]]
			#for these matches, check if they are also disrupted in the same patient.
			for matchedPair in matchingPairs:
				matchedPairStr = '_'.join([str(i) for i in matchedPair])
				
				matchedPairPatients = tadPairs[matchedPairStr]
				
				for patient in matchedPairPatients:
					if patient in pairPatients:
						matched = True
						#print(pairStr, ' has match in : ', matchedPairStr, ' patient: ', patient)
		
		#also check if there is not another affected TAD within 1 MB
		
		#if for the first part of the pair there is something overlapping within - 1MB, or for the second part + 1 MB, exclude this pair.
		# overlapSize = 1000000
		# #is there a TAD that starts or ends within the 1 MB?
		# startOverlap = pairChrSubset[(int(pair[1]) - overlapSize >= pairChrSubset[:,4]) + (int(pair[1]) <= pairChrSubset[:,5])]
		# endOverlap = pairChrSubset[(int(pair[5]) + overlapSize <= pairChrSubset[:,1]) + (int(pair[5]) >= pairChrSubset[:,2])]
		# 
		windowOverlap = False
		# if len(startOverlap) > 0 or len(endOverlap) > 0:
		# 	print(pairStr, ' has overlap: ')
		# 	print(startOverlap)
		# 	print(endOverlap)
		# 	windowOverlap = True
			
		if matched == False and windowOverlap == False:
			if pairStr not in tadPairsFiltered:
				tadPairsFiltered[pairStr] = pairPatients
	
	#Collect all random iterations here as well. Get these at once for easy access. For each gene/patient combination, list the zScores. 
	
						
	bins = 10 #have 10 on each side. 
	binZScores = dict()
	for binInd in range(0, bins*2):
			
		if binInd not in binZScores:
			binZScores[binInd] = []
	
	randomBinZScores = dict()
	for binInd in range(0, bins*2):
			
		if binInd not in randomBinZScores:
			randomBinZScores[binInd] = np.zeros([len(allShuffledZScores)])

	perTadPositivePatients = dict()

	for tad in tadPairs:
		perTadPositivePatients[tad] = []

		splitTad = tad.split('_')
		
		#Make a mapping for positions to the right bin.
		
		#determine the size and how large each bin should be
		binSizeTad1 = (float(splitTad[2]) - float(splitTad[1])) / bins

		#binSizeTad1 = (float(splitTad[2]) - (float(splitTad[1]) - offset)) / bins
		
		currentStart = float(splitTad[1]) #start at the TAD start
		#currentStart = float(splitTad[1]) - offset
		binStartsTad1 = [currentStart] #list at which position each bin should start.
		for binInd in range(0, bins):
			
			currentStart += binSizeTad1
			binStartsTad1.append(currentStart)
	
		#repeat for TAD 2
		binSizeTad2 = (float(splitTad[5]) - float(splitTad[4])) / bins
		#binSizeTad2 = ((float(splitTad[5]) + offset) - float(splitTad[4])) / bins
		currentStart = float(splitTad[4]) #start at the TAD start
		binStartsTad2 = [currentStart] #list at which position each bin should start.
		for binInd in range(0, bins):
			
			currentStart += binSizeTad2
			binStartsTad2.append(currentStart)
	
		#Go through the genes of the first TAD; find the genes that will be in this bin
		geneChrSubset = allGenes[allGenes[:,0] == splitTad[0]]
		
		for binInd in range(0, len(binStartsTad1)-1):
			
			
			#get the genes in this bin
			genes = geneChrSubset[(geneChrSubset[:,2] >= binStartsTad1[binInd]) * (geneChrSubset[:,1] <= binStartsTad1[binInd+1])]
			
			#get the z-scores of these genes
			allGeneZScores = []
			for gene in genes:
				geneName = gene[3].name
				
				if geneName in zScores[:,1]: #only add the gene if it has a match.
					
					
					geneZScores = zScores[zScores[:,1] == geneName]
					
					#keep the z-scores separate for each patient
					for patient in range(0, len(geneZScores[:,0])):
						
						if geneZScores[patient,0] not in perTadPositivePatients[tad]:
							perTadPositivePatients[tad].append(geneZScores[patient,0])
						allGeneZScores.append(float(geneZScores[patient,2]))

				
					
				#repeat for all shuffled p-values.
				
				#I need to know for this TAD how many DEGs there are. So that should be per iteration.
				# iterationInd = 0
				# for iterationZScores in allShuffledZScores:
				# 	shuffledGeneZScores = iterationZScores[iterationZScores[:,1] == geneName]
				# 	
				# 	for patient in range(0, len(shuffledGeneZScores[:,0])):
				# 		randomBinZScores[binInd][iterationInd] = randomBinZScores[binInd][iterationInd] + 1
				# 			
				# 	iterationInd += 1
			if len(allGeneZScores) > 0:
				binZScores[binInd].append(np.median(allGeneZScores))
	
		#now for TAD 2, start from where the TAD 1 indices left off. 
		geneChrSubset = allGenes[allGenes[:,0] == splitTad[3]]
		
		for binInd in range(0, len(binStartsTad2)-1):
			
			#get the genes in this bin
			genes = geneChrSubset[(geneChrSubset[:,2] >= binStartsTad2[binInd]) * (geneChrSubset[:,1] <= binStartsTad2[binInd+1])]
			
			#get the z-scores of these genes
			allGeneZScores = []
			for gene in genes:
				geneName = gene[3].name
				
				if geneName in zScores[:,1]:
					
					
					geneZScores = zScores[zScores[:,1] == geneName]
					
					#keep the z-scores separate for each patient
					for patient in range(0, len(geneZScores[:,0])):
						
						if geneZScores[patient,0] not in perTadPositivePatients[tad]:
							perTadPositivePatients[tad].append(geneZScores[patient,0])
						allGeneZScores.append(float(geneZScores[patient,2]))
						
			
							
				#I need to know for this TAD how many DEGs there are. So that should be per iteration.
				# iterationInd = 0
				# for iterationZScores in allShuffledZScores:
				# 	
				# 	shuffledGeneZScores = iterationZScores[iterationZScores[:,1] == geneName]
				# 	for patient in range(0, len(shuffledGeneZScores[:,0])):
				# 
				# 		randomBinZScores[binInd+bins][iterationInd] = randomBinZScores[binInd+bins][iterationInd] + 1
				# 	iterationInd += 1
			if len(allGeneZScores) > 0:
				binZScores[binInd+bins].append(np.median(allGeneZScores))
	#now, we need the additional bins around the sides of the TADs. Where does the effect end?
	#because we need to know what the affected patients are for this TAD, we need an extra loop here, where we use all patients from the previous loop to construct a positive and negative set. 
	
	#also use a map for the gene names
	geneNameConversionMap = dict()
	geneNameConversionFile = sys.argv[2]
	with open(geneNameConversionFile, 'r') as inF:
		
		lineCount = 0
		for line in inF:
			
			if lineCount < 1:
				lineCount += 1
				continue
			line = line.strip()
			splitLine = line.split("\t")
			ensgId = splitLine[3]
			splitEnsgId = ensgId.split('.') #we only keep everything before the dot
			geneName = splitLine[4]
			geneNameConversionMap[splitEnsgId[0]] = geneName
	
	#get the expression data
	expressionFile = sys.argv[3]

	expressionData = []
	samples = []
	with open(expressionFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip()
			if lineCount == 0:
				samples = ['']
				samples += line.split("\t")
	
				lineCount += 1
				continue
			splitLine = line.split("\t")
			fullGeneName = splitLine[0]
			if fullGeneName not in geneNameConversionMap:
				continue
			geneName = geneNameConversionMap[fullGeneName] #get the gene name rather than the ENSG ID
	
			data = splitLine[1:len(splitLine)] 
			fixedData = [geneName]
			fixedData += data
			expressionData.append(fixedData)
	
	expressionData = np.array(expressionData, dtype="object")
	
	#divide the region into 3 bins on each side.
	#so, get the coordinates on each side depending on where the TAD pair starts and ends
	#determine which genes are in these regions
	#add the additional bins.
	
	offset = 500000
	
	leftSideOffset = (float(splitTad[1]) - (float(splitTad[1]) - offset)) / bins

	currentStart = float(splitTad[1]) - offset #start at the TAD start
	binStartsLeft = [currentStart] #list at which position each bin should start.
	for binInd in range(0, 3):
		
		currentStart += leftSideOffset
		binStartsLeft.append(currentStart)
	
	#repeat for the right side
	rightSideOffset = ((float(splitTad[5]) + offset) - float(splitTad[5])) / bins

	currentStart = float(splitTad[5])
	binStartsRight = [currentStart] #list at which position each bin should start.
	for binInd in range(0, 3):
		
		currentStart += rightSideOffset
		binStartsRight.append(currentStart)
	
	#re-make the bin array so that we don't lose the old one in case this does not work.
	#append the additional positions
	
	binZScoresOffset = dict()
	for binInd in range(0, 26):
			
		if binInd not in binZScoresOffset:
			binZScoresOffset[binInd] = []
	
	for binInd in range(0, bins*2):
		binZScoresOffset[binInd+3] = binZScores[binInd]
		
	for tad in tadPairs:		

		
		for binInd in range(0, len(binStartsLeft)-1):
			
			#get the genes in this bin
			genes = geneChrSubset[(geneChrSubset[:,2] >= binStartsLeft[binInd]) * (geneChrSubset[:,1] <= binStartsLeft[binInd+1])]
			
			#get the z-scores of these genes
			
			allGeneZScores = []
			for gene in genes:
				geneName = gene[3].name
				
				#get the expression of this gene in the negative set
				negativeExpr = []
				positiveExpr = []
				
				if geneName not in expressionData[:,0]:
					continue
				
				geneExpr = expressionData[expressionData[:,0] == geneName][0]
				for sample in range(0, len(samples)):
					
					if samples[sample] == '':
						continue
					
					if samples[sample] not in perTadPositivePatients[tad]:
						negativeExpr.append(float(geneExpr[sample]))
					else:
						positiveExpr.append(float(geneExpr[sample]))
				
				for patient in positiveExpr:
					z = (float(patient) - np.mean(negativeExpr)) / float(np.std(negativeExpr))
					allGeneZScores.append(z)
			if len(allGeneZScores) > 0:
				binZScoresOffset[binInd].append(np.median(allGeneZScores))
				
		for binInd in range(0, len(binStartsRight)-1):
			
			#get the genes in this bin
			genes = geneChrSubset[(geneChrSubset[:,2] >= binStartsRight[binInd]) * (geneChrSubset[:,1] <= binStartsRight[binInd+1])]
			
			#get the z-scores of these genes
			
			allGeneZScores = []
			for gene in genes:
				geneName = gene[3].name
				
				if geneName not in expressionData[:,0]:
					continue
				
				negativeExpr = []
				positiveExpr = []
				geneExpr = expressionData[expressionData[:,0] == geneName][0]
				for sample in range(0, len(samples)):
					
					if samples[sample] == '':
						continue
					
					if samples[sample] not in perTadPositivePatients[tad]:
						negativeExpr.append(float(geneExpr[sample]))
					else:
						positiveExpr.append(float(geneExpr[sample]))
				
				for patient in positiveExpr:
					z = (float(patient) - np.mean(negativeExpr)) / float(np.std(negativeExpr))
					allGeneZScores.append(z)
			if len(allGeneZScores) > 0:
				binZScoresOffset[binInd+23].append(np.median(allGeneZScores))
						
	return binZScoresOffset, randomBinZScores


### based on DEGs

#Get the DEGs and their p-values

#plot the z-scores.
pValues = np.loadtxt('pValues.txt', dtype='object')
#pValues = np.loadtxt('tadDisr/pValues_shuffled_1.txt', dtype='object')
#pValues = np.loadtxt('pValues_shuffled_extra.txt', dtype='object')
randomPValuesDir = 'tadDisr/'
binScores, randomBinScores = getBinScores(pValues, randomPValuesDir)

print('plotting')

plt.figure()
allData = []
totalLen = 0
for binInd in range(0, len(binScores)):
	#print('bin:', binInd)
	#print(len(binScores[binInd]))
	#totalLen += len(binScores[binInd])
	#plt.plot(binInd+1, len(binScores[binInd]), marker='*')
	#plt.plot(binInd, len(randomBinScores[binInd]), marker='o') #this can later be a boxplot of 100 iterations

	allData.append(binScores[binInd])
	#allData.append(randomBinScores[binInd])


print(totalLen)

#plt.xticks(range(1,len(binScores)+1))
plt.boxplot(allData)
#plt.ylim(-5,25)
plt.show()

exit()



