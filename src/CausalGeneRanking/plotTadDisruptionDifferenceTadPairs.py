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
	
	svType = 'ITX'
	
	disruptedZScores = []
	shuffledZScores = []
	
	if svType in ['DEL', 'DUP', 'INV']:
		
		#First find which TAD the gene of the z-score is in.
		
		for pair in zScores:
			
			#get the position of the gene
			
			geneInfo = genes[genes[:,3] == pair[1]][0]
	
			#which TAD is the gene in?
			
			tadChrSubset = tadData[tadData[:,0] == geneInfo[0]]

			matchingTad = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])][0]

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
	
				matchingTad = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])][0]
	
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

			matchingTad = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])][0]

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

			matchingTad = tadChrSubset[(geneInfo[1] <= tadChrSubset[:,2]) * (geneInfo[2] >= tadChrSubset[:,1])][0]

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
	
	allData = [disruptedZScores, shuffledZScores]
	
	plt.boxplot(allData)
	plt.show()

	#do a quick t-test
	z = (np.mean(disruptedZScores) - np.mean(shuffledZScores)) / float(np.std(shuffledZScores))
	pValue = stats.norm.sf(abs(z))*2
	
	print(pValue)

	exit()
	return 0
	

def getBinScores(zScores, randomPValuesDir):
	
	splitZScores = []
	for zScore in zScores:
		splitScore = zScore[0].split("_")
		
		if zScore[7] == 'False': #only look at the ones which we actually mapped. 
			splitZScores.append([splitScore[0], splitScore[1], zScore[5]])
	
	zScores = np.array(splitZScores, dtype='object')

	#also load all the zScores of the random shuffles.
	shuffledZScoreFiles = glob.glob(randomPValuesDir + '/pValues*')
	
	allShuffledZScores = []
	shuffleCount = 0
	for shuffledFile in shuffledZScoreFiles:
		if shuffleCount > 0: #check only 1 shuffle for now. 
			continue
		print(shuffledFile)
		
		shuffledZScores = np.loadtxt(shuffledFile, dtype='object')
		splitShuffledZScores = []
		for zScore in shuffledZScores:
			splitScore = zScore[0].split("_")
			
			if zScore[7] == 'False':
			
				splitShuffledZScores.append([splitScore[0], splitScore[1], zScore[5]])
			
		splitShuffledZScores = np.array(splitShuffledZScores, dtype='object')
	
		allShuffledZScores.append(splitShuffledZScores)
		shuffleCount += 1

	allShuffledZScores = np.array(allShuffledZScores, dtype='object')

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
	plotExpressionComparison(zScores, allShuffledZScores, tadData, allGenes, filteredSVs)
	
	
	
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
		splitPairs.append([splitPair[0], splitPair[1], splitPair[2], splitPair[3], splitPair[4], splitPair[5]])
		
	splitPairs = np.array(splitPairs, dtype='object')
	
	tadPairsFiltered = dict()
	for pair in splitPairs:
		
		pairChrSubset = splitPairs[splitPairs[:,3] == pair[0]]
		pairStr = '_'.join(pair)
	
		pairPatients = tadPairs[pairStr]
		
		matched = False
		
		if pair[1] in pairChrSubset[:,5]:
			matchingPairs = pairChrSubset[pairChrSubset[:,5] == pair[1]]
			
			#for these matches, check if they are also disrupted in the same patient.
			for matchedPair in matchingPairs:
				matchedPairStr = '_'.join(matchedPair)
				
				matchedPairPatients = tadPairs[matchedPairStr]
				
				for patient in matchedPairPatients:
					if patient in pairPatients:
						print(pair, ' has match in : ', matchedPairStr, ' patient: ', patient)
						matched = True
						
		if pair[5] in pairChrSubset[:,1]:
			matchingPairs = pairChrSubset[pairChrSubset[:,1] == pair[5]]
			#for these matches, check if they are also disrupted in the same patient.
			for matchedPair in matchingPairs:
				matchedPairStr = '_'.join(matchedPair)
				
				matchedPairPatients = tadPairs[matchedPairStr]
				
				for patient in matchedPairPatients:
					if patient in pairPatients:
						matched = True
						print(pairStr, ' has match in : ', matchedPairStr, ' patient: ', patient)
			
		if matched == False:
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
	
	for tad in tadPairs:
		
		splitTad = tad.split('_')
		
		#Make a mapping for positions to the right bin.
		
		#determine the size and how large each bin should be
		binSizeTad1 = (float(splitTad[2]) - float(splitTad[1])) / bins
		
		currentStart = float(splitTad[1]) #start at the TAD start
		binStartsTad1 = [currentStart] #list at which position each bin should start.
		for binInd in range(0, bins):
			
			currentStart += binSizeTad1
			binStartsTad1.append(currentStart)
	
		#repeat for TAD 2
		binSizeTad2 = (float(splitTad[5]) - float(splitTad[4])) / bins
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
			
			for gene in genes:
				geneName = gene[3].name
				
				if geneName in zScores[:,1]: #only add the gene if it has a match.
					
					
					geneZScores = zScores[zScores[:,1] == geneName]
					
					#keep the z-scores separate for each patient
					for patient in range(0, len(geneZScores[:,0])):
						if geneZScores[patient,2] != 'nan' and geneZScores[patient,2] != 'inf':
							print('bin: ', binInd, geneName)
							print(geneZScores)

							binZScores[binInd].append(float(geneZScores[patient,2]))
					
				#repeat for all shuffled p-values.
				
				#I need to know for this TAD how many DEGs there are. So that should be per iteration.
				iterationInd = 0
				for iterationZScores in allShuffledZScores:
					shuffledGeneZScores = iterationZScores[iterationZScores[:,1] == geneName]
					
					for patient in range(0, len(shuffledGeneZScores[:,0])):
						if shuffledGeneZScores[patient,2] != 'nan' and shuffledGeneZScores[patient,2] != 'inf':
							randomBinZScores[binInd][iterationInd] = randomBinZScores[binInd][iterationInd] + 1
							
					iterationInd += 1
			
	
		#now for TAD 2, start from where the TAD 1 indices left off. 
		geneChrSubset = allGenes[allGenes[:,0] == splitTad[3]]
		
		for binInd in range(0, len(binStartsTad2)-1):
			
			#get the genes in this bin
			genes = geneChrSubset[(geneChrSubset[:,2] >= binStartsTad2[binInd]) * (geneChrSubset[:,1] <= binStartsTad2[binInd+1])]
			
			#get the z-scores of these genes
			
			for gene in genes:
				geneName = gene[3].name
				
				if geneName in zScores[:,1]:
					
					
					geneZScores = zScores[zScores[:,1] == geneName]
					
					#keep the z-scores separate for each patient
					for patient in range(0, len(geneZScores[:,0])):
						if geneZScores[patient,2] != 'nan' and geneZScores[patient,2] != 'inf':
							binZScores[binInd+bins].append(float(geneZScores[patient,2]))
							
				#I need to know for this TAD how many DEGs there are. So that should be per iteration.
				iterationInd = 0
				for iterationZScores in allShuffledZScores:
					
					shuffledGeneZScores = iterationZScores[iterationZScores[:,1] == geneName]
					for patient in range(0, len(shuffledGeneZScores[:,0])):
						if shuffledGeneZScores[patient,2] != 'nan' and shuffledGeneZScores[patient,2] != 'inf':
							randomBinZScores[binInd+bins][iterationInd] = randomBinZScores[binInd][iterationInd] + 1
					iterationInd += 1
				
						
	
	return binZScores, randomBinZScores


### based on DEGs

#Get the DEGs and their p-values

#plot the z-scores.
pValues = np.loadtxt('pValues_ranks.txt', dtype='object')
randomPValuesDir = 'tadDisr/'
binScores, randomBinScores = getBinScores(pValues, randomPValuesDir)



plt.figure()
allData = []
totalLen = 0
for binInd in range(0, len(binScores)):
	print('bin:', binInd)
	print(len(binScores[binInd]))
	totalLen += len(binScores[binInd])
	plt.plot(binInd+1, len(binScores[binInd]), marker='*')
	#plt.plot(binInd, len(randomBinScores[binInd]), marker='o') #this can later be a boxplot of 100 iterations

	allData.append(randomBinScores[binInd])

print(totalLen)

#plt.xticks(range(1,len(binScores)+1))
plt.boxplot(allData)
plt.show()

exit()


pValues = []
for binInd in range(0, len(binScores)):
	
	z = (np.mean(binScores[binInd]) - np.mean(randomBinScores[binInd])) / float(np.std(randomBinScores[binInd]))
	pValue = stats.norm.sf(abs(z))
	pValues.append(pValue)


print(pValues)
reject, pAdjusted, _, _ = multipletests(pValues, method='bonferroni') #fdr_bh or bonferroni
print(reject)
print(pAdjusted)



allData = []
for binInd in range(0, bins*2):
	allData.append(binScores[binInd])
	
	#Add significance stars to the plot
	if reject[binInd] == True:
		plt.scatter(binInd, 1800, marker='*')

plt.boxplot(allData)
plt.ylim([-100, 1900])
plt.show()
plt.clf()


allData = []
for binInd in range(0, bins*2):
	allData.append(randomBinScores[binInd])

plt.boxplot(allData)
plt.ylim([-100, 1900])
plt.show()

exit()

##### old code based on z-scores


#First collect the genes that have a z-score (filtered for mutation effects) and get their positions within the TAD
zScores = np.loadtxt('zScores.txt', dtype='object')
binScores = getBinScores(zScores)

randomZScores = np.loadtxt('zScores_random_degs.txt', dtype='object')
randomBinScores = getBinScores(randomZScores)

#Compute the p-value between distributions, which is significant?
pValues = []
for binInd in range(0, len(binScores)):
	
	print(np.mean(binScores[binInd]))

	print(np.mean(randomBinScores[binInd]))
	print(np.std(randomBinScores[binInd]))
	
	#print(binScores[binInd])
	#print(randomBinScores[binInd])
	
	z = (np.mean(binScores[binInd]) - np.mean(randomBinScores[binInd])) / float(np.std(randomBinScores[binInd]))
	pValue = stats.norm.sf(abs(z))
	pValues.append(pValue)


print(pValues)
reject, pAdjusted, _, _ = multipletests(pValues, method='bonferroni') #fdr_bh or bonferroni
print(reject)
print(pAdjusted)

#Make a series of boxplots

import matplotlib.pyplot as plt

allData = []
for binInd in range(0, bins*2):
	allData.append(binScores[binInd])

	#Add significance stars to the plot
	if reject[binInd] == True:
		plt.scatter(binInd, 1900, marker='*')


plt.boxplot(allData)
plt.ylim([-250, 1900])
plt.show()
	
	
	


