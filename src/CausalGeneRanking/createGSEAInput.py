"""
	Get the rule-based pairs, and make these into a .grp file.
	Format:
	#geneset name
	pair 1
	pair 2
	...
	
	Also make the .rnk file, containing the z-scores of z-scores for the tad-based idetified pairs
	Format:
	pair1	rank
	pair2	rank
	...

"""

import glob
import numpy as np

#1. Load the rule-based pairs
rulePairs = np.loadtxt('Output/RankedGenes/0/BRCA/nonCoding_geneSVPairs.txt__nonCodingPairsZScores.txt', dtype='object')

#2. Make the geneset, using patient ID and gene name only.

#split the names first to get the patient ID and gene name
splitPairs = []
for pair in rulePairs[:,0]:
	
	splitPair = pair.split('_')
	gene = splitPair[0]
	patient = splitPair[7]
	
	newPair = patient + '_' + gene
	splitPairs.append(newPair)
	

np.savetxt('ruleBasedGeneSet.grp', splitPairs, fmt='%s', delimiter='\t')


#3. Compute the rule-based z-scores of z-scores

def getZScoresOfZScores(zScores, shuffledZScoresDir):
	
	splitZScores = []
	for zScore in zScores:
		splitScore = zScore[0].split("_")
		splitZScores.append([splitScore[0], splitScore[1], float(zScore[5])])
	
	zScores = np.array(splitZScores, dtype='object')
	
	#also load all the zScores of the random shuffles.
	shuffledZScoreFiles = glob.glob(shuffledZScoresDir + '/pValues*')
	
	allShuffledZScores = []
	shuffleCount = 0
	genePatientShuffledZScores = dict()
	for shuffledFile in shuffledZScoreFiles:
	
		shuffledZScores = np.loadtxt(shuffledFile, dtype='object')
		splitShuffledZScores = []
		for zScore in shuffledZScores:
			splitScore = zScore[0].split("_")
			
			if splitScore[0] not in genePatientShuffledZScores:
				genePatientShuffledZScores[splitScore[0]] = dict()
			if splitScore[1] not in genePatientShuffledZScores[splitScore[0]]:
				genePatientShuffledZScores[splitScore[0]][splitScore[1]] = []
			
			splitShuffledZScores.append([splitScore[0], splitScore[1], float(zScore[5])])
			genePatientShuffledZScores[splitScore[0]][splitScore[1]].append(float(zScore[5]))
			
		splitShuffledZScores = np.array(splitShuffledZScores, dtype='object')
		
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
		
		if np.std(negScores) == 0:
			continue
		
		z = (zScore[2] - np.mean(negScores)) / np.std(negScores)
		
		

		zScoresOfZScores.append([zScore[0] + '_' + zScore[1], np.abs(z)])

	zScoresOfZScores = np.array(zScoresOfZScores, dtype='object')

	return zScoresOfZScores

zScores = np.loadtxt('pValues.txt', dtype='object')
shuffledZScoreDir = 'tadDisr/'
zScoresOfZScores = getZScoresOfZScores(zScores, shuffledZScoreDir)

#Then output these ranks as well.
np.savetxt('tadBasedRanks.rnk', zScoresOfZScores, fmt='%s', delimiter='\t')




