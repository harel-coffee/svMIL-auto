
import numpy as np
import matplotlib.pyplot as plt

#Load the z-score ranks for both cases
tadRanks = np.loadtxt('pValues.txt', dtype='object')
ruleRanks = np.loadtxt('Output/RankedGenes/0/BRCA/nonCoding_geneSVPairs.txt__nonCodingPairsZScores.txt', dtype='object')

tadRanksFiltered = []
ruleRanksFiltered = []

for rank in ruleRanks:
	
	splitPair = rank[0].split('_')
	pair = splitPair[7] + '_' + splitPair[0]
	
	if pair in tadRanks[:,0]:
		zScore = tadRanks[tadRanks[:,0] == pair][0][5]
		tadRanksFiltered.append(float(zScore))
		ruleRanksFiltered.append(float(rank[3]))
		
		
		if float(rank[3]) > 150:
			print(splitPair)
			print(zScore)

exit()
tadRanksFiltered = np.array(tadRanksFiltered)
ruleRanksFiltered = np.array(ruleRanksFiltered)

tadRanksFiltered = tadRanksFiltered[tadRanksFiltered < np.percentile(tadRanksFiltered,95)]
ruleRanksFiltered = ruleRanksFiltered[ruleRanksFiltered < np.percentile(ruleRanksFiltered,95)]

plt.scatter(np.sort(tadRanksFiltered), np.sort(ruleRanksFiltered))
plt.ylim([-5,180])
plt.xlim([-5,180])
plt.xlabel('Z-scores for TAD-based')
plt.ylabel('Z-scores for rule-based')
plt.show()