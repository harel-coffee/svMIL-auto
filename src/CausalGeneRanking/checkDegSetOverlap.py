"""
	Get the DEG set from the TAD disruption analysis and the one from linking SVs to genes. Which overlap? Are there things missing?

"""

import sys
import numpy as np

tadPairs = np.loadtxt('pValues_ranks.txt', dtype='object')
rulePairs = np.loadtxt('Output/RankedGenes/0/BRCA/nonCoding_geneSVPairs.txt__nonCodingPairsRanks.txt', dtype='object')

print(tadPairs)

tadPairsAnn = []
for pair in tadPairs:

	if pair[0] == 'CPCT02010447T_RNF5':
		print('tad pair')
		print(pair)
	else:
		continue
	if pair[7] == 'False':
		tadPairsAnn.append(pair[0] + '_' + pair[7])
		print(tadPairsAnn)
exit()
for pair in rulePairs:

	if pair[0] == 'CPCT02010447T_RNF5':
		print('rule pair')
	else:
		continue


	if pair[3] == 'False':
		rulePairAnn = pair[0] + '_' + pair[3]
		#check if the rule pair is also found in the TAD pairs.
		#everything found by rule-pairs should also be in TAD pairs.
		if rulePairAnn not in tadPairsAnn:
			
			print('pair not in tad pairs: ', rulePairAnn)
			print(pair)
	