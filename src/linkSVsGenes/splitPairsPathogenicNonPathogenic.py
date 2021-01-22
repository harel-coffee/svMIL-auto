
"""
	Script to read the rule-based SV-gene pairs and split these into pathogenic & non-pathogenic pairs based on the z-scores.

	We use this mainly for plotting fig 2e.

"""

import sys
import numpy as np

outDir = sys.argv[1]

#get the pairs identified with the rules, and the z-scores
svGenePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_', dtype='object')
zScores = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object')


#also load the mutation pairs.
#In the current zScore setup, some genes will be affected by an overlapping duplication AND CNV AMP,
#but sometimes we find no evidence that the duplication actually functions in a non-coding way.
#so for those genes, we need to remove the ones that have only CNV AMP and no non-coding duplication,
#because otherwise there would be too many false positives.
mutDir = outDir + '/patientGeneMutationPairs/'
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()

#then also split the pairs so that we can easily check for cases without non-coding dups.
splitSVGenePairs = []
for pair in svGenePairs:
	splitPair = pair[0].split('_')

	splitSVGenePairs.append(splitPair[7] + '_' + splitPair[0] + '_' + splitPair[12])

positivePairsFeatures = []
negativePairsFeatures = []

for pair in svGenePairs:

	splitPair = pair[0].split('_')
	shortPair = splitPair[7] + '_' + splitPair[0]

	if shortPair in zScores[:,0]:

		zPairInfo = zScores[zScores[:,0] == shortPair][0]

		if float(zPairInfo[5]) > 1.5 or float(zPairInfo[5]) < -1.5:
		#if float(zPairInfo[5]) > 2 or float(zPairInfo[5]) < -2:

			#only add to the positive set if there is no CNV amp without duplication.
			if splitPair[0] in cnvPatientsAmp[splitPair[7]] and shortPair + '_DUP' not in splitSVGenePairs:
				negativePairsFeatures.append(pair)
			else:
				positivePairsFeatures.append(pair)
		else:
			negativePairsFeatures.append(pair)


positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print(positivePairsFeatures.shape)
print(negativePairsFeatures.shape)

np.savetxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_nonPathogenicPairsFeatures.txt', negativePairsFeatures, fmt='%s', delimiter='\t')