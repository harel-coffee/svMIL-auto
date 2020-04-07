"""

	show the difference between positive and negative pairs.
	Which features are different?

"""

import numpy as np
import sys

#load pairs
allPairs = np.loadtxt(sys.argv[1], dtype='object')

degPairs = np.loadtxt(sys.argv[2], dtype='object')


svType = 'DEL'

positivePairs = []
negativePairs = []
positivePairLabels = []
negativePairLabels = []
for pair in allPairs:

	splitPair = pair[0].split('_')

	if splitPair[12] != svType:
		continue

	shortPair = splitPair[7] + '_' + splitPair[0]

	if shortPair not in degPairs[:,0]:
		continue

	pairInfo = degPairs[degPairs[:,0] == shortPair][0]

	if float(pairInfo[5]) > 1.5 or float(pairInfo[5]) < -1.5:

		positivePairs.append([float(i) for i in pair[1:]])
		positivePairLabels.append(pair[0])
	else:
		negativePairs.append([float(i) for i in pair[1:]])
		negativePairLabels.append(pair[0])

positivePairs = np.array(positivePairs)
negativePairs = np.array(negativePairs)

print(positivePairs.shape)
print(negativePairs.shape)

#compare features between the groups.

goodInstancesAvg = np.mean(positivePairs[:,0:60], axis=0)
print(goodInstancesAvg)
badInstancesAvg = np.mean(negativePairs[:,0:60], axis=0)

#try simple classification



import matplotlib.pyplot as plt
#
# barWidth = 0.35
# plt.bar(range(0, len(goodInstancesAvg)), goodInstancesAvg, width=barWidth, color='blue')
# #plt.xticks(range(0, len(xlabels)), xlabels, rotation=90)
# #plt.show()
# pos = np.arange(len(goodInstancesAvg))
# r2 = [i + barWidth for i in pos]
# plt.bar(r2, badInstancesAvg, color='orange', width=barWidth)
# #plt.xticks(r2, xlabels, rotation=90)
# plt.show()

## how often do we see biallelic loss in the positive vs negative set?

positiveCounts = dict()
negativeCounts = dict()
for pair in positivePairLabels:

	splitPair = pair.split('_')

	shortPair = splitPair[7] + '_' + splitPair[0]

	if shortPair not in positiveCounts:
		positiveCounts[shortPair] = 0
	positiveCounts[shortPair] += 1

for pair in negativePairLabels:

	splitPair = pair.split('_')

	shortPair = splitPair[7] + '_' + splitPair[0]

	if shortPair not in negativeCounts:
		negativeCounts[shortPair] = 0
	negativeCounts[shortPair] += 1

#compare distributions.
biallelicCountPos = 0
biallelicCountNeg = 0
biallelicDistPos = dict()
biallelicDistNeg = dict()
for pair in positiveCounts:

	if positiveCounts[pair] not in biallelicDistPos:
		biallelicDistPos[positiveCounts[pair]] = 0
	biallelicDistPos[positiveCounts[pair]] += 1
	if positiveCounts[pair] > 1:
		biallelicCountPos += 1

for pair in negativeCounts:

	if negativeCounts[pair] not in biallelicDistNeg:
		biallelicDistNeg[negativeCounts[pair]] = 0
	biallelicDistNeg[negativeCounts[pair]] += 1

	if negativeCounts[pair] > 1:
		biallelicCountNeg += 1

print('biallelic positive: ', biallelicCountPos / len(positivePairLabels))
print('biallelic negative: ', biallelicCountNeg / len(negativePairLabels))

for count in biallelicDistPos:
	print('count: ', count, ': ', biallelicDistPos[count] / len(positivePairLabels))

for count in biallelicDistNeg:
	print('count: ', count, ': ', biallelicDistNeg[count] / len(negativePairLabels))

print(biallelicDistPos)
print(biallelicDistNeg)

