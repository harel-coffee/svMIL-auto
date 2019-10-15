"""
	For each somatic DEG/non-DEG pair, get the nearest GL of the same type and similar size. 

"""
import sys
import numpy as np


#Split the DEG/non-DEG pairs
svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

print(svGenePairs.shape)

#Read the germline pairs
germlinePairs = np.loadtxt(sys.argv[3], dtype='object')

germlinePairsArray = []
for pair in germlinePairs:
	
	splitPair = pair[0].split("_")
	germlinePairsArray.append([splitPair[0], splitPair[1], int(splitPair[2]), int(splitPair[3]), splitPair[4], int(splitPair[5]), int(splitPair[6]), splitPair[7], splitPair[8]])

germlinePairsArray = np.array(germlinePairsArray, dtype='object')

svGenePairsArray = [] #make this an array so that all positions are easy to match
usedMatches = dict()
negativePairsFeatures = []
positivePairsFeatures = []
for pair in svGenePairs:
	
	splitPair = pair[0].split("_")
	svSize = np.abs(int(splitPair[6]) - int(splitPair[2]))
	
	#find the germline pair that is most similar to this one.
	if splitPair[8] == 'del': #match deletions
		subset = germlinePairsArray[germlinePairsArray[:,8] == 'deletion']
	if splitPair[8] == 'tandem': #match deletions
		subset = germlinePairsArray[germlinePairsArray[:,8] == 'duplication']	
	if splitPair[8] == 'invers': #match deletions
		subset = germlinePairsArray[germlinePairsArray[:,8] == 'inversion']		
	if splitPair[8] == 'transl': #match deletions
		continue
	#determine best size match by using absolute size difference
	glSizes = np.abs(subset[:,6] - subset[:,2])
	
	similarInd = np.argmin(np.abs(glSizes - svSize))
	bestGlMatch = subset[similarInd]
	
	matchStr = bestGlMatch[0] + '_' + bestGlMatch[1] + "_" + str(bestGlMatch[2]) + "_" + str(bestGlMatch[3]) + "_" + bestGlMatch[4] + "_" + str(bestGlMatch[5]) + "_" + str(bestGlMatch[6]) + "_" + bestGlMatch[7] + "_" + bestGlMatch[8]
	if matchStr not in usedMatches:
		usedMatches[matchStr] = 0
	usedMatches[matchStr] += 1	
	
	glMatchFeatures = germlinePairs[germlinePairs[:,0] == matchStr][0]
	
	#Set the features for the deg and non-deg pairs
	if pair[0] in degPairs[:,0]:
		features = [pair[0]] + list(glMatchFeatures[1:])
		positivePairsFeatures.append(features)
	else:
		features = [pair[0]] + list(glMatchFeatures[1:])
		negativePairsFeatures.append(features)
		
		
print(len(usedMatches))
total = 0
for match in usedMatches:
	total += usedMatches[match]
	
print(total)


positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print(positivePairsFeatures)
print(negativePairsFeatures)

	
np.savetxt(sys.argv[1] + '_degPairsFeatures_glMatched.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt(sys.argv[1] + '_nonDegPairsFeatures_glMatched.txt', negativePairsFeatures, fmt='%s', delimiter='\t')

