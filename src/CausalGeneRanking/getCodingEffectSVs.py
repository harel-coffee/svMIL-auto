import sys
import numpy as np


#1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way

def getSVsWithCodingEffects():
	
	codingPairs = np.loadtxt(sys.argv[1], dtype='object')
	degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')
	cosmicGenesFile = sys.argv[3]
	
	cosmicGenes = []
	with open(cosmicGenesFile, 'r') as f:
		lineCount = 0
		for line in f:
			if lineCount == 0:
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
		
			geneName = splitLine[0]
			cosmicGenes.append(geneName)
	
	print(codingPairs.shape)
	
	#1. Get the indices at which there is a COSMIC gene
	#first get all the genes separately
	genesOnly = []
	for pair in codingPairs:
		splitPair = pair.split("_")
		genesOnly.append(splitPair[0])
	genesOnly = np.array(genesOnly)
	
	#then find which indices are genes that are also in COSMIC.
	cosmicGeneInd = np.in1d(genesOnly,cosmicGenes).nonzero()
	
	#Extract the pairs
	cosmicPairs = codingPairs[cosmicGeneInd]
	print(cosmicPairs.shape)
	
	codingEffectSVs = dict()
	for pair in cosmicPairs:
		splitPair = pair.split("_")
		sv = "_".join(splitPair[1:])
		codingEffectSVs[sv] = 0
	
	#2. Do the same but then with DEG pairs
	for pair in degPairs[:,0]:
		splitPair = pair.split("_")
		sv = "_".join(splitPair[1:])
		codingEffectSVs[sv] = 0
	
	return codingEffectSVs
	
codingEffectSVs = getSVsWithCodingEffects()
print("Number of SVs filtered out with coding effects: ", len(codingEffectSVs))
np.savetxt(sys.argv[1] + '_codingEffectSVs2.txt', list(codingEffectSVs.keys()), delimiter='\t', fmt='%s')