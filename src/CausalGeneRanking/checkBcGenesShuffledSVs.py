"""
	Every time SVs are shuffled, see if there is indeed an SV in the same TAD as a BC gene. 

"""
from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from inputParser import InputParser
from genomicShuffler import GenomicShuffler

#get all BC genes
breastCancerGenesFile = sys.argv[1]
breastCancerGenes = []
with open(breastCancerGenesFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
		
		breastCancerGenes.append(line)

#get all cosmic genes
cosmicGenesFile = sys.argv[2]
cosmicGenes = []
with open(cosmicGenesFile, 'rb') as f:
	lineCount = 0
	for line in f:
		if lineCount == 0:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")

		geneName = splitLine[0]
		splitPos = splitLine[3].split(":")
		
		if splitPos[1] == '-':
			continue
		splitSplitPos = splitPos[1].split("-")
		chrom = 'chr' + splitPos[0]
		start = int(splitSplitPos[0].replace('"', ''))
		end = int(splitSplitPos[1].replace('"', ''))
		cosmicGenes.append([geneName, chrom, start, end])
cosmicGenes = np.array(cosmicGenes, dtype='object')

#show positions of BC genes
for gene in cosmicGenes:
	if gene[0] in breastCancerGenes:
		print(gene)
		
		
#Get the TADs that these genes are in
tads = InputParser().getTADsFromFile(sys.argv[3])

#Get SVs
codingEffectSVs = np.loadtxt(sys.argv[4], dtype='object')
somaticSVs = InputParser().getSVsFromFile(sys.argv[5], "all", codingEffectSVs)
genomicShuffler = GenomicShuffler()
#somaticSVs = genomicShuffler.shuffleSVs(somaticSVs)

geneSVCounts = dict()
for gene in cosmicGenes:
	if gene[0] in breastCancerGenes:
		geneSVCounts[gene[0]] = 0
		
		tadChrSubset = tads[tads[:,0] == gene[1]]
		
		#get which tad the gene is in: gene start < tad end, gene end > tad start
		matchingTads = tadChrSubset[(gene[2] <= tadChrSubset[:,2]) * (gene[3] >= tadChrSubset[:,1])]
		print(gene)
		print(matchingTads)
		
		#Check which SVs are in this TAD, and of which type these are.
		svsInTad = []
		
		for sv in somaticSVs:
			svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
			
			if sv[0] == sv[3]:
				
				for tad in matchingTads:
					if sv[0] == tad[0]:
						if sv[1] <= tad[2] and sv[5] >= tad[1]:
							svsInTad.append(svStr)
							geneSVCounts[gene[0]] += 1
			if sv[0] != sv[3]:
				
				
				for tad in matchingTads:
					#first check if there is a match on chr1
					if sv[0] == tad[0]:
						if sv[1] <= tad[2] and sv[2] >= tad[1]:
							svsInTad.append(svStr)
							geneSVCounts[gene[0]] += 1
					if sv[3] == tad[0]:
						if sv[4] <= tad[2] and sv[5] >= tad[1]:
							svsInTad.append(svStr)
							geneSVCounts[gene[0]] += 1
	
		print(np.array(svsInTad))
		
print(geneSVCounts)

#If we get all SVs, how often do we find SVs in the same TAD as these genes?
tadTrans = dict()
tadTransCounts = dict()
for sv in somaticSVs:
	svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
	
	if sv[0] != sv[3]:
				
		
		for tad in tads:
			if tad[3] not in tadTrans:
				tadTrans[tad[3]] = []
				tadTransCounts[tad[3]] = 0
			#first check if there is a match on chr1
			if sv[0] == tad[0]:

				if sv[1] <= tad[2] and sv[2] >= tad[1]:
					tadTrans[tad[3]].append(svStr)
					tadTransCounts[tad[3]] += 1
			if sv[3] == tad[0]:
				if sv[4] <= tad[2] and sv[5] >= tad[1]:
					tadTrans[tad[3]].append(svStr)
					tadTransCounts[tad[3]] += 1

print(tadTrans)
print(tadTransCounts)
	
import matplotlib.pyplot as plt

freqs = dict()
for value in tadTransCounts.values():
	if value not in freqs:
		freqs[value] = 0
	freqs[value] += 1

#plt.hist(tadTransCounts.values())
plt.bar(list(freqs.keys()), list(freqs.values()))
plt.show()


		
