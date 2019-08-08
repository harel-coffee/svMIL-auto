"""
	For each non-coding DEG gene, check if the gene has coding evidence in other samples in the form of coding SVs and coding SNVs

"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import sys
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import re

degGenes = np.loadtxt(sys.argv[1], dtype="object")
geneRanks = np.loadtxt(sys.argv[2], dtype="object")
snvFolder = sys.argv[3]



#For each gene, check if there are coding samples. This is ME, so everything in the deg set should be non-coding only.
codingEvidence = []
for geneExpr in degGenes:
	
	geneName = geneExpr[0].split("_")[0]
	
	geneRank = geneRanks[np.where(geneRanks[:,0] == geneName)][0]
	nonCodingSamples = geneRank[31].split(",")
	if geneRank[32] != "None":
		codingSamples = geneRank[32].split(",")
		
		codingEvidence.append([geneName, nonCodingSamples, codingSamples])
	else:
		codingEvidence.append([geneName, nonCodingSamples, []])
	
codingEvidence = np.array(codingEvidence, dtype="object")
print(codingEvidence)

print(degGenes.shape)
print(codingEvidence.shape)
	
#make a bar plot showing how many genes have evidence

evidenceCounts = dict()
for gene in codingEvidence:
	
	if len(gene[2]) not in evidenceCounts:
		evidenceCounts[len(gene[2])] = 0
	evidenceCounts[len(gene[2])] += 1	


# plt.bar(evidenceCounts.keys(), evidenceCounts.values())
# plt.show()

#Also collect coding SnV evidence
#Get all relevant samples from the SNV folder

snvFiles = [f for f in listdir(snvFolder) if isfile(join(snvFolder, f))]
snvEvidence = dict()
for geneExpr in degGenes:
	
	geneName = geneExpr[0].split("_")[0]
	sampleName = geneExpr[0].split("_")[1]
	splitSampleName = sampleName.split("brca")[1]
	snv = False 
	for snvFile in snvFiles:
		match = re.search(splitSampleName, snvFile)
		if match != None:
			#Read the file and check if the gene has an SNV
			with open(snvFolder + '/' + snvFile, 'r') as inF:
				for line in inF:
					splitLine = line.split("\t")
					snvFileGene = splitLine[0]
					if snvFileGene == geneName:
						snv = True
	
	if geneName not in snvEvidence:
		snvEvidence[geneName] = []
	
	if snv == True:
		snvEvidence[geneName].append(sampleName)

print(snvEvidence)
exit()
		


#Plot a bar with both the coding and non-coding counts, each gene on the x-axis
geneCodingEvidences = dict()
geneNonCodingEvidences = dict()
count = 0
for gene in codingEvidence:
	
	geneCodingEvidences[count] = len(gene[2])
	geneNonCodingEvidences[count] = len(gene[1])
	
	count += 1

plt.bar(list(geneCodingEvidences.keys()), list(geneCodingEvidences.values()), color='b')
plt.bar(list(geneNonCodingEvidences.keys()), list(geneNonCodingEvidences.values()), bottom = list(geneCodingEvidences.values()), color='g')
plt.xticks(list(geneCodingEvidences.keys()), codingEvidence[:,0], rotation=90)
plt.show()
	

