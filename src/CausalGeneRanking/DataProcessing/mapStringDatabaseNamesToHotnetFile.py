"""
	Make hierarchical HotNet input file, and make sure that we map the ensembl ID names from the STRING database file to HG gene names


"""

import sys
import numpy as np

stringFile = sys.argv[1]
geneMapping = sys.argv[2]
hotnetOutFileFolder = sys.argv[3]

geneNameMapping = dict()

with open(geneMapping, 'r') as inF:
	
	lineCount = 0
	
	for line in inF:
		if lineCount < 1:
			lineCount += 1
			continue
		
		line = line.strip()
		splitLine = line.split("\t")
		
		ensId = splitLine[2]
		hgId = splitLine[1]
		if ensId not in geneNameMapping:
			geneNameMapping[ensId] = hgId
			
interactions = []
geneIndex = dict()
indexCount = 1
genes = dict()
with open(stringFile, 'r') as stringF:
	lineCount = 0
	for line in stringF:
		if lineCount < 1:
			lineCount += 1
			continue
		
		splitLine = line.split(" ")
		
		genes[splitLine[0]] = 0
		genes[splitLine[1]] = 0
		
		gene1 = splitLine[0]
		#splitGene1 = gene1Full.split(".")
		#gene1 = splitGene1[1]
		#gene1 = gene1.replace("P", "G") #STRING uses the protein ID instead of the gene ID, difference is the G or P
		gene2 = splitLine[1]
		#splitGene2 = gene2Full.split(".")
		#gene2 = splitGene2[1]
		#gene2 = gene2.replace("P", "G")

		if gene1 not in geneNameMapping or gene2 not in geneNameMapping:
			continue
		
		if gene1 not in geneIndex:
			geneIndex[geneNameMapping[gene1]] = indexCount
			indexCount += 1
		if gene2 not in geneIndex:
			geneIndex[geneNameMapping[gene2]] = indexCount
			indexCount += 1

		gene1Name = geneNameMapping[gene1]
		gene2Name = geneNameMapping[gene2]
		
		interactions.append([gene1Name, gene2Name])



interactions = np.array(interactions)

#For HotNet, we need 2 input files based on the network: 1 with all the gene names and their identifier (numeric), and a file with interactions between these numeric IDs

with open(hotnetOutFileFolder + "/network_1_edge_list.tsv", 'w') as edgeListFile:
	for interaction in interactions:
		
		gene1Ind = geneIndex[interaction[0]]
		gene2Ind = geneIndex[interaction[1]]
		edgeListFile.write(str(gene1Ind) + "\t" + str(gene2Ind) + "\n")
		
with open(hotnetOutFileFolder + "/network_1_index_gene.tsv", 'w') as indexGeneFile:	
	for gene in geneIndex:
		
		indexGeneFile.write(str(geneIndex[gene]) + "\t" + gene + "\n")
		