"""
	Sample random non-cosmic genes to use in prioritization.

"""
from gene import Gene
import numpy as np
from random import randint

#1. Get the causal genes

causalGeneFile = '../../data/Genes/Census_allTue Apr 10 14_56_44 2018.tsv'

cosmicGenes = [] 
with open(causalGeneFile, 'r') as geneFile:
	
	lineCount = 0
	header = []
	for line in geneFile:
		splitLine = line.split("\t")
		#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
		if lineCount < 1:

			header = splitLine
			lineCount += 1
			continue
			
		#Obtain the gene name and gene position
		
		geneSymbolInd = header.index('Gene Symbol')
		genePositionInd = header.index('Genome Location')
		
		geneSymbol = splitLine[geneSymbolInd]
		genePosition = splitLine[genePositionInd]
		
		#Split the gene position into chr, start and end
		
		colonSplitPosition = genePosition.split(":")
		dashSplitPosition = colonSplitPosition[1].split("-")
		
		chromosome = colonSplitPosition[0]
		start = dashSplitPosition[0].replace('"',"") #apparently there are some problems with the data, sometimes there are random quotes in there
		end = dashSplitPosition[1].replace('"', "")
		
		if start == '' or end == '':
			continue
		
		gene = Gene(geneSymbol, chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
		
		cosmicGenes.append([chromosome, int(start), int(end), gene])
		
#Sort the genes
cosmicGenes = np.array(cosmicGenes, dtype='object')

sortedInd = np.lexsort((cosmicGenes[:,1], cosmicGenes[:,0])) #first sort by chromosome and then by start position. 
cosmicGenesSorted = cosmicGenes[sortedInd]

causalGenes = cosmicGenesSorted


nonCausalGeneFile = "../../data/Genes/allGenesAndIdsHg19.txt"

nonCausalGeneList = []
nonCausalGeneNameDict = dict() #dictionary to keep the names of genes that are already in our list and don't need to be sampled again. 

with open(nonCausalGeneFile, 'r') as geneFile:
	
	lineCount = 0
	for line in geneFile:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			lineCount += 1
			continue
		
		#Obtain the name, chromosome and positions of the gene. 
		
		geneID = splitLine[4]
		
		chrom = splitLine[0]
		splitChrom = chrom.split("chr")
		finalChrom = splitChrom[1]

		start = splitLine[1]
		end = splitLine[2]
		
		geneObj = Gene(geneID, finalChrom, int(start), int(end))
		
		#nonCausalGeneList.append([finalChrom, start, end, geneID, geneObj])
		if geneID not in nonCausalGeneNameDict:
			nonCausalGeneList.append([finalChrom, start, end, geneID])
			nonCausalGeneNameDict[geneID] = 0
			
nonCausalGenes = np.array(nonCausalGeneList)

causalGeneNames = dict()
for gene in causalGenes:
	causalGeneNames[gene[3].name] = 0


randomNonCausalGenes = []

#Include our specific set of genes for testing
#(ZFP57 and HLA-F and IFITM4P)
sampledGenes = dict()
for gene in nonCausalGenes:
	
	geneName = gene[3]
	if geneName == "ZFP57" or geneName == "HLA-F" or geneName == "IFITM4P":
		print geneName
		randomNonCausalGenes.append([gene[0], gene[1], gene[2], gene[3]])
		sampledGenes[geneName] = 0



#Then go through the list again, but take a random subset. Also make sure that the genes do not overlap with the COSMIC genes.
noOfGenesToSample = 700 #there are not more than 700, apparently. 
while len(sampledGenes) <= noOfGenesToSample:
#while True:
	
	randomIndex = randint(0, nonCausalGenes.shape[0]-1)
	
	geneName = nonCausalGenes[randomIndex,3]
	if geneName not in sampledGenes:

		#Add the gene to our list if not in cosmic
		
		if geneName not in causalGeneNames:
			
			randomNonCausalGenes.append([nonCausalGenes[randomIndex][0], nonCausalGenes[randomIndex][1], nonCausalGenes[randomIndex][2], nonCausalGenes[randomIndex][3]])
			sampledGenes[geneName] = 0
		
#print nonCausalGeneNameDict
randomNonCausalGenes = np.array(randomNonCausalGenes)


#Write the non causal genes to a file
outFile = "../../data/Genes/nonCosmicSubset.txt"

with open(outFile, 'w') as outF:
	
	#for gene in randomNonCausalGenes:
	for gene in nonCausalGenes:
		
		outF.write(gene[3] + "\t" + gene[0] + "\t" + gene[1] + "\t" + gene[2] + "\n")
		

