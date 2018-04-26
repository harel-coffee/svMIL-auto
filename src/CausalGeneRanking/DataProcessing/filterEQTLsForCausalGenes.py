"""
	The goal of this script is to read the multi-tissue eQTLs from GTEx and filter it by genes that are known to be causal.
	The ENSEMBL IDs first need to be mapped to gene names, and then linked back to the causal gene file.

"""

import sys
import numpy as np

eQTLFile = sys.argv[1]
geneListFile = sys.argv[2]
causalGeneFile = sys.argv[3]
filteredEQTLFile = sys.argv[4]

#1. Read causal genes

def readCausalGeneFile(causalGeneFile):
		
	cosmicGenes = [] #later change this into numpy array for quick overlap
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
			
			#gene = Gene(geneSymbol, chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
			gene = [geneSymbol, chromosome, int(start), int(end)]
			
			cosmicGenes.append(geneSymbol)
			
	
	#cosmicGenes = np.array(cosmicGeneList, dtype='object')
	
	return cosmicGenes

causalGenes = readCausalGeneFile(causalGeneFile)


#2. Read ensembl gene list and make lookup for causal genes only



def readEnsemblGeneFile(geneListFile, causalGenes):

	ensemblIDLookup = dict()
	
	with open(geneListFile, 'r') as geneFile:
		lineCount = 0
		for line in geneFile:
			line = line.strip()
			if lineCount < 1: #skip header line
				lineCount += 1
				continue
			splitLine = line.split("\t")
			
			geneId = splitLine[5]
			geneSymbol = splitLine[6]
			
			if geneSymbol in causalGenes:
				
				if geneId not in ensemblIDLookup:
					ensemblIDLookup[geneId] = ""
				
				ensemblIDLookup[geneId] = geneSymbol
	
	return ensemblIDLookup
			
ensemblIDLookup = readEnsemblGeneFile(geneListFile, causalGenes)


#3. Read eQTL file and only write the line back to the fitered file if it matches on a causal gene

def filterEQTLs(eQTLFile, ensemblIDLookup, filteredEQTLFile):

	with open(filteredEQTLFile, 'wb') as outFile:
		with open(eQTLFile, 'r') as f:
			lineCount = 0
			for line in f:
				line = line.strip()
				
				if lineCount < 1:
					newHeader = "chromosome\tstart\tend\tgene"
					outFile.write(newHeader)
					outFile.write("\n")
					lineCount += 1
					continue
				
				splitLine = line.split("\t")
				geneInfo = splitLine[0]
				splitGeneInfo = geneInfo.split(",")
				ensemblID = splitGeneInfo[1]
				splitEnsemblID = ensemblID.split(".")
				ensemblGeneID = splitEnsemblID[0]
				
				position = splitGeneInfo[0]
				splitPosition = position.split("_")
				chromosome = splitPosition[0]
				start = splitPosition[1]
				ref = splitPosition[2]
				alt = splitPosition[3]
				
				if len(ref) > len(alt): #deletion, start and end differ
					end = str(int(start) + (len(ref) - 1))
				else:
					end = str(int(start) + 1)
				
				if ensemblGeneID in ensemblIDLookup:
					
					newLine = chromosome + "\t" + start + "\t" + end + "\t" + ensemblIDLookup[ensemblGeneID]
					outFile.write(newLine)
					outFile.write("\n")
				
				
				
filterEQTLs(eQTLFile, ensemblIDLookup, filteredEQTLFile)			
			
		

