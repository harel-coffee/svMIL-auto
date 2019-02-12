"""
	The goal of this script is to check for lncRNAs and miRNAs how many eQTLs are in the GTEx database. 

"""

import sys
import re
import numpy as np

gffInFile = sys.argv[1]
eQTLInFile = sys.argv[2]
geneListFile = sys.argv[3]

lncRNAs = dict()
miRNAs = dict()
with open(gffInFile) as inF:
	
	for line in inF:
		
		if re.search("#", line) == None:
			line = line.strip()
			splitLine = line.split("\t")
			
			
			if splitLine[2] == "gene":
					
				info = splitLine[8]
				splitInfo = info.split(";")
				fieldDict = dict()
				
				miRNA = False
				name = None
				for field in splitInfo:
					splitField = field.split("=")
					
				
					if splitField[1] == "miRNA":
						miRNA = True
					if splitField[0] == "gene":
						name = splitField[1]
				if miRNA == True:
					if name not in miRNAs:
						miRNAs[name] = []
				
					#Complicated chromosome encoding
					splitLocusId = splitLine[0].split(".")
					
					splitSplitLocusId = splitLocusId[0].split("_")
					
					#Now replace all 0s in the second part
					chromId = splitSplitLocusId[1].replace("0", "")
					
					miRNAs[name].append([name, 'chr' + chromId, splitLine[3], splitLine[4]])
			
			
			
			if splitLine[2] == "lnc_RNA":
				
				info = splitLine[8]
				splitInfo = info.split(";")
				fieldDict = dict()
				for field in splitInfo:
					splitField = field.split("=")
					
					fieldDict[splitField[0]] = splitField[1]
				
				
					
				if fieldDict["gene"] not in lncRNAs:
					lncRNAs[fieldDict["gene"]] = []
				
				#Complicated chromosome encoding
				splitLocusId = splitLine[0].split(".")
				
				splitSplitLocusId = splitLocusId[0].split("_")
				
				#Now replace all 0s in the second part
				chromId = splitSplitLocusId[1].replace("0", "")
				
				
				
				
				lncRNAs[fieldDict["gene"]].append([fieldDict["gene"], 'chr' + chromId, splitLine[3], splitLine[4]])



#Get the ensembl IDs to get the right gene names for the eQTL file
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
		
		ensemblIDLookup[geneId] = geneSymbol


			
#Read the eQTLs
print "reading eQTLs"
eQTLList = dict()
with open(eQTLInFile, 'rb') as inF:
	indices = []
	lineCount = 0
	header = dict()
	for line in inF:
		line = line.strip()
		
		if lineCount < 1:
			
			splitHeader = line.split("\t")
			
			for col in range(0, len(splitHeader)):
				header[splitHeader[col]] = col #keep lookup for specific columns
				
			newHeader = "chromosome\tstart\tend\tgene"
			
			for headerEntry in header:
				typeMatch = re.search("pval", headerEntry, re.IGNORECASE) #using 'chr' will match intra, inter
				if typeMatch is None: #skip the most complicated kind for now
					continue
				indices.append(header[headerEntry])
			indices = np.array(indices)
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		
		#Filter by one tissue for now. Is the p-value significant?
		#pval_Breast_Mammary_Tissue for breast tissue
		#pval_Ovary for ovary
		# tissueInd = header['pval_Prostate']
		# 
		# if splitLine[tissueInd] == 'NA':
		# 	continue
		# 
		# pval = float(splitLine[tissueInd])
		# if pval > 5e-8:
		# 	continue
		
		#Check if there is at least 1 eQTL association for this gene
		eQTLCount = 0
		
		splitLine = np.array(splitLine)
		filteredSplitLine = splitLine[indices]
		filteredSplitLine = filteredSplitLine[filteredSplitLine != 'NA']
		
		pVals = np.array(filteredSplitLine, dtype = "float")
		eQTLCount = 0
		eQTLCount += pVals[pVals < 5e-8].shape[0]
		eQTLCount += pVals[pVals == 0].shape[0]
		

		geneInfo = splitLine[0]
		splitGeneInfo = geneInfo.split(",")
		ensemblID = splitGeneInfo[1]
		splitEnsemblID = ensemblID.split(".")
		ensemblGeneID = splitEnsemblID[0]
		
		
		if ensemblGeneID in ensemblIDLookup:
			eQTLList[ensemblIDLookup[ensemblGeneID]] = eQTLCount
			# 
			# if ensemblIDLookup[ensemblGeneID] in lncRNAs:
			# 
			# 	lncRNAsWithEQTLs += 1
			
print eQTLList
lncRNAsWithEQTLs = []
miRNAsWithEQTLs = []
for gene in eQTLList:
	if gene in lncRNAs:
		lncRNAsWithEQTLs.append(gene)
	if gene in miRNAs:
		miRNAsWithEQTLs.append(gene)

print lncRNAsWithEQTLs
print miRNAsWithEQTLs
print len(lncRNAsWithEQTLs)
print len(miRNAsWithEQTLs)

