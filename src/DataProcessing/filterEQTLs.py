"""
	The goal of this script is to read the multi-tissue eQTLs from GTEx and filter it by genes that are known to be causal.
	The ENSEMBL IDs first need to be mapped to gene names

"""

import sys
import numpy as np
import os

eQTLFile = sys.argv[1]
geneListFile = sys.argv[2]
outDir = sys.argv[3]
tissue = sys.argv[4]

if not os.path.exists(outDir):
   os.makedirs(outDir)

filteredEQTLFile = outDir + '/' + tissue + '_eQTLs.bed'

#1. Read ensembl gene list and make lookup for causal genes only

def readEnsemblGeneFile(geneListFile):

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
	
	return ensemblIDLookup
			
ensemblIDLookup = readEnsemblGeneFile(geneListFile)


#3. Read eQTL file and only write the line back to the fitered file if it matches on a causal gene

def filterEQTLs(eQTLFile, ensemblIDLookup, filteredEQTLFile):

	with open(filteredEQTLFile, 'w') as outFile:
		with open(eQTLFile, 'r') as f:
			lineCount = 0
			header = dict()
			for line in f:
				line = line.strip()
				
				if lineCount < 1:
					splitHeader = line.split("\t")
					for col in range(0, len(splitHeader)):
						header[splitHeader[col]] = col #keep lookup for specific columns

					lineCount += 1
					continue
				
				splitLine = line.split("\t")

				#Filter by one tissue for now. Is the p-value significant?
				#pval_Breast_Mammary_Tissue for breast tissue
				#pval_Ovary for ovary
				if tissue == 'hmec':
					tissueInd = header['pval_Breast_Mammary_Tissue']
				elif tissue == 'ov':
					tissueInd = header['pval_Ovary']
				elif tissue == 'luad':
					tissueInd = header['pval_Lung']
				elif tissue == 'coad':
					tissueInd = header['pval_Colon_Sigmoid']
				elif tissue == 'gm12878':
					tissueInd = header['pval_Whole_Blood']
				elif tissue == 'prostate':
					tissueInd = header['pval_Prostate']
				elif tissue == 'esophagus':
					tissueInd = header['pval_Esophagus_Mucosa']
				elif tissue == 'skin':
					tissueInd = header['pval_Skin_Sun_Exposed_Lower_leg']
				elif tissue == 'pancreas':
					tissueInd = header['pval_Pancreas']
				elif tissue == 'uterus':
					tissueInd = header['pval_Uterus']
				elif tissue == 'nervousSystem':
					tissueInd = header['pval_Brain_Cortex']
				
				##### for this one, requires v8 of GTEx
				elif tissue == 'kidney':
					tissueInd = header['pval_Kidney_Cortex']
				if splitLine[tissueInd] == 'NA':
					continue
				
				pval = float(splitLine[tissueInd])
				if pval > 0.05:
					if pval != 0:
						continue
				
				
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
					
					newLine = "chr" + chromosome + "\t" + start + "\t" + end + "\t" + ensemblIDLookup[ensemblGeneID]
					newLine = chromosome + "\t" + start + "\t" + end + "\t" + ensemblIDLookup[ensemblGeneID]
					outFile.write(newLine)
					outFile.write("\n")
				
				
				
filterEQTLs(eQTLFile, ensemblIDLookup, filteredEQTLFile)			
			
		

