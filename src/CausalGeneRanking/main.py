"""
	Set of scripts intended to do ranking of (causal) genes based on their neighborhood and the influence that SVs have in this neighborhood. 

	The idea is to look at all known causal genes. The neighborhood may consist of eQTLs that have an effect on this gene, or TADs directly neighboring the gene.
	If we annotate causal genes with these effects, then we can check for potential effects if anything in the neighborhood is affected by SVs.
	Using some magical statistics, we can then make a prediction on how much more likely the effects on the neighborhood are to be disruptive to the causal genes than expected by random chance.
	Later, this format can be extended to non-causal genes as well. 


	The setup of these scripts will initally be (it will not be too pretty and perfect at first):
	
	- The main script where the input (genes) are parsed and the relevant scripts are called
	- The getNeighborhood class, where things such as neighboring TADs, related eQTLs and overlapping SVs are obtained
	- The RankGenes class, where the genes are ranked by how likely they are causal for the disease depending on the SVs affecting their neighborhood
	- The output will be a matrix with all causal genes in the columns, and the relevant SVs in the rows.
	
	
	Using a gene-based approach will likely be quicker than an SV-based approach, and we can get the relevant SVs through the genes. If an SV is never affecting any of our features defined as interesting, there is no
	need to look at that SV at all. This idea may change along the way.
	

"""
import sys

from gene import Gene
from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking

#1. Read and parse the causal genes

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
			
			gene = Gene(geneSymbol, chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
			
			cosmicGenes.append(gene)
			
	
	#cosmicGenes = np.array(cosmicGeneList, dtype='object')
	
	return cosmicGenes

causalGeneFile = sys.argv[1]
causalGenes = readCausalGeneFile(causalGeneFile)

#2. Get the neighborhood for these genes
NeighborhoodDefiner(causalGenes) 

#3. Do simple ranking of the genes and report the causal SVs
GeneRanking(causalGenes)


		
	

