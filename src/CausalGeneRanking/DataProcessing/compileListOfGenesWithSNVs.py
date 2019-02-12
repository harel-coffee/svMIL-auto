"""
	The goal of this script is to take a list of genes and positions as input, and also a list of SNVs.
	The output is a list of genes that are affected by an SNV at least once (in at least 1 patient).

"""

import sys
import numpy as np
sys.path.insert(0, '../')
from inputParser import InputParser
import settings

#Read the genes
causalGenes = InputParser().readCausalGeneFile('../' + settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile('../' + settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes for now
genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#Read the SNV data
snvFile = settings.files['snvFile']
snvData = InputParser().getSNVsFromFile('../' + snvFile)


#Make list of genes with SNVs

for gene in genes:
	
	#SNVs should be after the start of the gene, but before the end of the gene.
	
	1
	
