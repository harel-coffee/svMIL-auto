"""
	Find out which genes are only found when looking at non-coding effects, which are mixed, and which are coding-only. 
	Later also add coding SNVs. 
"""

import sys
import numpy as np


nonCodingRanks = np.loadtxt(sys.argv[1], dtype="object")
mixedRanks = np.loadtxt(sys.argv[2], dtype="object")
geneExpression = "" #Which values to use? 

#Readfor each gene if it is linked to an SV. If yes, include it in the results.
effects = dict() #store for every gene what it is affected by

for rank in nonCodingRanks:
	
	gene = rank[0]
	samples = rank[31]
	if samples != "None": #Only if the gene has been linked to at least 1 SV, we include it in the effects table. 
		if gene not in effects:
			effects[gene] = []
		effects[gene].append("Non-coding")
	
for rank in mixedRanks:
	
	gene = rank[0]
	samples = rank[31]
	if samples != "None": #Only if the gene has been linked to at least 1 SV, we include it in the effects table. 
		if gene not in effects:
			effects[gene] = []
		effects[gene].append("Mixed")

for rank in mixedRanks:
	gene = rank[0]
	samples = rank[31]
	if float(rank[3]) > 0:
		if gene not in effects:
			effects[gene] = []
		effects[gene].append("Coding")
	if samples != "None":
		if gene not in effects:
			effects[gene] = []
		effects[gene].append("Non-coding")


mixedOnly = 0
ncOnly = 0
codingOnly = 0

effectsTable = np.empty([len(effects), 3], dtype="object") #Turn dictionary into readable tsv table
effectsTable[:,0] = effects.keys()
for geneInd in range(0, len(effects)):
	gene = effects.keys()[geneInd]
	if "Non-coding" in effects[gene]:
		effectsTable[geneInd,1] = 1
	else:
		effectsTable[geneInd,1] = 0
	
	if "Coding" in effects[gene]:
		effectsTable[geneInd,2] = 1
	else:
		effectsTable[geneInd,2] = 0
		
	if effectsTable[geneInd,1] == 1 and effectsTable[geneInd,2] == 0: #in non-coding, but not coding
		ncOnly += 1
		
	if effectsTable[geneInd,1] == 0 and effectsTable[geneInd,2] == 1: #in coding, but not in non-coding
		codingOnly += 1

	if effectsTable[geneInd,1] == 1 and effectsTable[geneInd,2] == 1: #in coding, but not in non-coding
		mixedOnly += 1

print mixedOnly
print ncOnly
print codingOnly

effectsTable = effectsTable[effectsTable[:,1].argsort()][::-1]
		
np.savetxt("Output/effectsTable.txt", effectsTable, fmt="%s", delimiter="\t")

		