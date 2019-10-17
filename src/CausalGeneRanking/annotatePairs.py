"""
	Annotate each DEG pair as being affected by an SNV in the same patient, or having a pathway effect in that patient. 

"""
import sys
import numpy as np
from os import listdir
from os.path import isfile, join

#Collect the DEG pairs
degData = np.loadtxt(sys.argv[1], dtype='object')

#Collect the SNVs.

snvDir = sys.argv[2]

#Collect all patients with mutations in these genes, even though we do not have matching expression data for all of these patients
snvs = []
allFiles = [f for f in listdir(snvDir) if isfile(join(snvDir, f))]

for currentFile in allFiles:
	
	if currentFile == "MANIFEST.txt":
		continue
	splitFileName = currentFile.split(".")
	patientID = splitFileName[0]

	#Load the contents of the file
	with open(snvDir + "/" + currentFile, 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip() #remove newlines
			if lineCount < 1: #only read the line if it is not a header line
				lineCount += 1
				continue

			splitLine = line.split("\t")
		
			geneName = splitLine[0]
			
			splitPatientID = patientID.split("-")
			shortID = 'brca' + splitPatientID[2]
			
			snvs.append([geneName + "_" + shortID, splitLine[8]])

#For which DEG pair do we see an SNV in the gene that could potentially explain the DEG status? 

snvs = np.array(snvs, dtype='object')

#Split the DEG pairs and match on the SNVs
snvPairs = []
for degPair in degData:
	
	splitPair = degPair[0].split("_")
	gene = splitPair[0]
	sample = splitPair[7]
	
	if gene + "_" + sample in snvs[:,0]:
		match = snvs[snvs[:,0] == gene + "_" + sample][0]
		print(degPair[0], " has SNV ", match[1])
		snvPairs.append([degPair[0], match[1]])
	
snvPairs = np.array(snvPairs, dtype='object')

np.savetxt(sys.argv[1] + '_snvs.txt', snvPairs, delimiter='\t', fmt='%s')

#Collect pathways. For each DEG pair, find if there is another 

pathwayFile = sys.argv[3]
geneToPathway = dict() #per gene, list which pathway it is in
pathwayToGenes = dict() #for each pathway, list the genes. 
with open(pathwayFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		
		pathway = splitLine[0]
		genes = splitLine[2:]
		
		pathwayToGenes[pathway] = genes
		for gene in genes:
			if gene not in geneToPathway:
				geneToPathway[gene] = []
			geneToPathway[gene].append(pathway)

#Get the pairs that are DEG with a coding SNV and coding SV. Check which pathways these are in, and if there is a DEG in the same patient in the non-coding DEG set

#for each DEG pair, check which pathway it is in. Then check for that pathway which other genes are there. Is there a DEG in the same patient that is in the same pathway?
formattedPairs = []
for degPair in degData:
	
	splitPair = degPair[0].split("_")
	gene = splitPair[0]
	sample = splitPair[7]
	
	formattedPairs.append([gene, sample, degPair[0]])
	
formattedPairs = np.array(formattedPairs, dtype='object')	

codingSVDegPairs = np.load(sys.argv[4], allow_pickle=True, encoding='latin1')
codingSNVDegPairs = np.load(sys.argv[5], allow_pickle=True, encoding='latin1')

allCodingPairs = []
for pair in codingSVDegPairs[:,0]:
	splitPair = pair.split("_")
	allCodingPairs.append([splitPair[0], splitPair[7]])
for pair in codingSNVDegPairs[:,0]:
	splitPair = pair.split("_")
	allCodingPairs.append([splitPair[0], splitPair[1]])

allCodingPairs = np.array(allCodingPairs, dtype='object')

pathwayAnnotatedPairs = []	
for pair in formattedPairs:
	gene = pair[0]
	sample = pair[1]
	
	if gene not in geneToPathway: #gene is not annotated to any pathway
		continue
	
	pathways = geneToPathway[gene]
	
	samePatientGenes = allCodingPairs[allCodingPairs[:,1] == sample] #get the other genes that are DEG in this patient.
	if len(samePatientGenes) < 1: #stop if there are no other DEGs in this patient
		continue
	
	#for the pathways that this gene is involved in, which other genes are there?
	matchedGenes = ""
	for pathway in pathways:
		otherGenes = pathwayToGenes[pathway]
		#is any of these genes DEG in the same patient?
		for otherGene in otherGenes:
			if otherGene in samePatientGenes and otherGene != gene:
				matchedGenes += pathway
				matchedGenes += "_"
				matchedGenes += otherGene
				matchedGenes += "_"
						
	if matchedGenes != "":
		pathwayAnnotatedPairs.append([pair[2], matchedGenes[:-1]])
			
pathwayAnnotatedPairs = np.array(pathwayAnnotatedPairs, dtype='object')
print(pathwayAnnotatedPairs)
print("Genes that have DEGs in other pathways: ", pathwayAnnotatedPairs.shape, " out of total : ", degData.shape)
	
np.savetxt(sys.argv[1] + "_pathwayAnnotation.txt", pathwayAnnotatedPairs, fmt='%s', delimiter='\t')	
