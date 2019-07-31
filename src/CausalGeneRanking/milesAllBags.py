import sys
import numpy as np
import matplotlib.pyplot as plt

# Combining DEGs, recurrence and germline

somaticScores = np.loadtxt(sys.argv[1], dtype="object")
germlineScores = np.loadtxt(sys.argv[2], dtype="object")
degLabels = np.loadtxt(sys.argv[3], dtype="object")
somaticRanks = np.loadtxt(sys.argv[4], dtype="object")
germlineRanks = np.loadtxt(sys.argv[5], dtype="object")

#get the deg genes
degs = []
for pair in degLabels[:,0]:
	splitPair = pair.split("_")
	gene = splitPair[0]
	degs.append(gene)

#1. How many genes are linked to germline SVs? And how many of these are also linked to somatic? Which are exclusive? 

somaticGenes = []
for geneSVPair in somaticScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
		
	splitScore = geneSVPair[0].split("_")
	
	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]
	

	if score == True:
		
		if splitScore[0] not in somaticGenes:
			somaticGenes.append(splitScore[0])


germlineGenes = []
for geneSVPair in germlineScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
		
	splitScore = geneSVPair[0].split("_")
	
	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	

	if score == True:
		
		if splitScore[0] not in germlineGenes:
			germlineGenes.append(splitScore[0])

print len(somaticGenes)
print len(germlineGenes)

print len(np.intersect1d(somaticGenes, germlineGenes))
somaticUniqueGenes = np.setdiff1d(somaticGenes, germlineGenes)
germlineUniqueGenes = np.setdiff1d(germlineGenes, somaticGenes)

#Of the somatic unique genes, which are recurrent and which are DEG?
recurrenceThreshold = 3
somaticRecurrentGenes = []
for gene in somaticRanks:
	
	samples = gene[31]
	
	if samples == "None":
		continue
	
	splitSamples = samples.split(",")
	sampleNum = len(splitSamples)
	
	if sampleNum >= recurrenceThreshold and gene[0] in somaticUniqueGenes:
		
		somaticRecurrentGenes.append(gene[0])

print "Number of somatic recurrent genes: ", len(somaticRecurrentGenes)
somaticUniqueDegs = np.intersect1d(somaticUniqueGenes, degs)
print "Number of somatic unique genes that are DEG: ", len(somaticUniqueDegs)

#Of the germline unique genes, which are recurrent and which are DEG? 
recurrenceThreshold = 3
germlineRecurrentGenes = []
for gene in germlineRanks:
	
	samples = gene[31]
	
	if samples == "None":
		continue
	
	splitSamples = samples.split(",")
	sampleNum = len(splitSamples)
	
	if sampleNum >= recurrenceThreshold and gene[0] in germlineUniqueGenes:
		
		germlineRecurrentGenes.append(gene[0])

print "Number of germline genes that are recurrent: ", len(germlineRecurrentGenes)
print "Number of germline unique genes that are DEG: ", len(np.intersect1d(germlineUniqueGenes, degs))

#How many SVs do we get when we only use the ones that are linked to the unique germline set?

svBagContents = dict()
genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
for geneSVPair in germlineScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
		
	splitScore = geneSVPair[0].split("_")
	
	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	

	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		geneName = splitGeneSVPairName[0]
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		sv = "_".join(splitGeneSVPairName)
		
		if sv not in genesPerBag:
			genesPerBag[sv] = []
		genesPerBag[sv].append(geneName)

		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)

negativeAllBags = [] #based on all germline genes
negativeBags = [] #based on genes unique to germline
recurrentNegativeBags = [] #based on genes unique to germline that are also recurrent
negativePairNames = []
for sv in svBagContents:
	
	for gene in genesPerBag[sv]:
		if gene in germlineGenes:
			negativeAllBags.append(svBagContents[sv])
			negativePairNames.append(gene + "_" + sv)
	
	bagLabel = 0
	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent
	negativeStatus = True
	negativeRecurrentStatus = True
	negativeAllStatus = True
	for gene in genesPerBag[sv]:	
		
		#get the patient ID
		splitSVName = sv.split("_")
		patientId = splitSVName[len(splitSVName)-1]
		
		
		#Using all genes unique to germline
		if gene not in germlineUniqueGenes:
			negativeStatus = False
		#Using all genes that are germline + recurrent (VERY likely negative)
		if gene not in germlineRecurrentGenes:
			negativeRecurrentStatus = False
			
	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent	
		
	if negativeRecurrentStatus == True:
		recurrentNegativeBags.append(svBagContents[sv])
		negativePairNames.append(gene + "_" + sv)
	
	if negativeAllStatus == True:
		negativeAllBags.append(svBagContents[sv])
		negativePairNames.append(gene + "_" + sv)
		
		
print "Number of SVs in negative set using only genes unique to germline: ", len(negativeBags)
print "Number of SVs in negative set using only genes unique to germline AND recurrent: ", len(recurrentNegativeBags)

#How many SVs do we get in the positive set?

svBagContents = dict()
genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
for geneSVPair in somaticScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
		
	splitScore = geneSVPair[0].split("_")
	
	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	

	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		geneName = splitGeneSVPairName[0]
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		sv = "_".join(splitGeneSVPairName)
		
		if sv not in genesPerBag:
			genesPerBag[sv] = []
		genesPerBag[sv].append(geneName)

		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)

positiveAllBags = [] #based on all somatic genes
positiveBags = [] #based on genes unique to somatic
recurrentPositiveBags = [] #based on genes unique to somatic that are also recurrent
degPositiveBags = [] #based on genes unique to somatic that are also DEGs
recurrentDegPositiveBags = []
positivePairNames = []
for sv in svBagContents:
	
	for gene in genesPerBag[sv]:
		if gene in somaticGenes:
			positiveAllBags.append(svBagContents[sv])
			positivePairNames.append(gene + "_" + sv)
			
	
	bagLabel = 0
	#An SV is only in the positive set if it does not affect ANY gene that is not unique to somatic/recurrent
	positiveStatus = True
	positiveRecurrentStatus = True
	positiveDegStatus = True
	positiveRecurrentDegStatus = True
	positiveAllStatus = True
	for gene in genesPerBag[sv]:
		
		
		#get the patient ID
		splitSVName = sv.split("_")
		patientId = splitSVName[len(splitSVName)-1]
		
		if gene not in somaticUniqueGenes:
			positiveStatus = False
		#Using all genes that are somatic + recurrent
		if gene not in somaticRecurrentGenes:
			positiveRecurrentStatus = False
		if gene not in somaticUniqueDegs:
			positiveDegStatus = False
		if gene not in somaticUniqueDegs or gene not in somaticRecurrentGenes:
			positiveRecurrentDegStatus = False
		
		
			
	#An SV is only in the negative set if it does not affect ANY gene that is not unique to germline/recurrent	
	
	if positiveStatus == True:
		positiveBags.append(svBagContents[sv])
	if positiveRecurrentStatus == True:
		recurrentPositiveBags.append(svBagContents[sv])
	if positiveDegStatus == True:
		degPositiveBags.append(svBagContents[sv])
		
	if positiveRecurrentDegStatus == True:
		recurrentDegPositiveBags.append(svBagContents[sv])
	
print "Number of SVs in positive set using only genes unique to somatic: ", len(positiveBags)
print "Number of SVs in positive set using only genes unique to somatic AND recurrent: ", len(recurrentPositiveBags)
print "Number of SVs in positive set using only genes unique to somatic AND DEGs: ", len(degPositiveBags)
print "Number of SVs in positive set using only genes unique to somatic AND DEGs AND recurrent: ", len(recurrentDegPositiveBags)

#Combine the bags and set labels

#Genes in somatic & germline
bags = positiveAllBags + negativeAllBags
posLabels = [1] * len(positiveAllBags)
negLabels = [0] * len(negativeAllBags)
labels = posLabels + negLabels
pairNames = positivePairNames + negativePairNames

#Genes unique to somatic & germline
# bags = positiveBags + negativeBags
# posLabels = [1] * len(positiveBags)
# negLabels = [0] * len(negativeBags)
# labels = posLabels + negLabels
# pairNames = positivePairNames + negativePairNames

#Genes unique to somatic & germline + recurrence
# bags = positiveBags + recurrentNegativeBags
# posLabels = [1] * len(positiveBags)
# negLabels = [0] * len(recurrentNegativeBags)
# labels = posLabels + negLabels
# pairNames = positivePairNames + negativePairNames

# #Genes unique to somatic + recurrence & germline + recurrence
# bags = recurrentPositiveBags + recurrentNegativeBags
# posLabels = [1] * len(recurrentPositiveBags)
# negLabels = [0] * len(recurrentNegativeBags)
# labels = posLabels + negLabels
# pairNames = positivePairNames + negativePairNames

# #Genes unique to somatic + DEG & germline + recurrence
# bags = degPositiveBags + recurrentNegativeBags
# posLabels = [1] * len(degPositiveBags)
# negLabels = [0] * len(recurrentNegativeBags)
# labels = posLabels + negLabels
# pairNames = positivePairNames + negativePairNames

# 
# #Genes unique to somatic + DEG + recurrent & germline + recurrence
# bags = recurrentDegPositiveBags + recurrentNegativeBags
# posLabels = [1] * len(recurrentDegPositiveBags)
# negLabels = [0] * len(recurrentNegativeBags)
# labels = posLabels + negLabels

# Combining DEG per-patient labels with recurrence, in OR fashion
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# degLabels = np.loadtxt(sys.argv[2], dtype="object")
# geneRanks = np.loadtxt(sys.argv[3], dtype="object")
# 
# #Determine threshold for recurrence
# recurrenceThreshold = 5
# 
# recurrentGenes = []
# 
# #How many genes match this threshold?
# for gene in geneRanks:
# 	
# 	samples = gene[31]
# 	
# 	if samples == "None":
# 		continue
# 	
# 	splitSamples = samples.split(",")
# 	sampleNum = len(splitSamples)
# 	
# 	if sampleNum >= recurrenceThreshold:
# 		
# 		recurrentGenes.append(gene[0])
# 
# print "Number of recurrent genes: ", len(recurrentGenes)
# 
# #get the deg genes
# degs = []
# for pair in degLabels[:,0]:
# 	splitPair = pair.split("_")
# 	gene = splitPair[0]
# 	degs.append(gene)
# 
# #Check overlap between recurrent genes and DEGs
# # degRecurrence = []
# # for gene in recurrentGenes:
# # 	degs.append(gene)
# 		
# 	#if gene in degs:
# 		#degRecurrence.append(gene)
# 
# # print np.unique(degs)
# # print len(np.unique(degs))
# # exit()
# # print degRecurrence
# # print len(degRecurrence)
# 
# #np.savetxt('Output/recurrentDegs.txt', degRecurrence, fmt='%s')
# 
# exit()
# 	
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 
# 	
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# 
# #Then construct the bags and add the right labels
# 
# matchedPairs = []
# bags = []
# labels = []
# pairNames = []
# pos = 0
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 		
# 		#The bag can only be negative if all instances are negative. 
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		# Recurrence + expression
# 		if gene + "_" + patientId in degLabels[:,0] or gene in recurrentGenes:
# 			bagLabel = 1
# 			matchedPairs.append(gene + "_" + patientId)
# 		
# 		#Only include a DEG if it is also recurrent
# 		# if gene + "_" + patientId in degLabels[:,0] and gene in recurrentGenes:
# 		# 	bagLabel = 1
# 		# 	matchedPairs.append(gene + "_" + patientId)
# 		# 
# 		# if gene in recurrentGenes and gene + "_" + patientId not in degLabels[:,0]:
# 		# 	bagLabel = 1
# 		# 	matchedPairs.append(gene + "_" + patientId)
# 		# 	
# 		#Recurrence only:
# 		# if gene in recurrentGenes:
# 		# 	bagLabel = 1
# 		# 	matchedPairs.append(gene + "_" + patientId)
# 
# 	if bagLabel > 0:
# 		pos += 1
# 	
# 	bags.append(svBagContents[sv])
# 	labels.append(bagLabel)
# 
# labels = np.array(labels)
# print "positive bags: ", len(np.where(labels == 1)[0])
# print "negative bags: ", len(np.where(labels == 0)[0])


# Running MILES with somatic only, using the > 3 patient DEG labels
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# degLabels = np.loadtxt(sys.argv[2], dtype="object")
# 
# 
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 	
# 
# #Then construct the bags and add the right labels
# 
# bags = []
# labels = []
# pairNames = []
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 		
# 		#The bag can only be negative if all instances are negative. 
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		if gene in degLabels[:,0]: #we already make bags for sv-gene pairs, so if the gene is not affected by an SV, the bag will not be positive. 
# 			bagLabel = 1
# 
# 		
# 	bags.append(svBagContents[sv])
# 	labels.append(bagLabel)
# labels = np.array(labels)
# print "positive bags: ", len(np.where(labels == 1)[0])
# print "negative bags: ", len(np.where(labels == 0)[0])
# exit()

# #Running MILES with somatic only, but using per patient-gene pair DEG labels
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# degLabels = np.loadtxt(sys.argv[2], dtype="object")
# 
# somaticScorePairs = []
# for score in somaticScores:
# 	splitScore = score[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]
# 	somaticScorePairs.append(pair)
# 	
# print len(np.setdiff1d(somaticScorePairs, degLabels[:,0]))
# print np.setdiff1d(degLabels[:,0], somaticScorePairs)
# 
# 
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 		
# 	splitScore = geneSVPair[0].split("_")
# 	
# 	pair = splitScore[0] + "_" + splitScore[len(splitScore)-1]	
# 	# if score == False and pair in degLabels[:,0]:
# 	# 	print pair
# 	# 	print geneSVPair
# 	# 	exit()
# 	
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# 
# #Then construct the bags and add the right labels
# 
# matchedPairs = []
# bags = []
# labels = []
# pairNames = []
# pos = 0
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 		
# 		#The bag can only be negative if all instances are negative. 
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		if gene + "_" + patientId in degLabels[:,0]:
# 			bagLabel = 1
# 			matchedPairs.append(gene + "_" + patientId)
# 		
# 
# 	if bagLabel > 0:
# 		pos += 1
# 	
# 	bags.append(svBagContents[sv])
# 	labels.append(bagLabel)
# 
# print np.setdiff1d(matchedPairs, degLabels[:,0])
# print np.setdiff1d(degLabels[:,0], matchedPairs)
# 
# 
# 
# labels = np.array(labels)
# print "positive bags: ", len(np.where(labels == 1)[0])
# print "negative bags: ", len(np.where(labels == 0)[0])
# exit()

# Germline vs somatic miles
# 
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# germlineScores = np.loadtxt(sys.argv[2], dtype="object")
# 
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# #A bag is positive if it comes from the somatic set, and negative if it comes from the germline set
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		# samples = somaticAllSampleScores[somaticAllSampleScores[:,0] == geneName, 31][0]
# 		# splitSamples = samples.split(",")
# 		# if len(splitSamples) < 7:
# 		# 	continue
# 		
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 	
# 
# #Then construct the bags
# 
# bags = []
# labels = []
# pairNames = []
# for sv in svBagContents:
# 	
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 	bags.append(svBagContents[sv])
# 	labels.append(1)
# 	
# svBagContents = dict()
# genesPerBag = dict()
# for geneSVPair in germlineScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		# samples = germlineAllSampleScores[germlineAllSampleScores[:,0] == geneName, 31][0]
# 		# splitSamples = samples.split(",")
# 		# if len(splitSamples) < 7:
# 		# 	continue
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 		
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# 
# #Then construct the bags
# 
# #bags = []
# for sv in svBagContents:
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 	bags.append(svBagContents[sv])
# 	labels.append(0)
# 
# labels = np.array(labels)
# print "positive bags: ", len(np.where(labels == 1)[0])
# print "negative bags: ", len(np.where(labels == 0)[0])

# # bags = bags[1:100]
# # labels = labels[1:100]
bags = np.array(bags)
instances = np.vstack(bags)

#np.random.shuffle(labels)

print "generating similarity matrix"

#Unfold the training bags so that we can compute the distance matrix at once to all genes
bagMap = dict()
reverseBagMap = dict()
geneInd = 0
for bagInd in range(0, bags.shape[0]):
	reverseBagMap[bagInd] = []
	for gene in bags[bagInd]:
		bagMap[geneInd] = bagInd
		reverseBagMap[bagInd].append(geneInd)
		
		geneInd += 1


similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
print "Number of bags: ", bags.shape[0]
for bagInd in range(0, bags.shape[0]):
	
	print bagInd
	
	#Get the indices of the instances that are in this bag
	instanceIndices = reverseBagMap[bagInd]
	
	instanceSubset = instances[instanceIndices,:]
	
	#Compute the pairwise distance matrix here
	minDistance = float("inf")
	minDistanceInd = 0
	for instanceInd in range(0, instanceSubset.shape[0]):
		instance = instanceSubset[instanceInd]
		distance = np.abs(instance - instances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 

		summedDistance = np.sum(distance,axis=1)

		currentMinDistance = np.min(summedDistance)
		if currentMinDistance < np.min(minDistance):
			minDistance = summedDistance
			minDistanceInd = instanceInd

	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
	similarityMatrix[bagInd] = minDistance

#Save similarity matrix
np.save("Output/similarityMatrices/somaticGermline.txt", similarityMatrix)
np.save("Output/similarityMatrices/somaticGermlineLabels.txt", labels)
np.save("Output/similarityMatrices/somaticGermlinePairNames.txt", pairNames)
exit()


# # 4. Train a classifier on the similarity space
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import auc, precision_recall_curve
# from sklearn.model_selection import StratifiedKFold
# 
# cv = StratifiedKFold(n_splits=10)
# print "training the classifier in similarity space"
# np.random.seed(500)
# 
# accs = []
# aucs = []
# coeffs = []
# predDiffs = []
# for train, test in cv.split(similarityMatrix, labels):
# 	
# 	rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
# 	rfClassifier.fit(similarityMatrix[train], labels[train]) #Use the bag labels, not the instance labels
# 
# 	predictions = rfClassifier.predict(similarityMatrix[test])
# 	precision, recall, thresholds = precision_recall_curve(labels[test], predictions)
# 	aucScore = auc(recall, precision)
# 	predsDiff = np.average(labels[test] == np.sign(predictions))
# 	#Now select the most important features with random forest
# 	importances = rfClassifier.feature_importances_
# 	std = np.std([tree.feature_importances_ for tree in rfClassifier.estimators_],
# 				 axis=0)
# 	indices = np.argsort(importances)[::-1]
# 	
# 	nonZeroIndices = []
# 	for index in indices:
# 		if importances[index] > 0:
# 			nonZeroIndices.append(index)
# 	
# 	#Get the genes at these indices
# 	positiveGenes = dict()
# 	positivePairs = []
# 	for index in nonZeroIndices:
# 		pairName = pairNames[index]
# 		positivePairs.append(pairName)
# 		splitPairName = pairName.split("_")
# 		geneName = splitPairName[len(splitPairName)-1]
# 		if geneName not in positiveGenes:
# 			positiveGenes[geneName] = 0
# 		positiveGenes[geneName] += 1
# 
# 	aucs.append(aucScore)
# 	coeffs.append(len(positiveGenes))
# 	predDiffs.append(predsDiff)
# 
# #Report the averages per alpha
# 
# print "Actual acc: ", np.mean(predDiffs)
# print "Mean AUC: ", np.mean(aucs)
# print "Mean coeffs: ", np.mean(coeffs)
# 
# 
# 
# exit()
# 
np.save("lasso5RecurrencePP/similarityMatrix.txt", similarityMatrix)
np.save("lasso5RecurrencePP/labels.txt", labels)

### Loading pre-made data to save time

#To save time, bags and labels have been stored on disk already and can be re-loaded
# bags = np.load("SomaticGermline/bags.txt.npy")
# labels = np.load("SomaticGermline/labels.txt.npy")
# pairNames = np.load("SomaticGermline/pairNames.txt.npy") #the sv-gene pair names of each bag entry
# similarityMatrix = np.load("SomaticGermline/similarityMatrix.txt.npy")

#Shuffle the labels
# np.random.shuffle(labels)


# similarityMatrix = similarityMatrix[1:100, 1:100]
# labels = labels[1:100]

# ## Using SVM
# from sklearn.metrics import auc, precision_recall_curve
# print "SVM performance: "
# 
# from sklearn.svm import LinearSVC
# clf = LinearSVC()
# print "Fitting classifier: "
# clf.fit(similarityMatrix, labels)
# print "Scoring classifier: "
# score = clf.score(similarityMatrix, labels)
# print score
# 
# #predictionProbs = clf.predict_proba(similarityMatrix)[:,1]
# preds = clf.predict(similarityMatrix)
# predsDiff = np.average(labels == np.sign(preds))
# print "mean score: ", predsDiff
# 
# precision, recall, thresholds = precision_recall_curve(labels, preds)
# aucScore = auc(recall, precision)
# 
# print "AUC: ", aucScore
# print "Used features: ", len(clf.coef_[0] != 0)
# 
# 
# #collect the pairs in this order
# pairsRanking = np.empty([len(clf.coef_[0]), 2], dtype="object")
# for coefInd in range(0, len(clf.coef_[0])):
# 	pairName = pairNames[coefInd]
# 	pairsRanking[coefInd, 0] = pairName
# 	pairsRanking[coefInd, 1] = clf.coef_[0][coefInd]
# 
# #Sort the features by importance/weights
# pairsRanking = pairsRanking[np.argsort(np.abs(pairsRanking[:,1]))[::-1]]	
# exit()
# 	
# Output the pairs in this order
# pairsRankingOut = "svm5RecurrencePP/pairsRanking_random.txt"
# with open(pairsRankingOut, 'w') as outF:
# 	for pairInd in range(0, pairsRanking.shape[0]):
# 		outF.write(pairsRanking[pairInd, 0] + "\t" + str(pairsRanking[pairInd, 1]) + "\n")

# ###Using SVM with cross validation
# print "SVM with RFECV performance: "
# from sklearn.feature_selection import RFECV
# from sklearn.metrics import auc, precision_recall_curve
# from sklearn.svm import LinearSVC
# clf = LinearSVC()
# selector = RFECV(clf, step=1, cv=5)
# print "Fitting classifier: "
# selector = selector.fit(similarityMatrix, labels)
# print "score: ", selector.score(similarityMatrix, labels)
# preds = selector.predict(similarityMatrix)
# predsDiff = np.average(labels == np.sign(preds))
# 
# print "mean score: ", predsDiff
# 
# precision, recall, thresholds = precision_recall_curve(labels, preds)
# aucScore = auc(recall, precision)
# 
# print "AUC: ", aucScore
# 
# selectedGenesInd = selector.support_
# selectedPairs = pairNames[selectedGenesInd]
# 
# print "Number of selected pairs: ", selectedPairs
# pairsRankingOut = "svmRFECVPerPatient/pairsRanking.txt"
# with open(pairsRankingOut, 'w') as outF:
# 	for pairInd in range(0, selectedPairs.shape[0]):
# 		outF.write(selectedPairs[pairInd] + "\n")

# output the similarity matrix and also the labels and bag indices so that we can do analysis after running once
# np.save("svmRFECVPerPatient/bags.txt", bags)
# np.save("svmRFECVPerPatient/pairNames.txt", pairNames) #the sv-gene pair names of each bag entry
# np.save("svmRFECVPerPatient/labels.txt", labels)
# np.save("svmSomaticGermline/similarityMatrix.txt", similarityMatrix) #no need to store this every time
# np.save("svmRFECVPerPatient/coef.txt", clf.coef_[0])


### Using Lasso
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits=10)
labels = np.array(labels)
alphas = [1e-15, 1e-10, 1e-8, 1e-4, 1e-3, 1e-2, 1, 5, 10, 20]
alphas = [1e-4, 1e-3, 1e-2, 1]
alphas = [1e-2]
# 
# accs = dict()
# aucs = dict()
# coeffs = dict()
# predDiffs = dict()
# 
# # Lasso but then with CV
# for currentAlpha in alphas:
# 	print "alpha: ", currentAlpha
# 	accs[currentAlpha] = []
# 	aucs[currentAlpha] = []
# 	coeffs[currentAlpha] = []
# 	predDiffs[currentAlpha] = []
# 	for train, test in cv.split(similarityMatrix, labels):
# 		
# 		lasso = Lasso(alpha=currentAlpha)
# 		lasso.fit(similarityMatrix[train],labels[train])
# 		
# 		#train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
# 		test_score=lasso.score(similarityMatrix[test],labels[test])
# 		coeff_used = np.sum(lasso.coef_!=0)
# 		preds = lasso.predict(similarityMatrix[test])
# 		predsDiff = np.average(labels[test] == np.sign(preds))
# 		
# 		precision, recall, thresholds = precision_recall_curve(labels[test], preds)
# 		aucScore = auc(recall, precision)
# 		
# 		accs[currentAlpha].append(test_score)
# 		aucs[currentAlpha].append(aucScore)
# 		coeffs[currentAlpha].append(coeff_used)
# 		predDiffs[currentAlpha].append(predsDiff)
# 
# 	#Report the averages per alpha
# 	
# 	print "Actual acc: ", np.mean(predDiffs[currentAlpha])
# 	print "Mean acc: ", np.mean(accs[currentAlpha])
# 	print "Mean AUC: ", np.mean(aucs[currentAlpha])
# 	print "Mean coeffs: ", np.mean(coeffs[currentAlpha])
# 	
# 	# np.save("lasso2Patients/acc.txt", np.mean(predDiffs[currentAlpha]))
# 	# np.save("lasso2Patients/preds.txt", np.mean(accs[currentAlpha]))
# 	# np.save("lasso2Patients/auc.txt", np.mean(aucs[currentAlpha]))
# 	# np.save("lasso2Patients/coeffs.txt", np.mean(coeffs[currentAlpha]))
# # 	
# exit()	

# #####Getting the concept genes without cross validation
currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
lasso.fit(similarityMatrix,labels)

#train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
test_score=lasso.score(similarityMatrix,labels)
coeff_used = np.sum(lasso.coef_!=0)
preds = lasso.predict(similarityMatrix)
predsDiff = np.average(labels == np.sign(preds))

precision, recall, thresholds = precision_recall_curve(labels, preds)
aucScore = auc(recall, precision)

import matplotlib.pyplot as plt

plt.step(recall, precision, color='b', alpha=0.2,
         where='post')
plt.fill_between(recall, precision, alpha=0.2, color='b')

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
#plt.show()
#plt.savefig('lasso2RecurrencePP/lasso.svg')


print "acc: ", test_score
print "predsDiff: ", predsDiff
print "auprc: ", aucScore
print "coeffs: ", coeff_used


geneIndices = np.where(lasso.coef_ !=0)[0]
positiveGenes = dict()
positivePairs = []
for index in geneIndices:
	pairName = pairNames[index]
	positivePairs.append(pairName)
	splitPairName = pairName.split("_")
	geneName = splitPairName[len(splitPairName)-1]
	if geneName not in positiveGenes:
		positiveGenes[geneName] = 0
	positiveGenes[geneName] += 1

print len(positivePairs)
print len(positiveGenes)
# 
# milesConceptGenesOut = "lasso2RecurrencePP/milesConceptGenes_nonRandom.txt"
# with open(milesConceptGenesOut, 'w') as outF:
# 	for gene in positiveGenes:
# 		outF.write(gene + "\t" + str(positiveGenes[gene]) + "\n")
exit()

#output the similarity matrix and also the labels and bag indices so that we can do analysis after running once
np.save("lasso2Patients/alpha.txt", currentAlpha)
np.save("lasso2Patients/bags.txt", bags)
np.save("lasso2Patients/pairNames.txt", pairNames) #the sv-gene pair names of each bag entry
np.save("lasso2Patients/labels.txt", labels)
np.save("lasso2Patients/similarityMatrix.txt", similarityMatrix)
# np.save("lassoPerPatient/conceptIndices.txt", geneIndices)

#Instead save the scores to file
# np.save("lassoPerPatient/acc_random.txt", test_score)
# np.save("lassoPerPatient/preds_random.txt", predsDiff)
# np.save("lassoPerPatient/auc_random.txt", aucScore)
# np.save("lassoPerPatient/coeffs_random.txt", coeff_used)

exit()



