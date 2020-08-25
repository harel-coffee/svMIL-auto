#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH -o workflow.out
#SBATCH -e workflow.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.m.nieboer@umcutrecht.nl

#1. Step 1: preparing everything to get the HMF_BRCA training data in place

#Load in a settings file with the paths to the data that all code will be run on.
## path here
settingsFolder='./settings/settings_HMF_BRCA/'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/HMF_BRCA'

### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=false #Only skip this step if all output has already been generated!

if $run; then
	runFolder='./linkSVsGenes/'
	#Map the SVs to genes. This also outputs bags for MIL.
	python "$runFolder/main.py" "" "False" "0" "$settingsFolder" "$outputFolder"

fi

### (REQUIRED) PART 3 - IDENTIFY PATHOGENIC SV-GENE PAIRS ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#first link mutations to patients. These are required to quickly check which patients
	#have which mutations in which genes
	#python "$runFolder/determinePatientGeneMutationPairs.py" "$settingsFolder" "$outputFolder"

	#identify which TADs are disrupted in these patients, and compute the
	#z-scores of the genes in these TADs. Filter out genes with coding mutations.
	python "$runFolder/computeZScoresDisruptedTads.py" "$settingsFolder" "$outputFolder" "False"

	#split the SV-gene pairs into pathogenic/non-pathogenic, which we use later on.
	runFolder='./linkSVsGenes/'
	python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"

fi

###  ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True" "$settingsFolder"
fi

#2. Step 2: preparing everything for PCAWG OV

#Load in a settings file with the paths to the data that all code will be run on.
## path here
settingsFolder='./settings/settings_PCAWG_OV/'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/PCAWG_OV'

### process the clustering ###

run=false
if $run; then
	runFolder='./DataProcessing/'
	type='gm12878'

	#Cluster histones
	python "$runFolder/genericClustering.py" "../data/histones/$type/"

	#Cluster DNASe
	python "$runFolder/genericClustering.py" "../data/dnase/$type/"

	#cluster RNA pol
	python "$runFolder/genericClustering.py" "../data/rnapol/$type/"

	#Cluster CTCF
	python "$runFolder/genericClustering.py" "../data/ctcf/$type/"

	#cluster eQTLs
	python "$runFolder/eQTLClustering.py" "../data/eQTLs/$type/"


fi


### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=false #Only skip this step if all output has already been generated!

if $run; then
	runFolder='./linkSVsGenes/'
	#Map the SVs to genes. This also outputs bags for MIL.
	python "$runFolder/main.py" "" "False" "0" "$settingsFolder" "$outputFolder"

fi

### (REQUIRED) PART 3 - IDENTIFY PATHOGENIC SV-GENE PAIRS ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#first link mutations to patients. These are required to quickly check which patients
	#have which mutations in which genes
	#python "$runFolder/determinePatientGeneMutationPairs.py" "$settingsFolder" "$outputFolder"

	#identify which TADs are disrupted in these patients, and compute the
	#z-scores of the genes in these TADs. Filter out genes with coding mutations.
	python "$runFolder/computeZScoresDisruptedTads.py" "$settingsFolder" "$outputFolder" "False"

	#split the SV-gene pairs into pathogenic/non-pathogenic, which we use later on.
	runFolder='./linkSVsGenes/'
	python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"

fi

###  ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'
	
	#first normalize the bags
	#python "$runFolder/normalizeBags.py" "$outputFolder"

	#first output the full similarity matrix to train the classifier on the whole dataset.
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True" "$settingsFolder"
fi

#run lopoCV and get the predictions
run=true #these steps only need to be done when outputting anything related to multiple instance learning

if $run; then
	runFolder='./multipleInstanceLearning/'

	#then generate the similarity matrices for all SVs
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False" "$settingsFolder"

fi

run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi


#3. Step 3: apply HMF-trained classifier to PCAWG OV.
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'
	outputFolderBRCA='output/HMF_BRCA'
	outputFolderPCAWG='output/PCAWG_OV'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	python "$runFolder/applyClassifier.py" "$outputFolderBRCA" "$outputFolderPCAWG"
fi



