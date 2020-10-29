#!/bin/bash
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH -o workflow.out
#SBATCH -e workflow.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.m.nieboer@umcutrecht.nl

### INSTRUCTIONS ###

#	This script intends to show the workflow used to generate all figures in the paper.

#	The workflow is listed in order. 'Required' indicates a step that is needed to generate
#	output that is used in all figures. Only skip this steps if you have already generated
#	the data, and want to re-use that in a figure.

#	If you want to skip a part of the workflow, change the 'false' to 'true' for each step.

#	Memory requirements vary per workflow step and based on your dataset size.
#	For the HMF data, at least 16 GB of memory is required to load the bags for MIL in memory.


### (REQUIRED) PART 1 - DATA AND PATHS ###

#Load in a settings file with the paths to the data that all code will be run on.
## path here
settingsFolder='./settings/settings_HMF_Ovary'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/HMF_Ovary'

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
	python "$runFolder/determinePatientGeneMutationPairs.py" "$settingsFolder" "$outputFolder"

	#identify which TADs are disrupted in these patients, and compute the
	#z-scores of the genes in these TADs. Filter out genes with coding mutations.
	python "$runFolder/computeZScoresDisruptedTads.py" "$settingsFolder" "$outputFolder" "False"
	
	#split the SV-gene pairs into pathogenic/non-pathogenic, which we use later on.
	runFolder='./linkSVsGenes/'
	python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"
	
fi

### PART 4 - SETTING UP FOR MULTIPLE INSTANCE LEARNING ###
run=false #these steps only need to be done when outputting anything related to multiple instance learning

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	python "$runFolder/normalizeBags.py" "$outputFolder"

	#then generate the similarity matrices for all SVs in the lopoCV case
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False" "$settingsFolder"

	#also generate the similarity matrices for all SVs, no splitting for lopoCV
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True" "$settingsFolder"


fi

### PART 5 - PER-PATIENT CV FOR MIL PERFORMANCE ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi

### CLEANUP LARGE SIMILARITY MATRICES ###
run=false

if $run; then

	rm "$outputFolder/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/similarityMatrix*"

fi
