#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
#$ -e workflow_err
#$ -o workflow_out

### INSTRUCTIONS ###

#	This script intends to show the workflow used to generate all figures in the paper.

#	The workflow is listed in order. 'Required' indicates a step that is needed to generate
#	output that is used in all figures. Only skip this steps if you have already generated
#	the data, and want to re-use that in a figure.

#	If you want to skip a part of the workflow, change the 'true' to 'false' for each step.

#	Memory requirements vary per workflow step and based on your dataset size.
#	For the HMF data, at least 16 GB of memory is required to load the bags for MIL in memory.


### (REQUIRED) PART 1 - DATA AND PATHS ###

#Load in a settings file with the paths to the data that all code will be run on.
## path here
settingsFolder='./settings/settings_HMF_BRCA_CRD/'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/HMF_BRCA_CRD'

### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=false #Only skip this step if all output has already been generated!

if $run; then
	runFolder='./linkSVsGenes/'
	#Map the SVs to genes. This also outputs bags for MIL.
	python "$runFolder/main.py" "" "False" "0" "$settingsFolder" "$outputFolder"

fi

### (REQUIRED) PART 3 - IDENTIFY PATHOGENIC SV-GENE PAIRS ###
run=true

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

	#then generate the similarity matrices for all SVs
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False"

fi

### FIGURE 2 - 2A-D ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#make sure to also have random z-scores ready.

	#run with super enhancers (PANEL B), and shuffled expression (PANEL C)
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "False" "se"
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "True" "se"

	#and run for all genes, no rules (PANEL A)
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "False" "False" "False" "all"
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "False" "False" "True" "all"


fi


### FIGURE 2 - 2E: HEATMAP ###
run=false

if $run; then
	runFolder='./linkSVsGenes/'

	#Map the SVs to genes. Specific for germline, settings for this are different
	settingsFolder='./settings/settings_HMF_BRCA_germline/'
	python "$runFolder/main.py" "germline" "False" "0" "$settingsFolder" "$outputFolder"

	#Repeat for shuffled SVs
	settingsFolder='./settings/settings_HMF_BRCA/'
	python "$runFolder/main.py" "random" "True" "0" "$settingsFolder" "$outputFolder"

	#split affected/non-affected pairs
	python "$runFolder/splitPairsPathogenicNonPathogenic.py" "$outputFolder"

	#Then make the heatmap plot
	python "$runFolder/plotPathogenicNonPathogenicFeatures.py" "$outputFolder"

fi

### FIGURE 3 - 3A: MIL PERFORMANCE CURVES PER SV TYPE, PER-PATIENT CV ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi

### FIGURE 3 - 3B-C: MIL PCAWG results ###

#Run respective workflows for PCAWG data: workflow_pcawg_ov.sh and workflow_pcawg_liver.sh

