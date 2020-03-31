#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=06:00:00
#$ -e workflow_err
#$ -o workflow_out

### INSTRUCTIONS ###

#	This script intends to show the workflow used to generate all figures in the paper.

#	The workflow is listed in order. 'Required' indicates a step that is needed to generate
#	output that is used in all figures. Only skip this steps if you have already generated
#	the data, and want to re-use that in a figure.

#	If you want to skip a part of the workflow, change the 'False' to 'True' for each figure.

#	Memory requirements vary per workflow step and based on your dataset size.
#	For the HMF data, at least 16 GB of memory is required to generate the bags for MIL.

#	Each step has its own .sh file that can be used to run this step on the cluster, which
#	is named by the step it executes.


### (REQUIRED) PART 1 - DATA AND PATHS ###

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
	python "$runFolder/main.py" "" "0" "N" "$settingsFolder" "$outputFolder"
fi

### (REQUIRED) PART 3 - IDENTIFY PATHOGENIC SV-GENE PAIRS ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#first link mutations to patients. These are required to quickly check which patients
	#have which mutations in which genes
	#python "$runFolder/determinePatientGeneMutationPairs.py" "$outputFolder"

	#identify which TADs are disrupted in these patients, and compute the
	#z-scores of the genes in these TADs. Filter out genes with coding mutations.
	python "$runFolder/computeZScoresDisruptedTads.py" '../../data/genes/allGenesAndIdsHg19.txt' '/hpc/compgen/users/mnieboer/data/pipeline/read_counts/brca_tmm.txt' '/hpc/compgen/users/mnieboer/data/somatics/' "$settingsFolder" "$outputFolder"
fi

### PART 4 - SETTING UP FOR MULTIPLE INSTANCE LEARNING ###
run=false #these steps only need to be done when outputting anything related to multiple instance learning

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	python "$runFolder/normalizeBags.py" "$outputFolder"

	#then generate the similarity matrices for all SVs
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False"

fi


### FIGURE 1 ###

### FIGURE 2 ###

### FIGURE 3 - MIL PERFORMANCE CURVES PER SV TYPE ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False"

fi


### FIGURE 4 - FEATURE IMPORTANCE AND RECURRENCE ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#Generate the plotting data, and plot the feature importances for ALL instances,
	#don't split into cosmic/non-cosmic and gains/losses
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "ALL"

	#plot for GAINS only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "ALL"

	#plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "ALL"
fi



### SUPPLEMENTARY FIGURE 1 ###

### SUPPLEMENTARY FIGURE 2 ###

### SUPPLEMENTARY FIGURE 3 ###

### SUPPLEMENTARY FIGURE 4 ###

### SUPPLEMENTARY FIGURE 5 ###

### SUPPLEMENTARY TABLE 1 - FEATURE ELIMINATION ###
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#First generate the similarity matrices based on the bags in which 1 feature is shuffled each time
	#python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "True"
	#Then get the performance on all of these similarity matrices
	python "$runFolder/runMILClassifier.py" "$outputFolder" "True"
fi




### ADDITIONAL, NON-FIGURE ###

#optimizing the MIL classifiers
#simple ML
