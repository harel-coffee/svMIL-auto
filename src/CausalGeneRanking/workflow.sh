#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
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

#the output needs to be fixed, to where settings can access it too.

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
fi

### PART 4 - SETTING UP FOR MULTIPLE INSTANCE LEARNING ###
run=false #these steps only need to be done when outputting anything related to multiple instance learning

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	#python "$runFolder/normalizeBags.py" "$outputFolder"

	#then generate the similarity matrices for all SVs
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False"

fi


### FIGURE 1 ###

### FIGURE 2 ###
run=false

if $run; then
	runFolder='./linkSVsGenes/'

	#Map the SVs to genes. Specific for germline, settings for this are different
	#settingsFolder='./settings/settings_HMF_BRCA_germline/'
	#python "$runFolder/main.py" "germline" "False" "0" "$settingsFolder" "$outputFolder"

	#Repeat for shuffled SVs
	#settingsFolder='./settings/settings_HMF_BRCA/'
	#python "$runFolder/main.py" "random" "False" "0" "$settingsFolder" "$outputFolder"

	#split affected/non-affected pairs
	inFile="$outputFolder/linkedSVGenePairs/nonCoding_geneSVPairs.txt_"
	zScoresFile="$outputFolder/tadDisruptionsZScores/zScores.txt"
	#python "outputDegNonDegPairsFeatures.py" "$inFile" "$zScoresFile"

	#Then make the heatmap plot
	python "plotCodingNonCodingFeatures.py" "$outputFolder"

fi

### FIGURE 3 - MIL PERFORMANCE CURVES PER SV TYPE, PER-PATIENT CV ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi


### FIGURE 4 - FEATURE IMPORTANCE AND RECURRENCE ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	#python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True"

	#Generate the plotting data, and plot the feature importances for ALL instances,
	#don't split into cosmic/non-cosmic and gains/losses
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "ALL" "$settingsFolder" "False"

	#plot for GAINS only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "ALL" "$settingsFolder" "False"

	#plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "ALL" "$settingsFolder" "False"

	#needs a step to output deg/non-deg pairs

	##Recurrence figure
	#python "$runFolder/recurrenceAnalysis.py" "$outputFolder"

fi



### SUPPLEMENTARY FIGURE 1 ###

### SUPPLEMENTARY FIGURE 2 ###

### SUPPLEMENTARY FIGURE 3 ###

### SUPPLEMENTARY FIGURE 4 ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	#python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True"

	#Generate the figure. The final param will make sure only the plot is saved
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "ALL" "$settingsFolder" "True"

fi

### SUPPLEMENTARY FIGURE 5 - FEATURE IMPORTANCES SPECIFIC FOR COSMIC GENES ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#Generate the plotting data, and plot the feature importances for ALL instances,
	#don't split into cosmic/non-cosmic and gains/losses
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "COSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "COSMIC" "$settingsFolder" "False"

	#plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "COSMIC" "$settingsFolder" "False"

	#then run for non-cosmic
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "NONCOSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "NONCOSMIC" "$settingsFolder" "False"

	#plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "NONCOSMIC" "$settingsFolder" "False"

fi


### SUPPLEMENTARY TABLE 1 - FEATURE ELIMINATION ###
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#First generate the similarity matrices based on the bags in which 1 feature is shuffled each time
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "True" "False" "True" "False" "False"
	#Then get the performance on all of these similarity matrices
	python "$runFolder/runMILClassifier.py" "$outputFolder" "True" "False" "False" "False" "False"
fi

### ADDITIONAL, NON-FIGURE ###

#optimizing the MIL classifiers
#simple ML


## leave one patient out CV with random labels
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "True"

fi


## Per-chromosome CV
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "True" "False" "False"

	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False" "True" "False" "False"

fi

## leave bags out CV
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "True" "False"

	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False" "False" "True"

fi

## checking cosmic in CV loop
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False"

fi

#cosmic analysis
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	python "$runFolder/checkCosmicEnrichment.py" "$outputFolder" "$settingsFolder"

fi

#comparing to simple ML
run=false

if $run; then
	runFolder='./'

	python "$runFolder/simpleML.py" "$outputFolder" "$settingsFolder"

fi

#test with promoters, strengths etc
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'
	outputFolder='./output/test'

	#normalize bags
	python "$runFolder/normalizeBags.py" "$outputFolder"

	#chrCV
	#python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "True" "False" "False"

	#python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False" "True" "False" "False"


fi


#tad plot
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#rules, cosmicRules, expression, random

	#first a run with rules and expression cutoff
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "True" "False"

	#also re-run z-scores with random expression to get plot with randomized expression.
	#python "$runFolder/computeZScoresDisruptedTads.py" "$settingsFolder" "$outputFolder" "True"
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "True" "True"

	#then a run with all genes and cutoff, but no rules
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "False" "True" "False"

	#and finally a run with cosmic rules and expression cutoff
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "True" "True" "False"

	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "True" "False" "False"

	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "True" "True" "True"

	#rules, random SVs, no cutoff.
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "True"

	#rules,
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "False"

	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "True" "False"

	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "True" "False"


fi
