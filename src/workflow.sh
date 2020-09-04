#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=24:00:00
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
settingsFolder='./settings/settings_HMF_BRCA_generic/'

#Create a folder in which all output for this data will be stored
#Different steps will create their own intermediate folders in here
outputFolder='output/HMF_BRCA_generic'

### (REQUIRED) PART 2 - LINK SVS TO GENES ###
run=true #Only skip this step if all output has already been generated!

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
run=true #these steps only need to be done when outputting anything related to multiple instance learning

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first normalize the bags
	python "$runFolder/normalizeBags.py" "$outputFolder"

	#then generate the similarity matrices for all SVs
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "True" "False" "False" "False" "$settingsFolder"

fi

### FIGURE 2 - 2A-D ###
run=false

if $run; then
	runFolder='./tadDisruptionsZScores/'

	#make sure to also have random z-scores ready.

	#run with super enhancers (PANEL B), and shuffled expression (PANEL C)
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "False" "se"
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "True" "False" "False" "True" "se"

	#and run for all genes, no rules (PANEL A)
	python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "False" "False" "False" "all"
	#python "$runFolder/plotDisruptedTadZScores.py" "$outputFolder" "$settingsFolder" "False" "False" "False" "True" "all"

	
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
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#test the classifier and output the MIL curves
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "False"

fi

### FIGURE 3 - 3B-C: MIL PCAWG results ###

#Run respective workflows for PCAWG data: workflow_pcawg_ov.sh and workflow_pcawg_liver.sh

### FIGURE 3 - 3D: METHOD COMPARISON PRE-PROCESSING ###
run=false

if $run; then
	runFolder='./methodComparison/'

	#Fix the VCFs so that VEP and SVScore can use them
	python "$runFolder/fixVCFs.py" "$settingsFolder"

	#First run VEP in parallel
	qsub -t 1-181 -tc 75 "$runFolder/runVEP.sh"

	#Then also run SVScore in parallel
	#This requires some additional setup; see notes in script
	qsub -t 1-181 -tc 75 "$runFolder/runSVScore.sh"

fi

### FIGURE 3 - 3D: METHOD COMPARISON GENERATING FIGURE ###
run=false

if $run; then
	runFolder='./methodComparison/'

	#Get the TPRs and FPRs of each method
	#For non-MIL-RF
	python "$runFolder/simpleML.py" "$outputFolder"

	#For VEP
	#Should have been a better output folder in hindsight
	#make sure it matches the one in runVEP.sh.
	vepOutDir='/hpc/compgen/users/mnieboer/Tools/ensembl-vep/outDir'
	python "$runFolder/getVEPTprFpr.py" "$vepOutDir"

	#For SVScore
	python "$runFolder/getSVScoreTprFpr.py" "$settingsFolder"

	#Then gather all the TPRs and FPRs output from these method and add them in this script.
	#The TPR and FPR for svMIL are output by the leave-one-patient-out CV (for figure 3A)
	#was filled in by hand for simplicity.
	python "$runFolder/plotResultsGraph.py" "$outputFolder"

fi


### FIGURE 4 - FEATURE IMPORTANCE AND RECURRENCE ###
run=true

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True" "$settingsFolder"

	#Generate the plotting data, and plot the feature importances for ALL instances,
	#don't split into gains/losses yet. We need these for normalization.
	#First we make just the plotting data. We make the actual plots themselves later, which is important
	#for figure S4 later on to normalize properly.
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "ALL" "$settingsFolder" "False"
	#
	##plot for GAINS only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "ALL" "$settingsFolder" "False"
	#
	##plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "ALL" "$settingsFolder" "False"
	#
	##make the actual plots
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "ALL" "ALL" "$settingsFolder" "False"
	#
	##plot for GAINS only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "GAIN" "ALL" "$settingsFolder" "False"
	#
	##plot for LOSSES only
	#python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "LOSS" "ALL" "$settingsFolder" "False"
	#
	##Recurrence figure
	##This depends on the pair labels from the ALL feature importance run.
	#python "$runFolder/recurrenceAnalysis.py" "$outputFolder"

fi


### SUPPLEMENTARY FIGURE 1 - NUMBER OF GENES PUTATIVELY DISRUPTED PER TAD ###
run=false

if $run; then
	runFolder='./linkSVsGenes/'

	python "$runFolder/getGeneCountPerTad.py" "$settingsFolder" "$outputFolder"

fi

### SUPPLEMENTARY FIGURE 2 - lopoCV random labels, chrSV and leave-bags out CV ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#Get the performance with random labels using leave-one-patient-out CV
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "True" "False" "False" "True"

	#Get the performance using leave-one-chromosome-out CV
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "True" "False" "False"
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False" "True" "False" "False"

	#Get the performance using leave-bags-out CV
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "True" "False"
	python "$runFolder/runMILClassifier.py" "$outputFolder" "False" "False" "False" "True" "False"

fi


### SUPPLEMENTARY FIGURE 3 - Full feature importance plots ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#first output the full similarity matrix to train the classifier on the whole dataset.
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True"

	#Generate the figure. The final param will make sure only the plot is saved
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "ALL" "$settingsFolder" "True"

fi

### SUPPLEMENTARY FIGURE 4 - FEATURE IMPORTANCES SPECIFIC FOR COSMIC GENES ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#Generate the plotting data, and plot the feature importances for ALL instances,
	#don't split into cosmic/non-cosmic and gains/losses

	#First generate just the data for the plots, so that the normalization works out here (p-values of cosmic loss are smaller than the all set)
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "COSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "COSMIC" "$settingsFolder" "False"

	#plot for LOSSES only (this will fail because there are no losses for cosmic genes.)
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "COSMIC" "$settingsFolder" "False"

	#Then make the actual plots
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "ALL" "COSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "GAIN" "COSMIC" "$settingsFolder" "False"

	#plot for LOSSES only (this will fail because there are no losses for cosmic genes.)
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "LOSS" "COSMIC" "$settingsFolder" "False"

	#then run for non-cosmic
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "ALL" "NONCOSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "GAIN" "NONCOSMIC" "$settingsFolder" "False"

	#plot for LOSSES only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "True" "LOSS" "NONCOSMIC" "$settingsFolder" "False"
	
	#and also make the plots
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "ALL" "NONCOSMIC" "$settingsFolder" "False"

	#plot for GAINS only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "GAIN" "NONCOSMIC" "$settingsFolder" "False"

	#plot for LOSSES only
	python "$runFolder/plotFeatureImportances.py" "$outputFolder" "False" "LOSS" "NONCOSMIC" "$settingsFolder" "False"

fi

### ADDITIONAL, NON-FIGURE ###

#optimizing the MIL classifiers
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	#First generate the similarity matrices based on the bags in which 1 feature is shuffled each time
	python "$runFolder/generateSimilarityMatrices.py" "$outputFolder" "False" "False" "False" "False" "True"

	#Then do optimization of the classifiers.
	python "$runFolder/optimizeClassifiers.py" "$outputFolder"

fi


### COSMIC ANALYSIS WITHIN LEAVE-ONE-PATIENT-OUT-CV ###
run=false

if $run; then
	runFolder='./multipleInstanceLearning/'

	python "$runFolder/checkCosmicEnrichment.py" "$outputFolder" "$settingsFolder"

fi
