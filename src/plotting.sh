#!/bin/bash
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH -o re.out
#SBATCH -e re.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.m.nieboer@umcutrecht.nl

# Script to create all paper figures for svMIL2.

# Assumes that preprocessing.sh has already been run.

#Change 'false' to 'true' for each section to generate that output.


##### REQUIRED RUNS #####
#1. In createSettingsFiles.py, set all correct paths.
#2. Run svMIL2 to create all output needed for figures.

run=false
if $run; then
	#These are all our cancer types
	types=('Breast' 'Colorectal' 'Esophagus' 'Lung' 'Kidney' 'NervousSystem' 'Ovary' 'Pancreas' 'Prostate' 'Skin' 'Uterus' 'UrinaryTract')

	#This will run all swaps, including the 'normal' non-swap run, e.g. Breast + hmec
	#To run only the non-swaps, see runSwaps.sh
	for type in ${types[@]}; do
		sh runSwaps.sh "$type"
	done
fi

#2. Run svMIL2 with CTCF loops

#To get the chromatin loops, run iTAD as described in the paper.
#Then in a settings file, use the final_iTADs_ctcf_score.bed12 output file in place
#of the tadFile.
#Run workflow_breast_ctcf.sh, workflow_lung_ctcf.sh and workflow_colorectal_ctcf.sh
#to get the results when using CTCF data.

#### PLOTTING #####


### FIGURE 1C, 4B ###
run=true

if $run; then
	runFolder='./plotting/'

	python "$runFolder/plotAUC.py"
fi

### FIGURE 2A, 2B, S4 ###
run=true

if $run; then
	runFolder='./plotting/'

	python "$runFolder/figure2.py"
fi

### FIGURE 3, 4C ###

run=true

if $run; then
	runFolder='./plotting/'

	python "$runFolder/figure3.py"
fi

### FIGURE 5 ###
run=true

if $run; then
	runFolder='./plotting/'

	python "$runFolder/makeSwapHeatmap.py"

fi

### FIGURE S2 ###

#1. First run the old svMIL model and output normalizedBags.pkl.
#2. Use these bags as input for plotVariances.py.

run=true

if $run; then
	runFolder='./plotting/'
	#specify directory with the normalizedBags.pkl
	outDir='./output/HMF_Breast_oldFeatures'

	python "$runFolder/plotVariances.py" "$outDir"

fi

### FIGURE S3, S5 and 4A ###
run=true

if $run; then
	runFolder='./plotting/'
	settingsFile='./settings/settings_HMF_Breast_hmec/'

	python "$runFolder/plotSVStats.py" "$settingsFile"

fi

### Figure S6 & S7 ###
run=true

if $run; then
	runFolder='./plotting/'
	python "$runFolder/plotFeatureImportancesWheel.py"

fi

### FIGURE S1 & S8 ###
#Generate S1 by running svMIL and svMIL2 on the breast samples, then combining
#the ROC curves in output/specifiedOutputFolder/rocCurves

#Generate S8 by running svMIL2 on the breast samples with a bag limit of 700 in
#the settings vs no limit (e.g. 70000), and combine the ROC curves
#in output/specifiedOutputFolder/rocCurves

