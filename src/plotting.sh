# Script to create all paper figures for svMIL2.

##### REQUIRED RUNS #####
#1. Run svMIL2 to create all output needed for figures.

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

#### PLOTTING #####

## TO DO ###

#
#run=true
#
#if $run; then
#	runFolder='./plotting/'
#	settingsFolder='./settings/settings_HMF_Breast' #we need access to ?
#
#	#figure 1 plots
#	python "$runFolder/plotFrequentCosmicGenes.py" "$settingsFolder"
#
#fi

### FIGURE 3 ###

run=false

if $run; then
	runFolder='./plotting/'

	python "$runFolder/figure3.py"
fi

### FIGURE 4A ###
run=false

if $run; then
	runFolder='./plotting/'

	python "$runFolder/makeSwapHeatmap.py"


fi

