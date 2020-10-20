#!/bin/bash
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH -o workflow.out
#SBATCH -e workflow.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.m.nieboer@umcutrecht.nl

run=false

if $run; then
	runFolder='./multipleInstanceLearning/'
	settingsFolder='./settings/settings_HMF_Breast'

	#figure 1 plots
	python "$runFolder/plotFrequentCosmicGenes.py" "$settingsFolder"

fi

run=true

if $run; then
	runFolder='./multipleInstanceLearning/'
	settingsFolder='./settings/settings_HMF_Breast'

	#figure 2

	#first run the classifier to generate full similarity matrices

fi