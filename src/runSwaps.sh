#run swaps, and wait each time until another job has finished.

#initial test with kidney
ownType='Kidney'

#1. Generate the settings files automatically
python createSettingsFiles.py "$ownType"

#2. Do the runs
#types=('hmec' 'ov' 'gm12878' 'coad' 'luad' 'urinaryTract' 'prostate' 'esophagus' 'skin' 'pancreas' 'uterus' 'nervousSystem' 'kidney')

types=('gm12878' 'pancreas')
RES='None'
for type in ${types[@]}; do

	if [ "$type" = "$ownType" ]; then
		continue
	fi

	runFolder='./workflows/'
	settingsFile="./settings/settings_HMF_${ownType}_${type}"
	outputFile="./output/HMF_${ownType}_${type}"

	#perform run
	if [ "$RES" = "None" ]; then
		RES=$(sbatch --parsable "$runFolder/"workflow_"$ownType".sh "$settingsFile" "$outputFile")

		continue
	fi

	#RES=$(sbatch -d afterok:${RES} --parsable "$runFolder/"workflow_"$ownType".sh "$settingsFile" "$outputFile")
	


done
