#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=03:00:00
#$ -e svscore_err
#$ -o svscore_out

#make sure to create this conda environment, using the dependencies as needed by SVScore.
conda activate python2

fileInd="$SGE_TASK_ID"

count=1
find '/hpc/compgen/users/mnieboer/data/somatics/' -name "*fixed.vcf" | while read line ; do

	if [ $count -eq $fileInd ]; then
		fileName=$(echo $line | cut -d '/' -f 9)
		outFile="$line.svScore.vcf"

		#run SVScore
		$(perl svscore.pl -i "$line" -v > "$outFile")
		exit
	fi


	
	count=$[count + 1]
	
done
