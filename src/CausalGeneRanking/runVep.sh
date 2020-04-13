#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=03:00:00
#$ -e vep_err
#$ -o vep_out

vepPath='/hpc/compgen/users/mnieboer/Tools/ensembl-vep/vep'
cacheDir='/hpc/compgen/users/mnieboer/Tools/ensembl-vep/cacheDir'
outDir='/hpc/compgen/users/mnieboer/Tools/ensembl-vep/outDir'

fileInd="$SGE_TASK_ID"
#fileInd=1

count=1
find '/hpc/compgen/users/mnieboer/data/somatics/' -name "*fixed.vcf" | while read line ; do

	if [ $count -eq $fileInd ]; then
		fileName=$(echo $line | cut -d '/' -f 9)

		#run VEP
		$(perl "$vepPath" --cache --dir_cache "$cacheDir" --format vcf --port 3337 --assembly "GRCh37" -i "$line" -o "$outDir/$fileName")
		exit
	fi


	
	count=$[count + 1]
	
done
