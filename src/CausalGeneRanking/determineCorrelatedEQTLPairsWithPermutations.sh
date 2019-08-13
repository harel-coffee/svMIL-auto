#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=00:30:00
#$ -e eQTLCorrelation_err
#$ -o eQTLCorrelation_out

#Run 1 normal run without permutations

#Then run the permutation runs 1000 times, create jobs for that. 


uuid=$(uuidgen)

permutationsYN="False"
shuffleEQTLs="False"

sh determineCorrelatedEQTLPairs.sh "$permutationsYN" "$shuffleEQTLs" "$uuid"

permutationsYN="False"
shuffleEQTLs="True"

qsub -t 1-100:1 -tc 100 determineCorrelatedEQTLPairs.sh "$permutationsYN" "$shuffleEQTLs" "$uuid" #use job array and inside this script SGE TASK ID to get the run number
