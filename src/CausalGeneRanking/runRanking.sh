#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=00:30:00
#$ -e ranking_normal_err
#$ -o ranking_normal_out

#Run 1 normal run without permutations


python main.py "$1" "$2" "$3" "$4" "$SGE_TASK_ID"
