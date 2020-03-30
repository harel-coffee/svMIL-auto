#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=20G
#$ -l h_rt=10:00:00
#$ -e optimize_err_dup_delOpt
#$ -o optimize_out_dup_delOpt

python main.py '../Output/RankedGenes/03032020/BRCA/bags.pkl' '../pValues_allGenes_smallTads.txt' 'DUP'
