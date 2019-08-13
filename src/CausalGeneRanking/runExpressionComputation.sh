#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=06:00:00
#$ -e expression_err
#$ -o expression_out


#run expression computation based on file order task id

python computeSVGenePairExpression.py "Output/geneSVPairs_somatic_me_12072019_shuffled.txt_$SGE_TASK_ID" "Output/geneCodingSVPairs_somatic_me_12072019_shuffled.txt_$SGE_TASK_ID" "../../data/expression/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
