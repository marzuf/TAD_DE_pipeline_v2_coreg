#!/bin/bash

#./run_AUC_coexprDist_withFam_sortNoDup.sh

maxLoad="100"

maxJobs=7

all_datasets=($( ls "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER" ))

parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript AUC_coexprDist_withFam_sortNoDup.R {} hgnc" -- ${all_datasets[@]}                
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript AUC_coexprDist_withFam_sortNoDup.R {} hgnc" -- ${all_datasets[@]}                
