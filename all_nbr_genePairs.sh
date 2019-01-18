#!/bin/bash

#./all_nbr_genePairs.sh

maxLoad="100"

maxJobs=7

all_datasets=($( ls "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/OUTPUT_FOLDER" ))


#all_datasets=( "GSE74927_neg_pos" )



parallel -i -j $maxJobs -l $maxLoad sh -c "echo Rscript nbr_genePairs_by_dist.R {}" -- ${all_datasets[@]}                
parallel -i -j $maxJobs -l $maxLoad sh -c "Rscript nbr_genePairs_by_dist.R {}" -- ${all_datasets[@]}                
