#!/bin/bash

#### input 
#i: CHR, 1,2,...,22
#S1_outpre: output of the step 1.
#S2_outpre: dir/prefix, output of the S2.

i=$1
S1_outpre=$2
S2_outpre=$3

####SAIGE GENE R
SPAtestsR=/home1/main-R-file/step2_SPAtests.R

####ukb wes bed bim fam, by chr
bedFile=/home1/sample_qc_final/ukb_wes_chr${i}_sample_qc_final.bed
bimFile=/home1/sample_qc_final/ukb_wes_chr${i}_sample_qc_final.bim
famFile=/home1/sample_qc_final/ukb_wes_chr${i}_sample_qc_final.fam

####GLMM File from step 1
GMMATmodelFile=${S1_outpre}${i}_10PC_both.rda
varianceRatioFile=${S1_outpre}${i}_10PC_both.varianceRatio.txt

####sparseGRMFile
sparseGRMFile=/home1/UKBWES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=/home1/UKBWES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt

####SNPEff file
groupFile=/home1/UKBWES/SAIGE_group_file/lof_missense_five/SnpEff_gene_group_chr${i}.txt

####Output
SAIGEOutputFile=${S2_outpre}${i}_both_10pca.txt


## run step 2
source activate RSAIGE

  Rscript $SPAtestsR \
     --bedFile=$bedFile \
     --bimFile=$bimFile \
     --famFile=$famFile \
     --SAIGEOutputFile=$SAIGEOutputFile \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --GMMATmodelFile=$GMMATmodelFile \
     --varianceRatioFile=$varianceRatioFile \
     --sparseGRMFile=$sparseGRMFile \
     --sparseGRMSampleIDFile=$sparseGRMSampleIDFile \
     --groupFile=$groupFile \
     --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous" \
     --maxMAF_in_groupTest=0.00001,0.0001,0.001,0.01 \
     --is_output_markerList_in_groupTest=TRUE \
     --LOCO=FALSE \
     --is_fastTest=TRUE
