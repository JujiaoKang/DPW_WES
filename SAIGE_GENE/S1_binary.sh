#!/bin/bash


#### input 
#i: CHR, 1,2,...,22
#S1_inputpheno: phenofile.txt, header is IID,pheno,age,gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
#S1_inputtype: binary 
#S1_outpre: dir/prefix

i=$1
S1_inputpheno=$2
S1_inputtype=$3
S1_outpre=$4

####SAIGE GENE R
fitNULLGLMMR=/home1/main-R-file/step1_fitNULLGLMM.R

####sparseGRMFile
sparseGRMFile=/home1/UKBWES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=/home1/UKBWES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_3rd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt

####ukb wes bed bim fam, by chr
plinkFile=/home1/sample_qc_final/ukb_wes_chr${i}_sample_qc_final


####pheno file
phenoFile=$S1_inputpheno
phenoCol=pheno
covarColList=age,gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10
sampleIDColinphenoFile=IID
traitType=$S1_inputtype

####Output: .rdata, Input of step 2
outputPrefix=${S1_outpre}${i}_10PC_both


## run step 1
source activate RSAIGE

Rscript $fitNULLGLMMR     \
    --sparseGRMFile=$sparseGRMFile \
    --sparseGRMSampleIDFile=$sparseGRMSampleIDFile \
    --plinkFile=$plinkFile \
    --useSparseGRMtoFitNULL=FALSE   \
    --useSparseGRMforVarRatio=TRUE \
    --phenoFile=$phenoFile \
    --phenoCol=$phenoCol \
    --covarColList=$covarColList \
    --sampleIDColinphenoFile=$sampleIDColinphenoFile \
    --isCovariateOffset=FALSE \
    --traitType=$traitType       \
    --nThreads=120    \
    --isCateVarianceRatio=TRUE \
    --outputPrefix=$outputPrefix \
    --IsOverwriteVarianceRatioFile=TRUE 
