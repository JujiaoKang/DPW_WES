#!/bin/bash
#PBS -l ncpus=1	
#PBS -l mem=100G		
#PBS -j oe

S1 1 /home1/udata/phenofile.txt quantitative /home1/AUDIT/AUDIT_T/GLMM/AUDIT |tee /home1/AUDIT/AUDIT_T/log/S1_chr1.log
S2 1 /home1/AUDIT/AUDIT_T/GLMM/AUDIT /home1/AUDIT/AUDIT_T/SPAtests/AUDIT |tee /home1/AUDIT/AUDIT_T/log/S2_chr1.log
