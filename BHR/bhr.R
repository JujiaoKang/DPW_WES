
###BHR_input_quantitative.R
library(bhr)
library(tidyr)
library(readr)
library(data.table)
###
single_result_path="/home1/BHR/single/result_sum"
###
setwd(single_result_path)
filelist <- list.files(single_result_path)
for (i in 1:length(filelist)) {
  pheno_path<-filelist[i]
  pheno_single<-as.data.frame(fread(pheno_path))
  pheno_name<-gsub("total.","",pheno_path)
  pheno_name<-gsub(".sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed_british_10pca_3rd_new_group_v2.txt","",pheno_name)
  
  
  func<-as.data.frame(fread("/home1/ukb_wes_chr_all_with_function_snpeff_chr_all_out_extracted_all_sorted.txt",header = F))
  colnames(func)<-c("MarkerID","symbol","gene","func")
  pheno_single2<-merge(pheno_single,func,by="MarkerID",all.x=T)
  
  
  ##plof
  BHR_input_plof<-subset(pheno_single2,pheno_single2$func=="1_pLOFs")
  BHR_input_plof<-BHR_input_plof[,c(16,2,3,14,9,7)]
  BHR_input_plof$phenotype_key<-pheno_name
  colnames(BHR_input_plof)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
  write_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  write.table(BHR_input_plof,write_path,sep = "\t",row.names = F,quote = F)
  
  ##missense
  BHR_input_missense<-subset(pheno_single2,pheno_single2$func=="2_missense")
  BHR_input_missense<-BHR_input_missense[,c(16,2,3,14,9,7)]
  BHR_input_missense$phenotype_key<-pheno_name
  colnames(BHR_input_missense)<-c("gene","chromosome","gene_position","N","beta","AF","phenotype_key")
  write_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  write.table(BHR_input_missense,write_path,sep = "\t",row.names = F,quote = F)
  
}


library(bhr)
library(tidyr)
library(data.table)
###
single_result_path="/home1/BHR/single/result_sum"
###
setwd(single_result_path)
filelist <- list.files(single_result_path)
for (i in 1:length(filelist)) {
  pheno_path<-filelist[i]
  pheno_name<-gsub("total.","",pheno_path)
  pheno_name<-gsub(".sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed_british_10pca_3rd_new_group_v2.txt","",pheno_name)
  
  ########################################################################################################################################################################################
  cat("################################################### plof")
  cat("\n################################################### ultra-rare##\n")
  ## ultra-rare
  
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-5 )
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_ultra_rare.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  ur_lof_h2=pheno_univariate_rare
  ur_lof=BHR_input_rare
  
  cat("\n################################################### rare1 1e-04-1e-05##\n")
  ## rare1 1e-04-1e-05
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-4 & BHR_input$AF>1e-5)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_1e-04-1e-05.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r0_lof_h2=pheno_univariate_rare
  r0_lof=BHR_input_rare
  
  cat("\n################################################### rare1 1e-03-1e-04##\n")
  ## rare1 1e-03-1e-04
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-3 & BHR_input$AF>1e-4)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_1e-03-1e-04.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r1_lof_h2=pheno_univariate_rare
  r1_lof=BHR_input_rare
  
  cat("\n################################################### rare2 1e-02-1e-03##\n")
  ## rare2 1e-02-1e-03
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-2 & BHR_input$AF>1e-3)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_1e-02-1e-03.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r2_lof_h2=pheno_univariate_rare
  r2_lof=BHR_input_rare
  
  
  cat("\n################################################### rare aggregate##\n")
  ## rare aggregate
  
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(r0_lof,r1_lof,r2_lof),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_rare.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r_lof_h2=pheno_univariate_rare
  
  cat("\n################################################### ultra-rare + rare aggregate##\n")
  ## ultra-rare + rare aggregate
  
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(ur_lof,r0_lof,r1_lof,r2_lof),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_plof_urr.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  urr_lof_h2=pheno_univariate_rare
  
  
  ######################################################################################################################################################################################################################
  cat("\n################################################### missense\n")
  cat("\n###################################################  ultra-rare##\n")
  ## ultra-rare
  
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-5 )
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_ultra_rare.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  ur_mis_h2=pheno_univariate_rare
  ur_mis=BHR_input_rare
  
  
  cat("\n###################################################  rare1 1e-04-1e-05##\n")
  ## rare1 1e-03-1e-05
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-4 & BHR_input$AF>1e-5)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_1e-04-1e-05.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r0_mis_h2=pheno_univariate_rare
  r0_mis=BHR_input_rare
  
  cat("\n###################################################  rare1 1e-03-1e-04##\n")
  ## rare1 1e-03-1e-05
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-3 & BHR_input$AF>1e-4)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_1e-03-1e-04.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r1_mis_h2=pheno_univariate_rare
  r1_mis=BHR_input_rare
  
  
  cat("\n###################################################  rare2 1e-02-1e-03##\n")
  ## rare2 1e-02-1e-03
  BHR_input_path<-paste("/home1/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/home1/Huashan1/wbs_data/software/BHR/data/ms_baseline_oe5.txt")
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<=1e-2 & BHR_input$AF>1e-3)
  pheno_univariate_rare <- BHR(mode = "univariate",
                               trait1_sumstats = BHR_input_rare,
                               annotations = list(baseline))
  total_h2<-pheno_univariate_rare$mixed_model$heritabilities[1,5]
  total_h2_se<-pheno_univariate_rare$mixed_model$heritabilities[2,5]
  intercept<-pheno_univariate_rare$qc$intercept
  intercept_se<-pheno_univariate_rare$qc$intercept_se
  lambda_gc<-pheno_univariate_rare$qc$lambda_gc
  result<-cbind(pheno_name,total_h2,total_h2_se,intercept,intercept_se,lambda_gc)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_1e-02-1e-03.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r2_mis_h2=pheno_univariate_rare
  r2_mis=BHR_input_rare
  
  cat("\n###################################################  rare aggregate##\n")
  ## rare aggregate
  
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(r0_mis,r1_mis,r2_mis),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_rare.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  r_mis_h2=pheno_univariate_rare
  
  cat("\n###################################################  ultra-rare + rare aggregate##\n")
  ## ultra-rare + rare aggregate
  
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(ur_mis,r0_mis,r1_mis,r2_mis),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_missense_urr.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  urr_mis_h2=pheno_univariate_rare
  
  
  ###############################################################################################################################################################
  cat("\n################################################### total##\n")
  ## ultra-rare
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(ur_lof,ur_mis),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_ultra_rare_aggregate.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  ag_ur_t_h2=pheno_univariate_rare
  
  
  cat("\n################################################### rare1 1e-3- 1e-5##\n")
  
  ## rare aggregate
  cat("\n################################################### rare aggregate##\n")
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(r0_lof,r1_lof,r2_lof,r0_mis,r1_mis,r2_mis),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_rare_aggregate.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  ag_r_t_h2=pheno_univariate_rare
  
  
  ## ultra-rare + rare aggregate
  cat("\n################################################### ultra-rare + rare aggregate##\n")
  pheno_univariate_rare <- BHR(mode = "aggregate",
                               trait_list=list(pheno_name),
                               ss_list = list(ur_lof,ur_mis,r0_lof,r1_lof,r2_lof,r0_mis,r1_mis,r2_mis),
                               annotations = list(baseline))
  
  total_h2<-pheno_univariate_rare$aggregated_mixed_model_h2
  total_h2_se<-pheno_univariate_rare$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,total_h2,total_h2_se,NA,NA,NA)
  write_path<-paste("/home1/BHR/BHR_output/",pheno_name,"_urr_t_aggregate.txt",sep = "")
  write.table(result,write_path,sep = "\t",row.names = F,quote = F)
  ag_urr_t_h2=pheno_univariate_rare
  
  
}

save.image("/home1/BHR.Rdata")

