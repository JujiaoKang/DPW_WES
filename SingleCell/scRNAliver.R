# source activate r41
# conda install -c conda-forge r-seurat
# conda install -c conda-forge r-vctrs
# 
# conda install -c conda-forge r-crayon
# conda install -c conda-forge r-pkgconfig
# conda update r-vctrs
# conda install -c conda-forge r-rcolorbrewer 
# 
# conda install -c "conda-forge/label/cf202003" r-vctrs
# conda install -c "conda-forge/label/cf202003" r-data.table
# 
# conda update r-vctrs
# conda install -c r r-data.table
# conda install -c conda-forge r-patchwork
# which R

##


rm(list=ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggsci)
data <- fread(file = "/home1/GSE115469_Data.csv",header =T,data.table=F)

head(data[,1:10])
row.names(data)=data$V1
data=data[,-1]
dim(data)
ss=row.names(data)
grep("ADH1",ss)
ss[grep("ADH1",ss)]
rm(ss)
data1 <- CreateSeuratObject(counts = data, project = "t1", min.cells = 3, min.features = 200)
data1
#An object of class Seurat 
#19984 features across 8439 samples within 1 assay 
#Active assay: RNA (19984 features, 0 variable features)

data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^MT-")


data1@meta.data$orig.ident = "a"
png("/home1/GSE115469_Vlnplot.png")
VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,group.by = "orig.ident")
dev.off()


data1 <- subset(data1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &  percent.mt < 7.5)
data1 <- NormalizeData(data1, normalization.method = "LogNormalize", scale.factor = 10000)
data1
#An object of class Seurat 
#19984 features across 8420 samples within 1 assay 
#Active assay: RNA (19984 features, 0 variable features)

data1 <- FindVariableFeatures(data1, selection.method = "vst", nfeatures = 5000) %>% ScaleData() %>% RunPCA(verbose=FALSE)
data1 <- FindNeighbors(data1, dims = 1:40)
data1 <- FindClusters(data1, resolution = 0.8)

data1 <- RunTSNE(data1, dims = 1:25)

new.cluster.ids <- c("Hepatocyte","Hepatocyte","Hepatocyte","alpha-beta_T_Cells/gamma-delta_T_Cells",
                     "Hepatocyte","alpha-beta_T_Cells/gamma-delta_T_Cells","Plasma_Cells","Macrophage","alpha-beta_T_Cells/gamma-delta_T_Cells",
                     "Central_venous_LSECs","Periportal_LSECs","Macrophage","NK-like_Cells","Macrophage","Macrophage","Portal_endothelial_Cells",
                     "NK-like_Cells","Hepatocyte","Mature_B_Cells","alpha-beta_T_Cells/gamma-delta_T_Cells","Hepatocyte","Erythroid_Cells","Cholangiocytes",
                     "Hepatocyte","Hepatic_Stellate_Cells","Cholangiocytes")
                     
names(new.cluster.ids) <- levels(data1)
data2 <- RenameIdents(data1, new.cluster.ids)
DimPlot(data2, reduction = "tsne", pt.size = 0.5) + NoLegend()
mycol <- c(pal_d3()(7),pal_aaas()(7),pal_uchicago()(7),pal_jama()(7),pal_nejm()(7))
p1<-DimPlot(data2, reduction = "tsne", cols = mycol, pt.size = 0.3)
p2<-FeaturePlot(data2,reduction = "tsne", features = c("SERPINC1","SRSF6","PHPT1","STAB2","PROS1","PROC"),col=c("#D9DADB","#6A5C9C","#274A88"),pt.size = 0.3)


pdf("/home1/GSE115469_tsne.pdf",width=8,height=5)
p1
dev.off()
pdf("/home1/GSE115469_feature.pdf",width=10,height=5)
p2
dev.off()

