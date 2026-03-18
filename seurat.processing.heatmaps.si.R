library(Seurat)
library(sctransform)
library(dplyr)
library(stringr)
library(ggplot2)

df = read.csv("cc.genes.df.csv")
df2 = read.delim("pbio.2004050.csv", sep = ",")
esr.genes = word(df2$Annotation, 4)



pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/seurat_in")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc@meta.data$nCount_RNA, base = 2))

#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

#expr.mfalpha2<- FetchData(object = pbmc, vars = c("MF(ALPHA)2"))
dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
#dim(pbmc)


pbmc= SCTransform(pbmc)
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
#unlike above, do not use just cc genes for pca
pbmc= RunPCA(pbmc, npcs = 12)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca")
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.25  )
dim(pbmc)


FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)


#current
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

dim(top10)


#new
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

dim(top10)

DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1")) 





#########################-----------------BYxRM SP t10

pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t10//filtered_feature_bc_matrix/seurat_in/")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc@meta.data$nCount_RNA, base = 2))

#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

#expr.mfalpha2<- FetchData(object = pbmc, vars = c("MF(ALPHA)2"))
dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
#dim(pbmc)


pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
dim(pbmc)

pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.2 )


dim(pbmc)



FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))


#mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 25) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1"))





###########
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj > 1) %>%
  slice_head(n = 25) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1"))


#output markers and deal with cluster 0 
df = pbmc.markers
for(i in 1:length(unique(df$cluster))){
  
  #assign(paste("cluster", i, sep = "."), subset(pbmc.markers, pbmc.markers$cluster == i))
  tdf <- subset(df, df$cluster== i-1)
  label <- paste0("3051.spt10.cluster",i-1,"ids.txt", sep="")
  write.table(tdf$gene, quote=FALSE, file=label, col.names=FALSE, row.names=FALSE)
  
}

write.table(rownames(pbmc@assays$SCT@scale.data), quote=FALSE, file="3051.spt10.background.genes.txt", col.names=FALSE, row.names=FALSE)







###############################------------BYxRM salt t0

pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t0/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)


expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
expr.dse2 <- FetchData(object = pbmc, vars = c("DSE2"))
#expr.mfalpha2<- FetchData(object = pbmc, vars = c("MFA2"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < 1)]
dim(pbmc)

hist(log(pbmc@meta.data$nCount_RNA, base = 2))

pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
#pca using cc genes
pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12) 
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3 )
dim(pbmc)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(top10$gene, cc_markers))



#assign labels to clusters, combining what we think are similar populations/clusters based on markers  
naclt03051.cc = c("G2_M + mating", "S", "G1_S", "S", "G1","G2_M w/o mating",  "M_G1 + mating", "M/G1 w/o mating")
names(naclt03051.cc) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, naclt03051.cc)




###############################################################################----------BYxRM salt t30
##########################3051 nacl t30
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t30/outs//filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc@meta.data$nCount_RNA, base = 2))


#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

#expr.mfalpha2<- FetchData(object = pbmc, vars = c(""))
dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
#dim(pbmc)


pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12) 
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.43 )

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(top10$gene, cc_markers))

#assign state labels to clusters based on marker expression 
naclt303051.cc = c("G1_S", "G1", "G2_M_no_mating",  "G2_M_mating" , "M_G1" ,"S", "G2_M_mating","S")
names(naclt303051.cc) <- levels(pbmc)

pbmc <- RenameIdents(pbmc, naclt303051.cc)






#######################------------3004 nacl t0
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t0_2/NaCl_0.7M_3004_rep1_t0/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)

#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
#expr.mfalpha2<- FetchData(object = pbmc, vars = c("MF(ALPHA)2"))
dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
#dim(pbmc)


pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
dim(pbmc)
#pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )


dim(pbmc)

hist(log(pbmc@meta.data$nCount_RNA, base = 2))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)

#mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(top10$gene, cc_markers))

#removal of cluster of indistinct cells that overlap with 'unclear' ones 
dim(pbmc)
cells_to_remove <- WhichCells(pbmc, idents = c(4))
pbmc <- subset(pbmc, cells = setdiff(Cells(pbmc), cells_to_remove))
dim(pbmc)

#repeat sctransform using subset object
pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3)
dim(pbmc)

FeaturePlot(pbmc, features = cc_markers)

#find markers and make another heatmap using final (cc) clusters
#pbmc = readRDS("3004.nacl.t0.alphaandcluster4filtered.RDS")



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(top10$gene, cc_markers))

#make state assignments using new clusters 
naclt03004.cc = c("G2_M + mating", "S", "M_G1 + Mating", "G2_M w/o mating", "S", "G1_S", "S", "M_G1 w/o Mating" )
names(naclt03004.cc) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, naclt03004.cc)




######################3004 nacl t30 

pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t30_2/NaCl_0.7M_3004_rep1_t30//filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc@meta.data$nCount_RNA, base = 2))

#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)



pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5,  )
dim(pbmc)


FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))


#mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(c(top10$gene, cc_markers))) 


###############################cluster 6 has the 'unclear' cells in this case 
#removal of cluster of indistinct cells that overlap with 'unclear' ones 
dim(pbmc)
cells_to_remove <- WhichCells(pbmc, idents = c(6))
pbmc <- subset(pbmc, cells = setdiff(Cells(pbmc), cells_to_remove))
dim(pbmc)


#repeat sctransform using subset object
pbmc= SCTransform(pbmc) 
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
pbmc= RunPCA(pbmc, npcs = 12)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.43,  )
dim(pbmc)



pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = c(c(top10$gene, cc_markers))) 



FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)



#here and in 3051 i am adding a high rp + lowish s cluster to s but it could be left out 
#nacl.t30.3004.cc = c("g1/S", "Stress", "M/G1", "G2/M", "S", "G2/M", "G2/M", "S", "G2/M", "M/G1", "G1")
nacl.t30.3004.cc = c("G1_S", "Stress", "M_G1 + mating", "G2_M + mating", "S", "G2_M w/o mating", "G2_M + mating", "S", "Stress", "M_G1 w/o mating", "G1")




names(nacl.t30.3004.cc) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, nacl.t30.3004.cc)
