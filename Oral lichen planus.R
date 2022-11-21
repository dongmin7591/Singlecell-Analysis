options(future.globals.maxSize=1024^2 * 8000)
plan("multisession",workers=4)
plan("sequential")
#! /usr/bin/env Rscript
getwd()
# Integrate single-cell RNA-seq dataset
#saveRDS(seurat.integrated, file = "../output/seurat.integrated3k_final.rds")

# script to integrate scRNA-Seq datasets to correct for batch effects
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

install.packages('Seurat',force = TRUE)
install.packages('httpuv')

# load libraries
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(sctransform)
library(dittoSeq)
# get data location

dirs <- list.dirs(path = 'Data/', recursive = F, full.names = F)

for (x in dirs) {
  name <- gsub('_filtered_feature_bc_matrix','',x)
  
  cts <- ReadMtx(mtx = paste0('./Data/',x,'/matrix.mtx.gz'),
                 features = paste0('./Data/', x,'/features.tsv.gz'),
                 cells = paste0('./Data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name,CreateSeuratObject(counts = cts))
}


# merge datasets

ids <- c("OLP_OLP1", "OLP_OLP2", "OLP_OLP3", "OLP_OLP4", "OLP_OLP5", "OLP_OLP6", "OLP_OLP7", "OLP_OLP8", "OLP_OLP9", "OLP_OLP10", "OLP_OLP11", "OLP_OLP12", "OLP_OLP13",
         "OLP_OLP14", "OLP_OLP15", "OLP_OLP16", 
         "Healthy_Healthy1", "Healthy_Healthy2", "Healthy_Healthy3", "Healthy_Healthy4", "Healthy_Healthy5", "Healthy_Healthy6", "Healthy_Healthy7", 
         "Healthy_Healthy8", "Healthy_Healthy9", "Healthy_Healthy10", "Healthy_Healthy11","Healthy_Healthy12", "Healthy_Healthy13", "Healthy_Healthy14",
         "Healthy_Healthy15", "Healthy_Healthy16")

merged_seurat<- merge(OLP1, y=c(OLP2, OLP3, OLP4, OLP5, OLP6, OLP7, OLP8, OLP9, OLP10, OLP11, OLP12, OLP13, OLP14, OLP15, OLP16,
                                Healthy1, Healthy2, Healthy3, Healthy4, Healthy5, Healthy6, Healthy7, Healthy8, Healthy9, Healthy10, Healthy11,
                                Healthy12, Healthy13, Healthy14, Healthy15, Healthy16),
                      add.cell.ids= ids,
                      project = 'OLP')

merged_seurat



##################################
##################################
##################################

# nFeature_RNA : The number of genes detected in each cell.
# nCount_RNA : The total number of molecules detected within a cell. or number of UMIs per cell.
# nFeature_RNA for a cell indicates that it may be dead/dying or an empty droplet. 
# High nCount_RNA and/or nFeature_RNA indicates that the "cell" may in fact be a doublet 
# (or multiplet). In combination with %mitochondrial reads, removing outliers from these 
# groups removes most doublets/dead cells/empty droplets, hence why filtering is a common 
# pre-processing step.

# QC & filtering


View(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample<- rownames(merged_seurat@meta.data)



# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Group','Patient','Barcode'),
                                    sep = '_')

#saveRDS(merged_seurat, file = "./OLP.rds")
# calculate mitochondrial percentage

View(merged_seurat@meta.data)

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3, split.by = 'Patient')

# before quality control
VlnPlot(merged_seurat, features =  "mitoPercent", group.by = 'Patient',fill.by = "feature",pt.size=0)
#
png("./png/beforeQCplot.png",width=9000,height=3000,res=500)
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), group.by = 'Patient',fill.by = "feature",pt.size=0)
dev.off()
#
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

png("./png/beforeQC_scatter plot.png",width=5000,height=3000,res=500)
plot1 + plot2
dev.off()

#explore QC



#Filtering

View(merged_seurat)

merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 500 & 
                                   nFeature_RNA < 6000 &
                                   mitoPercent < 10)


merged_seurat_filtered # 268917 >>> 234487   

# After quality control

png("./png/afterQC.png",width=9000,height=3000,res=500)
VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), group.by = 'Patient',fill.by = "feature",pt.size=0)
dev.off()
#
plot3 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot4 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4

png("./png/afterQC_scatter plot.png",width=5000,height=3000,res=500)
plot3 + plot4
dev.off()




# The NormalizeData step is basically just ensuring expression values across cells are on a comparable scale. 
# By default, it will divide counts for each gene by the total counts in the cell, 
# multiply that value for each gene by the scale.factor (10,000 by default), and then natural log-transform them.


# perform standard workflow steps to figure out if we see any batch effects

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <-  NormalizeData(object = merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:15)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:15)



DimPlot(merged_seurat_filtered, reduction = "umap", group.by='Patient')

#
png("./OLP/VizDimLoadings.png",width=3000,height=3000,res=500)
VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")
dev.off()

DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)

png("./OLP/DimHeatmap1.png",width=3000,height=3000,res=500)
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)
dev.off()


DimHeatmap(merged_seurat_filtered, dims = 1:15, cells = 500, balanced = TRUE)

png("./OLP/DimHeatmap15.png",width=5000,height=5000,res=500)
DimHeatmap(merged_seurat_filtered, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


#Determine the ‘dimensionality’ of the dataset
# Computes the nearest neighbors for a given dataset. Can also optionally (via ), 
# construct a shared nearest neighbor graph by calculating the neighborhood overlap (Jaccard index) 
# between every cell and its nearest
#merged_seurat_filtered<- FindNeighbors(object = merged_seurat_filtered, dims=1:10)
#merged_seurat_filtered<- FindClusters(object = merged_seurat_filtered, resolution = 0.5)


# Since kNN typically uses euclidian distance to find k nearest points from any given point, 
# using normalized features may select a different set of k neighbors than the ones chosen 
# when unnormalized features were used, hence the difference in accuracy.



# Doublet Finder
#pK Identification (no ground-truth)
sweep.res.list_nsclc <- paramSweep_v3(merged_seurat_filtered, PCs = 1:20, sct = FALSE )
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)


ggplot(bcmvn_nsclc, aes(pK, BCmetric, group = 1))+
  geom_point() +
  geom_line()

pK<- bcmvn_nsclc %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)

pK <- as.numeric(as.character(pK[[1]]))



# Homotypic doublet proportion estimate

annotations <- merged_seurat_filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*nrow(merged_seurat_filtered@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# run doubletFinder
merged_seurat_filtered <- doubletFinder_v3(merged_seurat_filtered,
                                           PCs = 1:20,
                                           pN = 0.25,
                                           pK = pK,
                                           nExp = nExp_poi.adj,
                                           reuse.pANN = FALSE, 
                                           sct = FALSE)

# visualize doublets
png("./OLP/multiplets.png",width=3000,height=3000,res=500)
DimPlot(merged_seurat_filtered_nofiltered, reduction = 'umap', group.by = "DF.classifications_0.25_0.19_4283")

dev.off()


merged_seurat_filtered_nofiltered<-merged_seurat_filtered

merged_seurat_filtered <- subset(merged_seurat_filtered, subset = DF.classifications_0.25_0.19_4283 == "Singlet")

# Plot

png("./OLP/UMAP_nofilter.png",width=3000,height=3000,res=500)
DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
dev.off()

p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')

#p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type', cols = c('red','green','blue'))
grid.arrange(p1,p2, ncol = 2, nrow = 2)


#1. harmony
# run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'Group', plot_convergence = FALSE,verbose=TRUE)


ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
ifnb.harmony <- JackStraw(ifnb.harmony, num.replicate = 100)
ifnb.harmony <- ScoreJackStraw(ifnb.harmony, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)

# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'Group',raster = FALSE)

before|after

#2. perform integration to correct for batch effects -----

obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Group')
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features

features <- SelectIntegrationFeatures(object.list = obj.list)


# find integration anchors (Canonical correlation analysis (CCA)), 
# reduction : CCA, Reciprocal PCA (rpca), Reciprocal LSI (rlsi)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features,
                                  dims = 1:15)


saveRDS(seurat.integrated,file = "./2nd/seurat.integrated.rds")
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)
?IntegrateData
# Scale data, run PCA and UMAP and visualize integrated data

seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)

saveRDS(seurat.integrated, file = './seurat.integrated.rds') 
#new
VizDimLoadings(seurat.integrated, dims = 1:2, reduction = "pca")
DimPlot(seurat.integrated, reduction = "pca", group.by = "Patient")




seurat.integrated <- JackStraw(seurat.integrated, num.replicate = 100)
seurat.integrated <- ScoreJackStraw(seurat.integrated, dims = 1:20)
JackStrawPlot(seurat.integrated, dims = 1:20,ymax = 0.8)
ElbowPlot(seurat.integrated)

png("./2nd/ElbowPlot.png",width=3000,height=3000,res=500)
dev.off()

seurat.integrated<- FindNeighbors(object = seurat.integrated, dims=1:15)
seurat.integrated<- FindClusters(object = seurat.integrated, resolution = 0.5)


seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:15)

###


p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient',raster = FALSE)

p2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Group',raster = FALSE)


p1|p2
png("./2nd//Patient_Group_UMAP.png",width=6000,height=3000,res=500)
grid.arrange(p1,p2, ncol = 2, nrow = 1)
dev.off()



# Visualization
# cell clustring
p3 <- DimPlot(seurat.integrated, reduction = "umap", label = TRUE)
p2 + p3




#그림 저장
png("./OLP/cell clustring.png",width=6000,height=3000,res=500)
grid.arrange(p2,p3, ncol = 2, nrow = 1)
dev.off()

png("./OLP/IntegratedDimHeatmap15.png",width=5000,height=5000,res=500)
DimHeatmap(seurat.integrated, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


saveRDS(seurat.integrated, file = "../scRNA/OLP/OLP.rds")



# SingleR: Label single-cell RNA data automatically

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")

BiocManager::install("celldex")



library(celldex)
library(SingleR)

data <- seurat.integrated

ref <- celldex::MonacoImmuneData()

results <- SingleR(test = as.SingleCellExperiment(data), ref = ref, labels =
                     ref$label.fine,assay.type.test=1)
table(results$labels)

plotScoreHeatmap(results)
#########################
# Annotation diagnostics
#########################
plotDeltaDistribution(results, ncol = 3)

data$pruned.fine<- results$pruned.labels
data<-subset(data,pruned.fine!="NA")
summary(is.na(data$pruned.fine))

table(data$pruned.fine)

data$label.fine<- results$labels


data$pruned.labels<- results$pruned.labels
data<-subset(data,pruned.labels!="TRUE")

DimPlot(data, reduction = 'umap', group.by = 'label.main', label = F,raster = FALSE )

# Basophils,Neutrophils 제거

data<-subset(data,label.main!="Basophils" & label.main!= "Neutrophils"& label.fine!="Low-density basophils" & label.fine!= "Low-density neutrophils")
table(test@meta.data$label.main)
# label.fine
data<-subset(data,label.fine!="Low-density basophils" & label.fine!= "Low-density neutrophils")

ref$label.main

png("./2nd//label.fine_UMAP.png",width=6000,height=3000,res=500)
DimPlot(data, reduction = 'umap', group.by = 'label.fine', label = F,label.size = 5,raster = FALSE) +ggtitle("Fine-grained annotation") 
dev.off()

png("./2nd//label.main_UMAP.png",width=4000,height=3000,res=500)
DimPlot(data, reduction = 'umap', group.by = 'label.main', label = F,label.size = 4,raster = FALSE) +ggtitle("Broad annotation") 
dev.off()

plotScoreHeatmap(results)
plotDeltaDistribution(results)

MainData@meta.data$seurat_clusters



dittoBarPlot(seurat.integrated, "label.main", group.by = "Group")
dittoBarPlot(seurat.integrated, "label.fine", group.by = "Group")

png("./2nd/dittoBarPlot_fine.png",width=4000,height=3000,res=500)
dev.off()


####################
## Subcluster is now saved in metadata. In this example it's location is "MainData@meta.data$singler_labels".
## Now use SetIdent to save your new cluster assignment to the main Ident in your object.


MainData<-SetIdent(MainData, value = MainData@meta.data$singler_labels)

DimPlot(MainData, reduction = "umap", label = TRUE, label.size = 6)

singler.markers.main <- FindAllMarkers(MainData, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

# findAll markers -----------------
BiocManager::install("SeuratData")

library(SeuratData)
library(ggplot2)
library(patchwork)


singler.markers.main<- FindAllMarkers(MainData,
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1,
                                      only.pos = TRUE,
)




B.cells.markers<- FindMarkers(MainData, ident.1 = "B cells", min.pct = 0.25)



# let's visualize top features
FeaturePlot(MainData, features = c('HLA-DRA','CD79A','CD74','IGLL5'), min.cutoff = 'q10')

RidgePlot(MainData, features = c('HLA-DRA','CD79A','CD74','IGLL5'), ncol = 2)

VlnPlot(MainData, features = c('HLA-DRA','CD79A','CD74','IGLL5'))

DoHeatmap(subset(MainData, downsample = 100), features = c('HLA-DRA','CD79A','CD74','IGLL5'), size = 3)
# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))



View(B.cells.markers)

# findMarkers between conditions ---------------------
MainData$celltype.cnd <- paste0(MainData$seurat_annotations,'_', MainData$stim)
View(MainData@meta.data)
Idents(MainData) <- MainData$celltype.cnd

DimPlot(MainData, reduction = 'umap', label = TRUE)
MainData@meta.data$singler_labels
# find markers
B_CD8+Tcells <- FindMarkers(MainData, ident.1 = 'B cells', ident.2 = 'CD8+ T cells')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)

FeaturePlot(MainData, features = c('FCGR3A', 'AIF1', 'IFIT1'), min.cutoff = 'q10')



######################################
######################################
######################################

# Cell chat


devtools::install_github("sqjin/CellChat")
install.packages("igraph")

install.packages("devtools")
require(devtools)
install_version("igraph", version = "1.3.3", repos = "http://cran.us.r-project.org")

library(igraph)

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# expression data

df<-MainData

DefaultAssay(df) <- 'RNA'

data.input <- GetAssayData(df, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(df)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

#Preprocessing the expression data for cell-cell communication analysis

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# 
table(cellchat@idents)

#Inference of cell-cell communication network

cellchat <- computeCommunProb(cellchat)

df.net <- subsetCommunication(cellchat)


# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))



par(mfrow = c(1,2), xpd=T) 
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#df.net 에서 pathway 확인
table(df.net$pathway_name)

pathways.show <- c("MIF") 

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot

png("./OLP/MIFsig.png",width=2000,height=2000,res=500)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap

png("./OLP/MIFsig_heatmap.png",width=2500,height=2500,res=500)
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()


png("./OLP/MIFsig_net.png",width=2000,height=1500,res=500)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")


plotGeneExpression(cellchat, signaling = "MIF")
#> Registered S3 method overwritten by 'spatstat.geom':
#>   method     from
#>   print.boxx cli
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")

netAnalysis_dot(cellchat, pattern = "outgoing")
