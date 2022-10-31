library(Seurat)
library(tidyverse)
install.packages("BiocManager") # Needed to install all Bioconductor packages
BiocManager::install("MAST")

library(MAST)

seurat.integrated@meta.data
seurat.integrated$celltype.cnd <- paste0(seurat.integrated$singler_labels,'_', seurat.integrated$Group)



View(seurat.integrated@meta.data)
Idents(seurat.integrated) <- seurat.integrated$celltype.cnd

table(seurat.integrated@meta.data$celltype.cnd)
seurat.integrated@assays
#1. Th17 cells

Th17.cells.Markers <- FindMarkers(seurat.integrated, ident.1 = 'Th17 cells_Healthy', ident.2 = 'Th17 cells_OLP',test.use = "MAST")

Th17.cells.Healthy<-Th17.cells.Markers[Th17.cells.Markers$p_val_adj < 0.05 & Th17.cells.Markers$avg_log2FC > 0.5,]
Th17.cells.OLP<-Th17.cells.Markers[Th17.cells.Markers$p_val_adj < 0.05 & Th17.cells.Markers$avg_log2FC < -0.5,]


hh<-sort(row.names(Th17.cells.Healthy)) 
p1<-FeaturePlot(seurat.integrated, row.names(Th17.cells.OLP), split.by = 'Group', min.cutoff = 'q10')
p2<-FeaturePlot(seurat.integrated, features = hh[5:8], split.by = 'Group', min.cutoff = 'q10')

grid.arrange(p1, p2,)

p1|p2


write.csv(Th17.cells.OLP,file = "./Th17.cells.OLP.csv")

#2. Central memory CD8 T cells

Central.memory.CD8.T.cells.Markers <- FindMarkers(seurat.integrated, ident.1 = 'Central memory CD8 T cells_Healthy', ident.2 = 'Central memory CD8 T cells_OLP',test.use = "DESeq2")

Central.memory.CD8.T.cells.Healthy<-Central.memory.CD8.T.cells.Markers[Central.memory.CD8.T.cells.Markers$p_val_adj < 0.05 & Central.memory.CD8.T.cells.Markers$avg_log2FC > 0.5,]
Central.memory.CD8.T.cells.OLP<-Central.memory.CD8.T.cells.Markers[Central.memory.CD8.T.cells.Markers$p_val_adj < 0.05 & Central.memory.CD8.T.cells.Markers$avg_log2FC < -0.5,]

write.csv(Central.memory.CD8.T.cells.OLP,file = "./Central.memory.CD8.T.cells.OLP.csv")

#3. Th1 cells


Th1.cells.Markers <- FindMarkers(seurat.integrated, ident.1 = 'Th1 cells_Healthy', ident.2 = 'Th1 cells_OLP',test.use = "MAST")

Th1.cells.Healthy<-Th1.cells.Markers[Th1.cells.Markers$p_val_adj < 0.05 & Th1.cells.Markers$avg_log2FC > 0.5,]
Th1.cells.OLP<-Th1.cells.Markers[Th1.cells.Markers$p_val_adj < 0.05 & Th1.cells.Markers$avg_log2FC < -0.5,]


write.csv(Th1.cells.OLP,file = "./Th1.cells.OLP.csv")


Group_DEGs <- FindMarkers(seurat.integrated, ident.1 = 'Healthy', ident.2 = 'OLP',test.use = "MAST")

Group_DEGs.Healthy<-Group_DEGs[Group_DEGs$p_val_adj < 0.05 & Group_DEGs$avg_log2FC > 0.5,]
Group_DEGs.OLP<-Group_DEGs[Group_DEGs$p_val_adj < 0.05 & Group_DEGs$avg_log2FC < -0.5,]


