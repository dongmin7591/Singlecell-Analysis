
seurat.integrated<-readRDS(file = "./Final.seurat.integrated.rds")

# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle", "scater"))

remotes::install_github("Oshlack/speckle", build_vignettes = TRUE, 
                        dependencies = "Suggest")

library(devtools)
devtools::install_github("phipsonlab/speckle")

browseVignettes("speckle")

library(speckle)
library(limma)
library(ggplot2)

# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
cell_info <- speckle_example_data()
head(cell_info)

# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = cell_info$clusters, sample = cell_info$samples, 
          group = cell_info$group)

# Plot cell type proportions
plotCellTypeProps(clusters=cell_info$clusters, sample=cell_info$samples)

install.packages("gt")
library(speckle)
library(limma)
library(edgeR)
library(pheatmap)
library(gt)

source("./code/convertData.R")



props <- getTransformedProps(clusters = seurat.integrated$singler_labels, 
                             sample = seurat.integrated$Patient)



png("./png/Cell.type.proportions.png",width=5000,height=2500,res=500)

p1 <- plotCellTypeProps(clusters = seurat.integrated$singler_labels, sample = seurat.integrated$Patient) + theme(axis.text.x = element_text(angle = 45))+ ggtitle("Cell type proportions") + 
  theme(plot.title = element_text(size = 20, hjust = 0))

p1 + theme_linedraw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(size=11,angle = 90),axis.text.y = element_text(size= 11), axis.title.y=element_text(size= 12))

png("./png/plotCellTypePropsMeanVar.png",width=5000,height=5000,res=500)

plotCellTypePropsMeanVar(counts)
dev.off()

plotCellTypePropsMeanVar(counts)

output.logit <- propeller(clusters=seurat.integrated$singler_labels, sample=seurat.integrated$Patient, group=seurat.integrated$Group, transform="logit")
output.logit

des.dose <- model.matrix(~grp.covid)
des.dose

fit.plot <- lmFit(props.covid$TransformedProps,des.dose)
fit.plot <- eBayes(fit.plot, robust=TRUE)

output.asin <- propeller(clusters=seurat.integrated$singler_labels, sample=seurat.integrated$Patient, group=seurat.integrated$Group, transform="asin")
output.asin

props.covid <- getTransformedProps(clusters=seurat.integrated$singler_labels, sample=seurat.integrated$Patient,
                                   transform="logit")
  
write.csv(output.logit, file = "./propeller.csv")

par(mfrow=c(1,1))
topTable(fit.plot,coef=2)




par(mfrow=c(1,3))
#1. Th17 cells
grp.covid <- rep(c("Healthy","OLP"), c(16, 16))
stripchart(as.numeric(props.covid$Proportions["Th17 cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)+
title("Th17 cells", cex.main=1.5, adj=0)+
text(1.5,0.06, labels = "Adj.Pval = 7.80E-04",cex=1.2)
#2. Central memory CD8 T cells
grp.covid <- rep(c("Healthy","OLP"), c(16, 16))
stripchart(as.numeric(props.covid$Proportions["Central memory CD8 T cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)+
title("Central memory CD8 T cells", cex.main=1.5, adj=0)+
text(1.5,0.03, labels = "Adj.Pval = 1.57E-02",cex=1.2)

#3. Th1 cells
grp.covid <- rep(c("Healthy","OLP"), c(16, 16))
stripchart(as.numeric(props.covid$Proportions["Th1 cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)+
title("Th1 cells", cex.main=1.5, adj=0)+
text(1.5,0.06, labels = "Adj.Pval = 4.91E-02",cex=1.2)
     # Ploting


png("./png/stripchart.png",width=5000,height=4000,res=500)

dev.off()        

