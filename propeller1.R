# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle", "scater"))

remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE, 
                        dependencies = "Suggest")

library(speckle)
library(limma)
library(ggplot2)

props <- getTransformedProps(clusters = seurat.integrated$label.fine, 
                             sample = seurat.integrated$Patient)

p1 <- plotCellTypeProps(clusters = seurat.integrated$label.fine, sample = seurat.integrated$Patient) + theme(axis.text.x = element_text(angle = 90))+ ggtitle("Cell type proportions") + 
  theme(plot.title = element_text(size = 18, hjust = 0))
p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90))


p2 <- plotCellTypeProps(clusters = seurat.integrated$label.main, sample = seurat.integrated$Patient) + theme(axis.text.x = element_text(angle = 90))+ ggtitle("Broad cell type proportions") + 
  theme(plot.title = element_text(size = 18, hjust = 0))
p2 + theme_bw() + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle = 90))





png("./2nd/plotCellTypeProps_main.png",width=3000,height=2500,res=500)
dev.off()



counts <- table(seurat.integrated$label.fine, seurat.integrated$Patient)
baselineN <- rowSums(counts)
N <- sum(baselineN)
baselineprops <- baselineN/N

seurat.integrated$final_ct_fine <- factor(seurat.integrated$label.fine, levels=names(sort(baselineprops, decreasing = TRUE)))

counts <- table(seurat.integrated$final_ct_fine, seurat.integrated$Patient)
baselineN <- rowSums(counts)
N <- sum(baselineN)
baselineprops <- baselineN/N

props <- getTransformedProps(clusters = seurat.integrated$final_ct_fine, 
                             sample = seurat.integrated$Patient)

cols <- ggplotColors(nrow(props$Proportions))
m <- match(rownames(props$Proportions),levels(factor(seurat.integrated$label.fine)))

par(mfrow=c(1,1))
par(mar=c(7,5,2,2))
plot(jitter(props$Proportions[,1]), col = cols[m], pch=16, ylim=c(0,max(props$Proportions)),
     xaxt="n", xlab="", ylab="Cell type proportion", cex.lab=1.5, cex.axis=1.5)
for(i in 2:ncol(props$Proportions)){
  points(jitter(1:nrow(props$Proportions)),props$Proportions[,i], col = cols[m],
         pch=16)
}
axis(side=1, at=1:nrow(props$Proportions), las=2, 
     labels=rownames(props$Proportions))
title("Cell type proportions estimates for 12 individuals")


plotCellTypeMeanVar(counts)
plotCellTypePropsMeanVar(counts)


png("./2nd/plotCellTypePropsMeanVar_fine.png",width=3000,height=3000,res=500)
dev.off()

output.logit_fine <- propeller(clusters=seurat.integrated$label.fine, sample=seurat.integrated$Patient, group=seurat.integrated$Group, transform="logit")
output.logit_fine

output.logit_main <- propeller(clusters=seurat.integrated$label.main, sample=seurat.integrated$Patient, group=seurat.integrated$Group, transform="logit")
output.logit_main

write.csv(output.logit_main,file = "./2nd/output.logit_main.csv")


par(mfrow=c(2,3))
#1. Naive CD8 T cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Naive CD8 T cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Naive CD8 T cells", cex.main=1.5, adj=0)
text(1.5,0.3, labels = "Adj.Pval = 0.001036")

#2. Th17 cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Th17 cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Th17 cells", cex.main=1.5, adj=0)
text(1.5,0.07, labels = "Adj.Pval = 0.001893")

#3. Th1 cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Th1 cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Th1 cells", cex.main=1.5, adj=0)
text(1.5,0.06, labels = "Adj.Pval = 0.001036")

#4. Naive CD4 T cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Naive CD4 T cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Naive CD4 T cells", cex.main=1.5, adj=0)
text(1.5,0.25, labels = "Adj.Pval = 0.012398")

#5. Central memory CD8 T cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Central memory CD8 T cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Central memory CD8 T cells", cex.main=1.5, adj=0)
text(1.5,0.03, labels = "Adj.Pval = 0.019697")

#6. Progenitor cells

grp.covid <- rep(c("Healthy","OLP"), c(16,16))
stripchart(as.numeric(props$Proportions["Progenitor cells",])~grp.covid,
           vertical=TRUE, pch=16, method="jitter",
           col = c(4,"tomato",2),cex=1.5, 
           ylab="Proportions",cex.axis=1.25, cex.lab=1.5)
title("Progenitor cells", cex.main=1.5, adj=0)
text(1.5,0.003, labels = "Adj.Pval = 0.021801")

png("./2nd/stripchart_fine.png",width=3500,height=4500,res=500)
dev.off()        

