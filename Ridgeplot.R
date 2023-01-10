#GSE150327. Salivary Gland
setwd("/home/dongmin/R projects/Singlecell-Analysis/GSE150327/")


png("P1.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(p1.epithelium, features = RSC.marker)
dev.off()


png("P1.SG.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(p1.epithelium, features = SGprogenitor.marker)
dev.off()


png("P30.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(P30.epithelium, features = RSC.marker)
dev.off()


png("P30.SG.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(P30.epithelium, features = SGprogenitor.marker)
dev.off()

png("Adult.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(Adult.epithelium, features = RSC.marker)
dev.off()


png("Adult.SG.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(Adult.epithelium, features = SGprogenitor.marker)
dev.off()


png("smgP.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(smgP.epithelium, features = RSC.marker)
dev.off()

png("smgP.SG.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(smgP.epithelium, features = SGprogenitor.marker)
dev.off()
#GSE113466. Salivary Gland
setwd("/home/dongmin/R projects/Singlecell-Analysis/GSE113466/")

png("GSE113466.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(MP9SG, features = RSC.marker)
dev.off()

#GSE164403
setwd("/home/dongmin/R projects/Singlecell-Analysis/GSE164403")

png("tissue.RSC.RidgePlot.png", width =16000, height = 3500, res = 1200 )
RidgePlot(tissue, features = RSC.marker)
dev.off()

png("organoid.RSC.RidgePlot.png", width =14000, height = 7000, res = 1200 )
RidgePlot(organoid, features = RSC.marker,group.by = "CultureCondition")
dev.off()






p1.epithelium
p1.epithelium.markers <- FindAllMarkers(p1.epithelium, only.pos = T,logfc.threshold = 0.25, max.cells.per.ident = 500)
Undefined.mito<- subset(x = p1.epithelium.markers, subset = cluster == "Undefined mitotic cells")
Undefined.mito.filtered<-Undefined.mito[Undefined.mito$p_val_adj < 0.05 & Undefined.mito$avg_log2FC > 0.5,]

write.csv(Undefined.mito.filtered,file="/home/dongmin/R projects/Singlecell-Analysis/GSE150327/Undefined.mito.filtered.csv")


