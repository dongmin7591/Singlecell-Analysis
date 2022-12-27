

Idents(seurat.all)<-seurat.all$IntercellDB

install.packages('enrichR')
library(enrichR)

dbs <- c("GO_Biological_Process_2021", "GO_Cellular_Component_2021", "GO_Molecular_Function_2021", "KEGG_2021_Human")
DEenrich<-list()

bp<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Biological_Process_2021"),
  num.pathway = 5,
  return.gene.list = F)

cc<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Cellular_Component_2021"),
  num.pathway = 5,
  return.gene.list = F)
mf<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Molecular_Function_2021"),
  num.pathway = 5,
  return.gene.list = F)

kegg<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("KEGG_2021_Human"),
  num.pathway = 5,
  return.gene.list = F)


png("./2nd/IntercellDB/healthyVSOLP/DEenriched/bp_Th17 cells.png",width=7000,height=1500,res=500)
bp
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/cc_Th17 cells.png",width=7000,height=1500,res=500)
cc
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/mf_Th17 cells.png",width=7000,height=1500,res=500)
mf
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/kegg_Th17 cells.png",width=7000,height=1500,res=500)
kegg
dev.off()


bp_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Biological_Process_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)

cc_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Cellular_Component_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)
mf_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Molecular_Function_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)

kegg_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th17 cells_OLP",
  ident.2 = "Th17 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("KEGG_2021_Human"),
  num.pathway = 5,
  return.gene.list = TRUE)

bp_gene$pos$group<-"upregulated"
bp_gene$neg$group<-"downregulated"

cc_gene$pos$group<-"upregulated"
cc_gene$neg$group<-"downregulated"

mf_gene$pos$group<-"upregulated"
mf_gene$neg$group<-"downregulated"

kegg_gene$pos$group<-"upregulated"
kegg_gene$neg$group<-"downregulated"

bp_gene_rbind<-rbind(bp_gene$pos,bp_gene$neg)
cc_gene_rbind<-rbind(cc_gene$pos,cc_gene$neg)
mf_gene_rbind<-rbind(mf_gene$pos,mf_gene$neg)
kegg_gene_rbind<-rbind(kegg_gene$pos,kegg_gene$neg)

write.csv(bp_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th17 cells_bp.csv")
write.csv(cc_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th17 cells_cc.csv")
write.csv(mf_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th17 cells_mf.csv")
write.csv(kegg_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th17 cells_kegg.csv")


bp<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Biological_Process_2021"),
  num.pathway = 5,
  return.gene.list = F)

cc<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Cellular_Component_2021"),
  num.pathway = 5,
  return.gene.list = F)
mf<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Molecular_Function_2021"),
  num.pathway = 5,
  return.gene.list = F)

kegg<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("KEGG_2021_Human"),
  num.pathway = 5,
  return.gene.list = F)


png("./2nd/IntercellDB/healthyVSOLP/DEenriched/bp_Th1 cells.png",width=7000,height=1500,res=500)
bp
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/cc_Th1 cells.png",width=7000,height=1500,res=500)
cc
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/mf_Th1 cells.png",width=7000,height=1500,res=500)
mf
dev.off()



png("./2nd/IntercellDB/healthyVSOLP/DEenriched/kegg_Th1 cells.png",width=7000,height=1500,res=500)
kegg
dev.off()


bp_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Biological_Process_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)

cc_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Cellular_Component_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)
mf_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("GO_Molecular_Function_2021"),
  num.pathway = 5,
  return.gene.list = TRUE)

kegg_gene<-DEenrichRPlot(
  seurat.all,
  ident.1 = "Th1 cells_OLP",
  ident.2 = "Th1 cells_Healthy",
  balanced = TRUE,
  logfc.threshold = 0.3,
  assay = "RNA" ,
  max.genes=100,
  test.use = "MAST",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = c("KEGG_2021_Human"),
  num.pathway = 5,
  return.gene.list = TRUE)

bp_gene$pos$group<-"upregulated"
bp_gene$neg$group<-"downregulated"

cc_gene$pos$group<-"upregulated"
cc_gene$neg$group<-"downregulated"

mf_gene$pos$group<-"upregulated"
mf_gene$neg$group<-"downregulated"

kegg_gene$pos$group<-"upregulated"
kegg_gene$neg$group<-"downregulated"

bp_gene_rbind<-rbind(bp_gene$pos,bp_gene$neg)
cc_gene_rbind<-rbind(cc_gene$pos,cc_gene$neg)
mf_gene_rbind<-rbind(mf_gene$pos,mf_gene$neg)
kegg_gene_rbind<-rbind(kegg_gene$pos,kegg_gene$neg)

write.csv(bp_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th1 cells_bp.csv")
write.csv(cc_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th1 cells_cc.csv")
write.csv(mf_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th1 cells_mf.csv")
write.csv(kegg_gene_rbind,file = "./2nd/IntercellDB/healthyVSOLP/DEenriched/Th1 cells_kegg.csv")

