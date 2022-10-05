# memory partitioning
#reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]
sce.all.list <- SplitObject(seu, split.by = "Patient")
sce.all.list
DefaultAssay(seu) <- "RNA"
library(DoubletFinder)
#sce.all.filt=sce.all.list[[1]]
phe_lt <- lapply(names(sce.all.list), function(x){
  sce.all.filt=sce.all.list[[x]]
  sce.all.filt = FindVariableFeatures(sce.all.filt)
  sce.all.filt = ScaleData(sce.all.filt)
  sce.all.filt = RunPCA(sce.all.filt, npcs = 20)
  sce.all.filt = RunTSNE(sce.all.filt, npcs = 20)
  sce.all.filt = RunUMAP(sce.all.filt, dims = 1:10)
  
  nExp <- round(ncol(sce.all.filt) * 0.04) # expect 4% doublets
  sce.all.filt <- doubletFinder_v3(sce.all.filt, 
                                   pN = 0.25, pK = 0.09, 
                                   nExp = nExp, PCs = 1:10)
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(sce.all.filt@meta.data)[grepl("DF.classification", 
                                                   colnames(sce.all.filt@meta.data))]

  #doublet
  phe=sce.all.filt@meta.data
  phe
})



kpCells=unlist(lapply(phe_lt, function(x){
  #table(x[,ncol(x)])
  rownames(x[ x[,ncol(x)]=='Singlet', ])
}))
kp = colnames(seu) %in% kpCells
table(kp)
sce=seu[,kp]



