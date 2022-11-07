
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat


;
df<-seurat.integrated

DefaultAssay(df) <- 'RNA'

data.input <- GetAssayData(df$OLP, assay = "RNA", slot = "data") # normalized data matrix
Idents(df$OLP)<-df$OLP@meta.data$singler_labels_main
labels <- Idents(df$OLP)
Group<-df$OLP$Group
meta <- data.frame(group = Group, row.names = names(labels),labels=labels) # create a dataframe of the cell labels
cellchat<- createCellChat(object = data.input, meta = meta, group.by = "labels")

unique(meta$labels)
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
cellchat <- computeCommunProb(cellchat,nboot = 20, do.fast = TRUE,)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Inference of cell-cell communication network

df.net <- subsetCommunication(cellchat)
table(df.net$pathway_name)

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
table(df.net$LS$pathway_name)

cellchat.OLP<-cellchat


################
################
################


data.input <- GetAssayData(df$Healthy, assay = "RNA", slot = "data") # normalized data matrix
Idents(df$Healthy)<-df$Healthy@meta.data$singler_labels_main
labels <- Idents(df$Healthy)
Group<-df$Healthy$Group
meta <- data.frame(group = Group, row.names = names(labels),labels=labels) # create a dataframe of the cell labels
cellchat<- createCellChat(object = data.input, meta = meta, group.by = "labels")

unique(meta$labels)
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
cellchat <- computeCommunProb(cellchat,nboot = 20, do.fast = TRUE,)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Inference of cell-cell communication network

df.net <- subsetCommunication(cellchat)
table(df.net$pathway_name)

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
table(df.net$LS$pathway_name)

cellchat.Healthy<-cellchat

object.list <- list(Healthy = cellchat.Healthy,OLP = cellchat.OLP)
cellchat  <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

png("./png/cellchat_compareInteractions_main.png",width=2000,height=2500,res=500)
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
} 




# Compare the major sources and targets in 2D space.


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}

for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

png("./png/cellchat.rankNet.png",width=4000,height=4000,res=500)
dev.off()


#
table(df.net$OLP$pathway_name)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4+ T cells",comparison = c(1, 2), signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Monocytes", comparison = c(1, 2),signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3<- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8+ T cells", comparison = c(1, 2),signaling.exclude = c("MIF"))

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

