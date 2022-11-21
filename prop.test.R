
?dittoBarPlot
png("./png/dittoBarplot_Celltype2.png",width=4250,height=2000,res=500)
dittoBarPlot(data, "singler_labels", group.by = "Group")
dev.off()

png("./png/dittoBarplot_Celltype1.png",width=4250,height=2000,res=500)
dittoBarPlot(data, "singler_labels_main", group.by = "Group")
dev.off()



t<-dittoBarPlot(data, "singler_labels_main", group.by = "Group")


DimPlot(object = Integrateddata, reduction = 'umap',
        group.by = 'seurat_clusters', pt.size = 1, label=FALSE,
        cols = c('0'='#F68282','4'='#31C53F','8'='#1FA195','6'='#B95FBB','3'='#D4D915',
                 '7'='#28CECA','1'='#ff9a36','2'='#2FF18B','10'='#aeadb3','11'='#faf4cf',
                 '5'='#CCB1F1','9'='#25aff5','12'='#A4DFF2','15'='#4B4BF7','13'='#AC8F14',
                 '14'='#E6C122'))



count=c(0.32085465,0.67914535)
total=c(0.292598487,0.707401513)

data=rbind(count,total)

chisq.test(data, correct = FALSE)


prop.test(c(40891,28313), c(127444,96764),alternative = "two.sided",correct = TRUE)

?prop.test


