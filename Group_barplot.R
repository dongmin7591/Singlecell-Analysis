bp<-read.csv(file = "./InterCellDB/Barplot/6955.csv")


q<-which(bp$sending%in%"Non classical monocytes")
p<-which(bp$sending%in%"Myeloid dendritic cells")

ncMonocyte<-bp[q,]
mDC<-bp[p,]

colnames(ncMonocyte)

# Change the colors manually
p <- ggplot(data=ncMonocyte, aes(x=receiving, y=power, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))+ theme(axis.text.x = element_text(angle = 90))


# Change the colors manually
p <- ggplot(data=ncMonocyte, aes(x=receiving, y=power, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_cowplot()+labs(title="Non classical monocytes (sending cells)", x="signal receiving cells", y = "Power")+
  theme(
    plot.title = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text( size=14))+theme(axis.text.x = element_text(size=12))
# Use custom colors
p + scale_fill_manual(values=c('#999999','#E69F00'))+ theme(axis.text.x = element_text(angle = 90))

png("./InterCellDB/Healthy/GO:0006955/ncMonocyte_Power_barplot.png",width=4000,height=4000,res=500)
dev.off()
