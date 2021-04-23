library(DESeq2)
library(ggplot2)
library(RColorBrewer)


setwd("~/Documents/MPI/Dynamics/manuscript_multi/Publication_2021/Chlamy_Project_Feb2021/")

Count<-read.csv("TableS2_RawCounts_Average.csv", row.names = 1, sep = ";")
col_data= read.csv("Infotable_Chlamy.csv", sep=";")
count_data=Count[, as.character(col_data$Sample)]


#DESeq norm
dds_counts=DESeqDataSetFromMatrix(countData = count_data, colData =col_data, design = ~ Sample )
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized=TRUE)

rlog_counts <- rlog(dds_counts, blind = TRUE)
rlog.norm.counts <- assay(rlog_counts)

pc <-prcomp(t(rlog.norm.counts))
pc_df <- as.data.frame(pc$x)
#rownames(Infotable) <- Infotable$Library.Name
pc_df$Sample <- col_data$Sample
pc_df$Predation <- col_data$Predation
pc_df$Treatment <- col_data$Treatment

pdf("pca_plot_new.pdf", width = 10, height = 8)
eigs <- pc$sdev^2
eigs[1] / sum(eigs)
eigs[2] / sum(eigs)

P <- ggplot(data = pc_df, mapping = aes(x=PC1, y=PC2, color=Predation, shape=Treatment )) +
  geom_point(size=4) + geom_text(aes(label=Treatment), nudge_y = -5000) + xlab("PC1 (45.29% variance)") + 
  ylab("PC2 (27.22% variance)")
P <- P + scale_color_manual(values = c("deepskyblue","darkorange1", "darkolivegreen"))
P <- P + scale_size_area(max_size=4)
P <- P + scale_x_continuous(limits=c(-70000, 70000))
P <- P + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 22)) + theme_bw()
P
dev.off()

#Differential expression with DESeq2
#Subset accordingly
col_data1= subset(col_data, col_data$Predation %in% c('Predation'))
count_data1=count_data[, as.character(col_data1$Sample)]

col_data2= subset(col_data, col_data$Treatment %in% c('Rotifer', 'Nitrogen'))
count_data2=count_data[, as.character(col_data2$Sample)]

dds_counts=DESeqDataSetFromMatrix(countData = count_data2, colData =col_data2, design = ~ Treatment)
dds_counts=dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts=estimateSizeFactors(dds_counts)
dds_norm=DESeq(dds_counts)
dds_normr=results(dds_norm)
table(dds_normr$padj < 0.05)

write.csv(DifExp3, file="DifExpression.csv")

# Heatmap
data<- read.csv('DifExpression.csv')
data2<- subset(data, data$padj < 0.05)
Down <- subset(data2, data2$log2FoldChange<0)
Up <-subset(data2, data2$log2FoldChange>0)

norm=t(rlog.norm.counts)
count_data1=norm[, as.character(data2$X)]
HM=t(count_data1)

data_distance=as.dist(1-cor(t(HM),method='spearman'))
data_hclust=hclust(data_distance)
AveR=HM[c(data_hclust$order),]

condition_colors <- c(rep("deepskyblue", 3), rep("darkorange1", 3))

pdf("HeatMap_Predation.pdf", width = 10, height = 8)
heatmap.2(as.matrix(AveR), 
          Rowv = NA, Colv = NA,
          dendrogram = 'none', 
          ColSideColors = condition_colors,
          cexRow=0.4, cexCol =1,
          col=rev(brewer.pal(11,"RdBu")),
          scale='row', trace='none',
          labCol=c("Control", "Rotifer", "Nitrogen", "Control", "Rotifer", "Nitrogen"),
          density.info=c("none"),
          margin=c(5,5),
          lhei = c(1,5))
dev.off()


#Table for publishing
Norm<-as.data.frame(rlog.norm.counts)
x<-Norm[,1:3]
Norm$AvePredation<-rowMeans(x)
x<-Norm[,4:6]
Norm$AveNoPredation<-rowMeans(x)
test<-c('MV3', 'MV6')
x<-Norm[test]
Norm$AveNitrogen<-rowMeans(x)
write.csv(Norm, file="TableS3_NormCounts_Average.csv")

