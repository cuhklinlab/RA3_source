suppressPackageStartupMessages({
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots) 
library(RColorBrewer)
library(vegan)
library(R.matlab)
library(JASPAR2016)
library(hash)
})
set.seed(2020)

load("donor_BM0828_peak300.RData")
pvaluematrix=resultPeakSelection$pvalue
Peaks=read.table("donor_BM0828_peak300.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
peaks=as.vector(as.matrix(Peaks))
rownames(pvaluematrix)=peaks

cluster_peak_num = 1000


peakselected = c()
for (i in 1:ncol(pvaluematrix)) {
  peakselected = c(peakselected,head(order(pvaluematrix[,i], decreasing = FALSE), cluster_peak_num))
}
peakselected = unique(sort(peakselected))

chr=c()
p1=c()
p2=c()
for (i in 1:length(Peaks$V1[peakselected])) {
  strings=unlist(strsplit(Peaks$V1[peakselected][i], "_"))
  chr=c(chr,strings[1])
  p1=c(p1,strings[2])
  p2=c(p2,strings[3])
}

peakrange=data.frame(chr,p1,p2)
peakrange$p1=as.numeric(as.character(peakrange$p1))
peakrange$p2=as.numeric(as.character(peakrange$p2))
head(peakrange)
peaks.gr = GRanges(peakrange[,1], IRanges(peakrange[,2], peakrange[,3]))



cluster_assign <- read.csv("Our_donor_BM0828_peak300_dim5_clusters.tsv", sep='\t')
colData <- data.frame(true_label=(cluster_assign$label))

metaData <- data.frame(true_label=(cluster_assign$label), pred_label=cluster_assign$louvain)
metaDataInd = with(metaData, order(pred_label, true_label))
colData <- colData[metaDataInd,]

donor_matrix=readMat("donor_BM0828_peak300.mat")
matrix_use=donor_matrix$count.mat
matrix_select=matrix_use[peakselected,] # matrix_select=matrix_use[topPeaks,] #matrix_select=matrix_use[peakselected,cell_order]
matrix_select = matrix_select[,metaDataInd]
colnames(matrix_select)=1:ncol(matrix_select)


frag_counts = SummarizedExperiment(assays=SimpleList(counts=matrix_select), rowRanges=peaks.gr, colData=colData)
frag_counts = addGCBias(frag_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
motifs <- getJasparMotifs()
motifs.matched = matchMotifs(motifs, frag_counts, genome = BSgenome.Hsapiens.UCSC.hg19)
dev = computeDeviations(object = frag_counts, annotations = motifs.matched)
dev.scores = deviationScores(dev)
variability = computeVariability(dev)



top_motifs = variability$name[head(order(variability$variability, decreasing = TRUE), 50)]
names(top_motifs) = rownames(variability[head(order(variability$variability, decreasing = TRUE), 50), ])
top_devs = dev.scores[which(rownames(dev.scores) %in% names(top_motifs)), ]
rownames(top_devs) = top_motifs[match(rownames(top_devs), names(top_motifs))]
library(gplots) 
library(RColorBrewer)
scalebluered <- colorRampPalette(c("blue", "white", "red"), space = "rgb")(256)
cols = brewer.pal(7, "Dark2")[1:7]
rowcols_louvain = cols[cluster_assign$louvain+1]    #color for louvain clustering #rowcols = cols[InSilicoSCABC$cluster_assignments]
rowcols_louvain = rowcols_louvain[metaDataInd]

legend_cellType_color <- read.csv("cellTypeAll1_U_color.csv", header=FALSE)
legend_cellType_color <- as.matrix(legend_cellType_color)
legend_cellType_color_rgb <- c()
for (i in 1:nrow(legend_cellType_color)){
  legend_cellType_color_rgb <- c( legend_cellType_color_rgb, rgb(legend_cellType_color[i,1],
                                                                 legend_cellType_color[i,2],
                                                                 legend_cellType_color[i,3]) ) 
}
cellTypeAll_U <- read.csv("cellTypeAll1_U.csv", header=FALSE)
cellTypeAll_U <- c(as.matrix(cellTypeAll_U))
color_map <- hash(cellTypeAll_U, legend_cellType_color_rgb)
rowcols_true <- values(color_map, keys=cluster_assign$label )
rowcols_true = rowcols_true[metaDataInd]
label = cellTypeAll_U

library(vegan)
d = vegdist(top_devs, method = "euclidean", na.rm = TRUE)
col.clus = hclust(d, "centroid")



pdf(file=sprintf("chrom_heatmap_%d.pdf",cluster_peak_num))
heatmap.2(t(top_devs), dendrogram='column', Rowv=FALSE, Colv=as.dendrogram(col.clus), trace='none', col = scalebluered, density.info = "none", RowSideColors = rowcols_true, margin = c(10, 1))
legend("bottomleft", legend = label, 
       col = legend_cellType_color_rgb, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, pch = 15)

heatmap.2(t(top_devs), dendrogram='column', Rowv=FALSE, Colv=as.dendrogram(col.clus), trace='none', col = scalebluered, density.info = "none", RowSideColors = rowcols_louvain, margin = c(10, 1))
legend("bottomleft", legend = c(paste0("cluster ", 1:7)), 
       col = cols, border=FALSE, bty="n", y.intersp = 0.7, cex=0.7, pch = 15)

dev.off()