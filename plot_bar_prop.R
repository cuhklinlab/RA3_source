# ~/anaconda3/bin/R
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(rafalib)
library(Rtsne)
library(umap)
library(hash)
library(gridExtra)

theme_set(theme_gray() +
            theme(
              axis.line = element_line(size=0.5),
              panel.background = element_rect(fill=NA,size=rel(20)),
              panel.grid.minor = element_line(colour = NA),
              axis.text = element_text(size=16),
              axis.title = element_text(size=18)
            )
)

res_table <- read.csv("./input_of_R/clustering_scores_Q2_RA3.csv")
res_table$peak.selection.k <- as.character(res_table$peak.selection.k)
res_table$H1drop_pct <- as.character(res_table$H1drop_pct)

res_table$DR.methods <- factor(res_table$DR.methods, levels=c("scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC","RA3"))



interval = seq(0.1,0.9,0.2)

pdf(file=sprintf("../results/plot/1021/Q2_bar.pdf"), width=10, height=6)

H2type <- "MgOc"

for (itv in interval){
  data <- which(res_table$H1drop_pct==itv & res_table$H2type==H2type & res_table$Clustering.methods%in%c("Louvain clustering","scABC") & res_table$peak.selection.k==3)
  data <- res_table[data,]
  # ari <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=ARI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$ARI,0),1)) + labs(title=H2type) )
  # hom <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=Homogeneity, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$Homogeneity,0),1)) + labs(title=H2type) )
  # nmi <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=NMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$NMI,0),1)) + labs(title=H2type) )
  ami <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=AMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$AMI,0),1)) + labs(title=sprintf('%s_H1dropPct%.2f',H2type,itv)) )

  grid.arrange(ami, nrow = 1)
  # grid.arrange(ari, hom, nmi, ami, nrow = 1)
}


dev.off()


