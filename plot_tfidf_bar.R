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

res_table <- read.csv("./input_of_R/clustering_scores_RA3_20201112.csv")
# head(res_table)
res_table$peak.selection.k <- as.character(res_table$peak.selection.k)
data_names = c('ALL_blood','forebrain_half','tissue4_forebrain')
res_table$Dataset <- factor(res_table$Dataset, levels=data_names)

res_table$DR.methods <- factor(res_table$DR.methods, levels=c("PCA","PCAnorm","tfidf"))



pdf(file=sprintf("../results/plot/1021/Q4_bar.pdf"), width=12, height=5)

data <- which(res_table$Dataset%in%data_names & res_table$peak.selection.k==3)
data <- res_table[data,]
# ari <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=ARI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$ARI,0),1)) + labs(title=H2type) )
# hom <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=Homogeneity, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$Homogeneity,0),1)) + labs(title=H2type) )
# nmi <- (ggplot(data) + geom_bar(aes(x=H1drop_pct, y=NMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$NMI,0),1)) + labs(title=H2type) )
ami <- (ggplot(data) + geom_bar(aes(x=Dataset, y=AMI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0),limits=c(min(data$AMI,0),1)) )

grid.arrange(ami, nrow = 1)

dev.off()


