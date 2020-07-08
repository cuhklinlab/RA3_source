library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
library(rafalib)
library(Rtsne)
library(umap)
library(hash)
library(gridExtra)
display.brewer.all(colorblindFriendly=TRUE)
theme_set(theme_gray() +
            theme(
              axis.line = element_line(size=0.5),
              panel.background = element_rect(fill=NA,size=rel(20)),
              panel.grid.minor = element_line(colour = NA),
              axis.text = element_text(size=16),
              axis.title = element_text(size=18)
            )
)

res_table <- read.csv("./input_of_R/res_pbmc_0503.csv")
res_table$peak.selection.k <- as.character(res_table$peak.selection.k)


plot_bar <- function(dataName){
  # dataName='10xpbmc5kPeak'
  print(ggplot(res_table[which(res_table$Dataset==dataName & res_table$Clustering.methods %in% c('Louvain clustering','scABC')),]) + geom_bar(aes(x=peak.selection.k, y=RAGI, fill=DR.methods), stat="identity", position="dodge", width=0.9) + scale_fill_manual(values=colorRampPalette(brewer.pal(8, name="Paired"))(15)) + scale_y_continuous(expand=c(0,0)) + ggtitle(dataName) )
  # + scale_fill_brewer(palette="Paired")
  # grid.arrange(RAGI, percent1, percent2, nrow = 1)
}
pdf(file="../results/plot/0503/bar_ragi_all.pdf", width=9, height=6)
for (dataName in unique(res_table$Dataset)){
plot_bar(dataName)
}
dev.off()




res_table <- read.csv("./input_of_R/res_pbmc_0503.csv")
res_table$peak.selection.k <- as.character(res_table$peak.selection.k)
used_methods <- c("scABC","SCALE","Scasat","cisTopic","Cusanovich2018","SnapATAC",'RA3')
res_table <- res_table[which(res_table$DR.methods %in% used_methods), ]
res_table$DR.methods <- factor(res_table$DR.methods, levels=used_methods)

res_table <- res_table[which(res_table$Dataset %in% c('10xpbmc5kPeak')), ]
res_table$Dataset <- factor(res_table$Dataset, levels=c('10xpbmc5kPeak'))

pdf(file="../results/plot/0503/bar_ragi.pdf", width=12, height=6)

print(ggplot(res_table[which(res_table$peak.selection.k=='3' & res_table$Clustering.methods %in% c('Louvain clustering','scABC')),]) + geom_bar(aes(x=Dataset, y=RAGI, fill=DR.methods), stat="identity", position="dodge", width=0.9)+ scale_fill_brewer(palette="Paired") + scale_y_continuous(expand=c(0,0))  )

dev.off()

